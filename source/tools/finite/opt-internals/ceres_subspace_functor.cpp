//
// Created by david on 2019-07-15.
//

#include "ceres_subspace_functor.h"
#include <state/class_finite_state.h>

using namespace tools::finite::opt::internals;

template<typename Scalar>
tools::finite::opt::internals::ceres_subspace_functor<Scalar>::ceres_subspace_functor(
        const class_finite_state & state,
        const class_simulation_status & sim_status,
        const Eigen::MatrixXcd & eigvecs_,
        const Eigen::VectorXd  & eigvals_)
        :
        ceres_base_functor(state,sim_status),
        eigvecs(eigvecs_),
        eigvals(eigvals_)
{
    tools::log->trace("Constructing subspace functor");
    if constexpr(std::is_same<Scalar,double>::value){
        H2 = state.get_multi_hamiltonian2_subspace_matrix(eigvecs).real();
//        H2 = (eigvecs.adjoint().real() * state.get_multi_hamiltonian2_matrix().real().template selfadjointView<Eigen::Upper>() * eigvecs.real());
    }
    if constexpr(std::is_same<Scalar,std::complex<double>>::value){
        H2 = state.get_multi_hamiltonian2_subspace_matrix(eigvecs);
//        H2 = (eigvecs.adjoint() * state.get_multi_hamiltonian2_matrix().template selfadjointView<Eigen::Upper>() * eigvecs);
    }
    double sparcity = (H2.array().cwiseAbs2() != 0.0).count()/(double)H2.size();
    tools::log->debug("H_local2 nonzeros: {:.8f} %", sparcity*100);
    num_parameters = eigvals.size();
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){num_parameters *= 2;}
}



template<typename Scalar>
bool tools::finite::opt::internals::ceres_subspace_functor<Scalar>::Evaluate(const double* v_double_double,
                                                                             double* fx,
                                                                             double* grad_double_double) const
{
    Scalar vH2v,vHv;
    Scalar ene,var;
    double vv, log10var;
    double norm_func,norm_grad;
    VectorType Hv, H2v;
    int vecSize = NumParameters();
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){vecSize = NumParameters()/2;}
    Eigen::Map<const VectorType> v (reinterpret_cast<const Scalar*>(v_double_double)   , vecSize);
    vv = v.squaredNorm();
    norm = std::sqrt(vv);

#pragma omp parallel
    {
#pragma omp sections
        {
#pragma omp section
            {
                Hv  = eigvals.asDiagonal() * v;
                vHv = v.dot(Hv);
            }
#pragma omp section
            {
                H2v = H2.template selfadjointView<Eigen::Upper>()*v;
                vH2v = v.dot(H2v);
            }
        }
    }
    ene             = vHv/vv;
    var             = vH2v/vv - ene*ene;
    if (std::real(var)      < 0.0           ) tools::log->warn("Variance is negative                        : {:.16f} + i {:.16f}" , std::real(var)    , std::imag(var));

    energy         = std::real(ene) / length;
    variance       = std::abs(var)  / length;
    norm_offset    = std::abs(vv) - 1.0 ;
    std::tie(norm_func,norm_grad) = windowed_func_grad(norm_offset,0.0);

    log10var       = std::log10(variance);
    if (fx != nullptr){
        fx[0] = log10var +  norm_func;
    }

    if (grad_double_double != nullptr){
        auto vv_1  = std::pow(vv,-1);
        auto var_1 = 1.0/var/std::log(10); // Possibly abs(var) is required here
        Eigen::Map<VectorType>  grad (reinterpret_cast<      Scalar*>(grad_double_double), vecSize);
        grad = var_1 * vv_1 * (H2v  - v  * vH2v - 2.0 * ene * (Hv - v * ene))
               +  norm_grad * v;
    }

//    tools::log->trace("Variance: {:<24.18f} Energy: {:<24.18f} norm: {:<24.18f} norm_func: {:<24.18f} norm_grad: {:<24.18f} fx: {:<24.18f} ",
//            std::log10(variance),energy,norm,norm_func,norm_grad,fx[0] );

    if(std::isnan(log10var) or std::isinf(log10var)){
        tools::log->warn("log10 variance is invalid");
        tools::log->warn("vv              = {:.16f} + i{:.16f}" , std::real(vv)  , std::imag(vv));
        tools::log->warn("vH2v            = {:.16f} + i{:.16f}" , std::real(vH2v) ,std::imag(vH2v) );
        tools::log->warn("vHv             = {:.16f} + i{:.16f}" , std::real(vHv)  ,std::imag(vHv)  );
        tools::log->warn("var             = {:.16f} + i{:.16f}" , std::real(var)  ,std::imag(var));
        tools::log->warn("ene             = {:.16f} + i{:.16f}" , std::real(ene)  ,std::imag(ene));
        tools::log->warn("log10(var/L)    = {:.16f}" , std::log10(variance/length) );
        std::cout << "v: \n " << v << std::endl;
        throw std::runtime_error("Subspace functor failed at counter = " + std::to_string(counter) );
    }

    counter++;
    return true;
}



template class tools::finite::opt::internals::ceres_subspace_functor<double>;
template class tools::finite::opt::internals::ceres_subspace_functor<std::complex<double>>;
