//
// Created by david on 2019-05-31.
//
#include <mps_tools/finite/opt.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <algorithms/class_simulation_state.h>

template<typename Scalar>
mpstools::finite::opt::internals::subspace_functor<Scalar>::subspace_functor(
        const class_superblock &superblock,
        const class_simulation_state &sim_state,
        const Eigen::MatrixXcd & eigvecs_,
        const Eigen::VectorXd  & eigvals_)
        :
        base_functor(superblock,sim_state),
        eigvecs(eigvecs_),
        eigvals(eigvals_)
{
    if constexpr(std::is_same<Scalar,double>::value){
        H2 = (eigvecs.adjoint().real() * superblock.get_H_local_sq_matrix<Scalar>().template selfadjointView<Eigen::Upper>() * eigvecs.real());
    }
    if constexpr(std::is_same<Scalar,std::complex<double>>::value){
        H2 = (eigvecs.adjoint() * superblock.get_H_local_sq_matrix<Scalar>().template selfadjointView<Eigen::Upper>() * eigvecs);
    }
}


template<typename Scalar>
double mpstools::finite::opt::internals::subspace_functor<Scalar>::operator()(const Eigen::VectorXd &v_double_double, Eigen::VectorXd &grad_double_double) {
    Scalar vH2v,vHv;
    Scalar ene,var;
    double vv,fx,log10var;
    double norm_func,norm_grad;
    VectorType Hv, H2v;
    int vecSize = v_double_double.size();
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){vecSize = v_double_double.size()/2;}
    auto v    = Eigen::Map<const VectorType> (reinterpret_cast<const Scalar*>(v_double_double.data())   , vecSize);
    auto grad = Eigen::Map<      VectorType> (reinterpret_cast<      Scalar*>(grad_double_double.data()), vecSize);
    vv    = v.squaredNorm();
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
//    double loss_of_precision = std::log10(std::abs(ene*ene));
//    double expected_error    = std::pow(10, -(14-loss_of_precision));
//    if (std::imag(ene)      > expected_error) mpstools::log->warn("Energy has imaginary component              : {:.16f} + i {:.16f}" , std::real(ene)    , std::imag(ene));
//    if (std::imag(vH2v/vv)  > expected_error) mpstools::log->warn("Hamiltonian squared has imaginary component : {:.16f} + i {:.16f}" , std::real(vH2v/vv), std::imag(vH2v/vv));
//    if (std::imag(var)      > expected_error) mpstools::log->warn("Variance has imaginary component            : {:.16f} + i {:.16f}" , std::real(var)    , std::imag(var));
    if (std::real(var)      < 0.0           ) mpstools::log->warn("Variance is negative                        : {:.16f} + i {:.16f}" , std::real(var)    , std::imag(var));

    energy         = std::real(ene);
    variance       = std::abs(var);
    variance       = variance < 1e-15  ? 1e-15 : variance;
    norm_offset    = std::abs(vv) - 1.0 ;
    norm_func      = windowed_func_pow(norm_offset,0.1);
    norm_grad      = windowed_grad_pow(norm_offset,0.1);
    log10var       = std::log10(variance);

    fx = log10var
            +  norm_func;

    auto vv_1  = std::pow(vv,-1);
    auto var_1 = 1.0/var/std::log(10);
    grad = var_1 * vv_1 * (H2v  - v  * vH2v - 2.0 * ene * (Hv - v * ene))
           +  norm_grad * v;

//    vecSize = grad.size();
//    if constexpr (std::is_same<Scalar,std::complex<double>>::value){vecSize = 2*grad.size();}
//    grad_double_double  = Eigen::Map<Eigen::VectorXd> (reinterpret_cast<double*> (grad.data()), vecSize);
//

//    std::cout   << std::setprecision(12) << std::fixed
//            << " Variance: "   << std::setw(18)   << std::log10(variance/length)
//            << " Variance: "   << std::setw(18)   << std::log10((vH2v/vv - vHv * vHv/vv/vv)/(double)length)
//            << " Energy : "    << std::setw(18)   << energy/length
//            << " norm : "      << std::setw(18)   << norm
//            << " normsq : "    << std::setw(18)   << vv
//            << " fx : "        << std::setw(18)   << fx
//            << std::endl;


    if(std::isnan(log10var) or std::isinf(log10var)){
        mpstools::log->warn("log10 variance is invalid");
        mpstools::log->warn("vv              = {:.16f} + i{:.16f}" , std::real(vv)  , std::imag(vv));
        mpstools::log->warn("vH2v            = {:.16f} + i{:.16f}" , std::real(vH2v) ,std::imag(vH2v) );
        mpstools::log->warn("vHv             = {:.16f} + i{:.16f}" , std::real(vHv)  ,std::imag(vHv)  );
        mpstools::log->warn("var             = {:.16f} + i{:.16f}" , std::real(var)  ,std::imag(var));
        mpstools::log->warn("ene             = {:.16f} + i{:.16f}" , std::real(ene)  ,std::imag(ene));
        mpstools::log->warn("log10(var/L)    = {:.16f}" , std::log10(variance/length) );
        std::cout << "v: \n " << v << std::endl;
        throw std::runtime_error("Subspace functor failed at counter = " + std::to_string(counter) );
    }

    counter++;
    return fx;
}



template class mpstools::finite::opt::internals::subspace_functor<double>;
template class mpstools::finite::opt::internals::subspace_functor<std::complex<double>>;