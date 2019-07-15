//
// Created by david on 2019-05-31.
//
#include <state/tools/finite/opt.h>
#include <state/class_finite_state.h>
#include <state/class_environment.h>
#include <model/class_model_base.h>
#include <simulation/class_simulation_status.h>

template<typename Scalar>
tools::finite::opt::internals::subspace_functor<Scalar>::subspace_functor(
        const class_finite_state & state,
        const class_simulation_status & sim_status,
        const Eigen::MatrixXcd & eigvecs_,
        const Eigen::VectorXd  & eigvals_)
        :
        base_functor(state,sim_status),
        eigvecs(eigvecs_),
        eigvals(eigvals_)
{
    if constexpr(std::is_same<Scalar,double>::value){
        H2 = (eigvecs.adjoint().real() * state.get_multi_hamiltonian2_matrix().real().template selfadjointView<Eigen::Upper>() * eigvecs.real());
    }
    if constexpr(std::is_same<Scalar,std::complex<double>>::value){
        H2 = (eigvecs.adjoint() * state.get_multi_hamiltonian2_matrix().template selfadjointView<Eigen::Upper>() * eigvecs);
    }
    double sparcity = (H2.array().cwiseAbs2() != 0.0).count()/(double)H2.size();
    tools::log->debug("H_local2 nonzeros: {:.8f} %", sparcity*100);
}


template<typename Scalar>
double tools::finite::opt::internals::subspace_functor<Scalar>::operator()(const Eigen::VectorXd &v_double_double, Eigen::VectorXd &grad_double_double) {
    Scalar vH2v,vHv;
    Scalar ene,var;
    double fx,log10var;
    double norm_func,norm_grad;
    VectorType Hv, H2v;
    int vecSize = v_double_double.size();
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){vecSize = v_double_double.size()/2;}
    auto v    = Eigen::Map<const VectorType> (reinterpret_cast<const Scalar*>(v_double_double.data())   , vecSize);
    double vv = v.squaredNorm();

    norm = std::sqrt(vv);

    auto grad = Eigen::Map<      VectorType> (reinterpret_cast<      Scalar*>(grad_double_double.data()), vecSize);
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
//    if (std::imag(ene)      > expected_error) tools::log->warn("Energy has imaginary component              : {:.16f} + i {:.16f}" , std::real(ene)    , std::imag(ene));
//    if (std::imag(vH2v/vv)  > expected_error) tools::log->warn("Hamiltonian squared has imaginary component : {:.16f} + i {:.16f}" , std::real(vH2v/vv), std::imag(vH2v/vv));
//    if (std::imag(var)      > expected_error) tools::log->warn("Variance has imaginary component            : {:.16f} + i {:.16f}" , std::real(var)    , std::imag(var));
    if (std::real(var)      < 0.0           ) tools::log->warn("Variance is negative                        : {:.16f} + i {:.16f}" , std::real(var)    , std::imag(var));

    energy         = std::real(ene);
    variance       = std::abs(var);
    variance       = variance < 1e-15  ? 1e-15 : variance;
    norm_offset    = std::abs(vv) - 1.0 ;
    norm_func      = windowed_func_pow(norm_offset,0.0);
    norm_grad      = windowed_grad_pow(norm_offset,0.0);
    log10var       = std::log10(variance);

    fx = log10var
            +  norm_func;

    auto vv_1  = std::pow(vv,-1);
    auto var_1 = 1.0/var/std::log(10);
    grad = var_1 * vv_1 * (H2v  - v  * vH2v - 2.0 * ene * (Hv - v * ene))
           +  norm_grad * v;
//    grad = var_1 * (H2v  - 2.0 * ene * Hv )
//           +  norm_grad * v;
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
    return fx;
}



template class tools::finite::opt::internals::subspace_functor<double>;
template class tools::finite::opt::internals::subspace_functor<std::complex<double>>;