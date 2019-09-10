//
// Created by david on 2019-07-15.
//

#ifdef _OPENMP
#include <omp.h>
#define EIGEN_USE_THREADS
#include <unsupported/Eigen/CXX11/Tensor>
Eigen::ThreadPool       tp (Eigen::nbThreads());
Eigen::ThreadPoolDevice dev(&tp,Eigen::nbThreads());
#else
#include <unsupported/Eigen/CXX11/Tensor>
Eigen::DefaultDevice dev;
#endif

#include "ceres_direct_functor.h"
#include <state/class_finite_state.h>



using namespace tools::finite::opt::internals;

template<typename Scalar>
ceres_direct_functor<Scalar>::ceres_direct_functor(
        const class_finite_state & state,
        const class_simulation_status & sim_status)
        : ceres_base_functor(state,sim_status)
{
    tools::log->trace("Constructing subspace functor");

    #ifdef _OPENMP
        tools::log->trace("Parallelizing with {} threads", omp_get_max_threads());
    #endif
    tools::log->trace("Generating multi components");
    energy_reduced      = tools::finite::measure::multisite::energy(state);
    auto state_reduced = tools::finite::measure::reduced::get_state_with_energy_reduced_mpo(state);

    if constexpr (std::is_same<Scalar,double>::value){
        mpo               = state_reduced.get_multimpo().real();
        auto & envL_cplx  = state_reduced.get_ENVL (state_reduced.active_sites.front());
        auto & envR_cplx  = state_reduced.get_ENVR (state_reduced.active_sites.back());
        auto & env2L_cplx = state_reduced.get_ENV2L(state_reduced.active_sites.front());
        auto & env2R_cplx = state_reduced.get_ENV2R(state_reduced.active_sites.back());

        envL  = envL_cplx.block.real();         envR  = envR_cplx.block.real();
        env2L = env2L_cplx.block.real();        env2R = env2R_cplx.block.real();
    }

    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        mpo               = state_reduced.get_multimpo();
        auto & envL_cplx  = state_reduced.get_ENVL (state_reduced.active_sites.front());
        auto & envR_cplx  = state_reduced.get_ENVR (state_reduced.active_sites.back());
        auto & env2L_cplx = state_reduced.get_ENV2L(state_reduced.active_sites.front());
        auto & env2R_cplx = state_reduced.get_ENV2R(state_reduced.active_sites.back());
        envL  = envL_cplx.block;         envR  = envR_cplx.block;
        env2L = env2L_cplx.block;        env2R = env2R_cplx.block;
    }

    dsizes        = state_reduced.active_dimensions();
    tools::log->trace("Finished building multicomponents");


    Hv_tensor.resize(dsizes);
    H2v_tensor.resize(dsizes);
    num_parameters = dsizes[0] * dsizes[1] * dsizes[2];
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){ num_parameters *= 2;}
}










template<typename Scalar>
bool ceres_direct_functor<Scalar>::Evaluate(const double* v_double_double,
                                            double* fx,
                                            double* grad_double_double) const {
    t_op->tic();
    Scalar ene,var;
    Scalar vHv, vH2v;
    double vv,log10var;
    double energy_func, energy_grad;
    double norm_func,norm_grad;
    int vecSize = NumParameters();
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){vecSize = NumParameters()/2;}
    Eigen::Map<const VectorType> v (reinterpret_cast<const Scalar*>(v_double_double)   , vecSize);
    auto lambdas =    Eigen::VectorXd::Ones(1) ;
    vv    = v.squaredNorm();
    norm  = std::sqrt(vv);
    get_H2v(v);
    get_Hv(v);

    auto Hv      = Eigen::Map<VectorType>(Hv_tensor.data() ,Hv_tensor.size());
    auto H2v     = Eigen::Map<VectorType>(H2v_tensor.data(),H2v_tensor.size());
    print_path   = false;
    vHv          = v.dot(Hv);
    vH2v         = v.dot(H2v);


    ene             = vHv/vv;
    var             = vH2v/vv - ene*ene;
//    var             = vH2v/vv;
//    double loss_of_precision = std::log10(std::abs(ene*ene));
//    double expected_error    = std::pow(10, -(13-loss_of_precision));
//    if (std::imag(ene)      > expected_error) tools::log->warn("Energy has imaginary component              : {:.16f} + i {:.16f}" , std::real(ene)    , std::imag(ene));
//    if (std::imag(vH2v/vv)  > expected_error) tools::log->warn("Hamiltonian squared has imaginary component : {:.16f} + i {:.16f}" , std::real(vH2v/vv), std::imag(vH2v/vv));
//    if (std::imag(var)      > expected_error) tools::log->warn("Variance has imaginary component            : {:.16f} + i {:.16f}" , std::real(var)    , std::imag(var));
    if (std::real(var)      < 0.0           ) tools::log->warn("Counter = {}. Variance is negative:  {:.16f} + i {:.16f}" , counter, std::real(var)    , std::imag(var));
    // Make sure var is valid
    var = std::real(var) <= 0.0 ? 1e-20 : std::real(var);

//    energy         = std::real(ene)/length;
    energy         = std::real(ene + energy_reduced) / length;
    variance       = std::abs(var)/length;
    energy_dens    = (energy - energy_min ) / (energy_max - energy_min);
    energy_offset  = energy_dens - energy_target_dens;
    energy_func    = windowed_func_pow(energy_offset,energy_window);
    energy_grad    = windowed_grad_pow(energy_offset,energy_window);


    norm_offset    = std::abs(vv) - 1.0 ;
    std::tie(norm_func,norm_grad) = windowed_func_grad(norm_offset,0.0);
    log10var       = std::log10(variance);
    if(fx != nullptr){
        fx[0] = log10var
                + energy_func * lambdas(0)
                + norm_func;
    }


    if (grad_double_double != nullptr){
        auto vv_1  = std::pow(vv,-1);
        auto var_1 = 1.0/var/std::log(10);  // Possibly abs(var) is required here
        Eigen::Map<VectorType>  grad (reinterpret_cast<      Scalar*>(grad_double_double), vecSize);
        grad = var_1 * vv_1 * (H2v  - v  * vH2v - 2.0 * ene * (Hv - v * ene))
               + lambdas(0) * energy_grad * vv_1 * (Hv - v * ene)
               +  norm_grad * v;
    }



    std::cout   << std::setprecision(12) << std::fixed
                << " Variance: "   << std::setw(18)   << log10var
                << " Energy : "    << std::setw(18)   << energy
                << " Energy t : "  << std::setw(18)   << energy_target
                << " Energy d : "  << std::setw(18)   << energy_dens
                << " Energy td : " << std::setw(18)   << energy_target_dens
                << " Energy o : "  << std::setw(18)   << energy_offset
                << " norm o : "    << std::setw(18)   << norm_offset
                << " lambda 0: "   << std::setw(18)   << lambdas(0)
                << " fx : "        << std::setw(18)   << fx[0]
                << std::endl;


    if(std::isnan(log10var) or std::isinf(log10var)){
        tools::log->warn("log10 variance is invalid");
        tools::log->warn("vv              = {:.16f} + i{:.16f}" , std::real(vv)  , std::imag(vv));
        tools::log->warn("vH2v            = {:.16f} + i{:.16f}" , std::real(vH2v) ,std::imag(vH2v) );
        tools::log->warn("vHv             = {:.16f} + i{:.16f}" , std::real(vHv)  ,std::imag(vHv)  );
        tools::log->warn("var             = {:.16f} + i{:.16f}" , std::real(var)  ,std::imag(var));
        tools::log->warn("ene             = {:.16f} + i{:.16f}" , std::real(ene)  ,std::imag(ene));
        tools::log->warn("log10(var/L)    = {:.16f}" , std::log10(variance/length) );
        tools::log->warn("energy offset   = {:.16f}" , energy_offset );
        tools::log->warn("norm   offset   = {:.16f}" , norm_offset );
        tools::log->warn("lambda 0        = {:.16f}" , lambdas(0));
        throw std::runtime_error("Direct functor failed at counter = " + std::to_string(counter) );
    }



    counter++;
    t_op->toc();
    return true;
}




template<typename Scalar>
void ceres_direct_functor<Scalar>::get_H2v (const VectorType &v)const{
    t_vH2->tic();
    size_t log2chiL  = std::log2(dsizes[1]);
    size_t log2chiR  = std::log2(dsizes[2]);
    size_t log2spin  = std::log2(dsizes[0]);
    Eigen::Tensor<Scalar,3> vH2;
    if (log2spin > log2chiL + log2chiR){
        if (log2chiL > log2chiR){
            if (print_path) tools::log->trace("get_H2v path: log2spin > log2chiL + log2chiR  and  log2chiL > log2chiR ");
            Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), dsizes).shuffle(Textra::array3{1,0,2});
            H2v_tensor.device(dev) =
                    theta
                            .contract(env2L, Textra::idx({0}, {0}))
                            .contract(mpo  , Textra::idx({0,3}, {2,0}))
                            .contract(env2R, Textra::idx({0,3}, {0,2}))
                            .contract(mpo  , Textra::idx({2,1,4}, {2,0,1}))
                            .shuffle(Textra::array3{2,0,1});
        }

        else{
            if (print_path) tools::log->trace("get_H2v path: log2spin > log2chiL + log2chiR  and  log2chiL <= log2chiR ");
            Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), dsizes).shuffle(Textra::array3{2,0,1});
            H2v_tensor.device(dev) =
                    theta
                            .contract(env2R, Textra::idx({0}, {0}))
                            .contract(mpo  , Textra::idx({0,3}, {2,1}))
                            .contract(env2L, Textra::idx({0,3}, {0,2}))
                            .contract(mpo  , Textra::idx({2,4,1}, {2,0,1}))
                            .shuffle(Textra::array3{2,1,0});
        }

    }else{
        if (print_path) tools::log->trace("get_H2v path: log2spin <= log2chiL + log2chiR");
        Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), dsizes).shuffle(Textra::array3{1,0,2});
        H2v_tensor.device(dev) =
                theta
                        .contract(env2L, Textra::idx({0}, {0}))
                        .contract(mpo  , Textra::idx({0,3}, {2,0}))
                        .contract(mpo  , Textra::idx({4,2}, {2,0}))
                        .contract(env2R, Textra::idx({0,2,3}, {0,2,3}))
                        .shuffle(Textra::array3{1, 0, 2});
    }

    t_vH2->toc();
}



template<typename Scalar>
void ceres_direct_functor<Scalar>::get_Hv (const VectorType &v)const{
    t_vH->tic();
    size_t log2chiL  = std::log2(dsizes[1]);
    size_t log2chiR  = std::log2(dsizes[2]);
//            size_t log2spin  = std::log2(multiComponents.dsizes[0]);
    if (log2chiL > log2chiR){
        if (print_path) tools::log->trace("get_Hv path: log2chiL > log2chiR ");

        Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), dsizes).shuffle(Textra::array3{1,0,2});
        Hv_tensor.device(dev) =
                theta
                        .contract(envL, Textra::idx({0}, {0}))
                        .contract(mpo , Textra::idx({0,3}, {2,0}))
                        .contract(envR, Textra::idx({0,2}, {0, 2}))
                        .shuffle(Textra::array3{1, 0, 2});
    }else{
        if (print_path) tools::log->trace("get_Hv path: log2chiL <= log2chiR ");

        Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), dsizes).shuffle(Textra::array3{2,0,1});
        Hv_tensor.device(dev) =
                theta
                        .contract(envR, Textra::idx({0}, {0}))
                        .contract(mpo , Textra::idx({0,3}, {2,1}))
                        .contract(envL, Textra::idx({0,2}, {0,2}))
                        .shuffle(Textra::array3{1, 2, 0});
    }

    t_vH->toc();
}


//template<typename Scalar>
//std::pair<typename ceres_direct_functor<Scalar>::VectorType, Scalar>
//ceres_direct_functor<Scalar>::get_Hv_vHv(const VectorType &v)const{
//    auto Hv = get_Hv(v);
//    t_vHv->tic();
//    auto vHv = v.dot(Hv);
//    t_vHv->toc();
//    return std::make_pair(Hv,vHv);
//}
//
//
//template<typename Scalar>
//std::pair<typename ceres_direct_functor<Scalar>::VectorType, Scalar>
//ceres_direct_functor<Scalar>::get_H2v_vH2v(const VectorType &v)const{
//    auto H2v = get_H2v(v);
//
//    return std::make_pair(H2v,vH2v);
//}



template class tools::finite::opt::internals::ceres_direct_functor<double>;
template class tools::finite::opt::internals::ceres_direct_functor<std::complex<double>>;



