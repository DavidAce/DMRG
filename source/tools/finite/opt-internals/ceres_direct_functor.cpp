//
// Created by david on 2019-07-15.
//

#include "ceres_direct_functor.h"
#include <state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <math/nmspc_math.h>

using namespace tools::finite::opt::internal;

template<typename Scalar>
ceres_direct_functor<Scalar>::ceres_direct_functor(
        const class_state_finite & state,
        const class_simulation_status & sim_status)
        : ceres_base_functor(state,sim_status)
{
    tools::log->trace("Constructing direct functor");

    #ifdef _OPENMP
        tools::log->trace("Parallelizing with {} threads", omp.num_threads);
    #endif
    tools::log->trace("Generating multi components");
    energy_reduced  = state.get_energy_reduced();
    if constexpr (std::is_same<Scalar,double>::value){
        mpo               = state.get_multimpo().real();
        auto & envL_cplx  = state.get_ENVL (state.active_sites.front());
        auto & envR_cplx  = state.get_ENVR (state.active_sites.back());
        auto & env2L_cplx = state.get_ENV2L(state.active_sites.front());
        auto & env2R_cplx = state.get_ENV2R(state.active_sites.back());

        envL  = envL_cplx.block.real();         envR  = envR_cplx.block.real();
        env2L = env2L_cplx.block.real();        env2R = env2R_cplx.block.real();
    }

    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        mpo               = state.get_multimpo();
        auto & envL_cplx  = state.get_ENVL (state.active_sites.front());
        auto & envR_cplx  = state.get_ENVR (state.active_sites.back());
        auto & env2L_cplx = state.get_ENV2L(state.active_sites.front());
        auto & env2R_cplx = state.get_ENV2R(state.active_sites.back());
        envL  = envL_cplx.block;         envR  = envR_cplx.block;
        env2L = env2L_cplx.block;        env2R = env2R_cplx.block;
    }

    dsizes        = state.active_dimensions();
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
    tools::common::profile::t_op->tic();
    Scalar ene,ene2,var;
    Scalar vHv, vH2v;
    double vv,log10var;
    double norm_func,norm_grad;
    int vecSize = NumParameters();
    if constexpr (std::is_same<Scalar,std::complex<double>>::value){vecSize = NumParameters()/2;}
    Eigen::Map<const VectorType> v (reinterpret_cast<const Scalar*>(v_double_double)   , vecSize);
    vv    = v.squaredNorm();
    norm  = std::sqrt(vv);

    if(!grad_double_double or grad_double_double == nullptr){
        tools::log->warn("Gradient ptr is null at step {}", counter);
    }
    get_H2v(v);
    get_Hv(v);


    auto Hv      = Eigen::Map<VectorType>(Hv_tensor.data() ,Hv_tensor.size());
    auto H2v     = Eigen::Map<VectorType>(H2v_tensor.data(),H2v_tensor.size());


    print_path   = false;
    vHv          = v.dot(Hv);
    vH2v         = v.dot(H2v);


    // Do this next bit carefully to avoid negative variance when numbers are very small
    ene             = vHv/vv;
    ene2            = vH2v/vv;
    if (std::real(ene2) < 0.0 ) tools::log->debug("Counter = {}. ene2 is negative:  {:.16f} + i {:.16f}" , counter, std::real(ene2) , std::imag(ene2));
    ene2             = std::real(ene2) <  0.0 ? std::abs(ene2)                         : std::real(ene2);
    ene2             = std::real(ene2) == 0.0 ? std::numeric_limits<double>::epsilon() : std::real(ene2);

    var             = ene2 - ene*ene;
//    tools::log->info("log10var/L = {:<24.18f} + i{:<24.18f} | ene = {:<24.18f} + i {:<24.18f}",
//            std::log10(std::abs(std::real(var))/length),  std::log10(std::abs(std::imag(var))/length),
//            std::real(ene), std::imag(ene)
//            );
    if (std::real(var)  < 0.0 ) tools::log->debug("Counter = {}. var  is negative:  {:.16f} + i {:.16f}" , counter, std::real(var)  , std::imag(var));
    var             = std::real(var) <  0.0 ? std::abs(var)                          : std::real(var);
    var             = std::real(var) == 0.0 ? std::numeric_limits<double>::epsilon() : std::real(var);


    energy         = std::real(ene + energy_reduced) / length;
    variance       = std::abs(var)/length;
    norm_offset    = std::abs(vv) - 1.0 ;
    std::tie(norm_func,norm_grad) = windowed_func_grad(norm_offset,0.05);
    log10var       = std::log10(variance);

    if(fx != nullptr){
        fx[0] = log10var + norm_func;
    }

    Eigen::Map<VectorType>  grad (reinterpret_cast<      Scalar*>(grad_double_double), vecSize);
    if (grad_double_double != nullptr){
        auto vv_1  = std::pow(vv,-1);
        auto var_1 = 1.0/var/std::log(10);
        grad = var_1 * vv_1 * (H2v - 2.0*ene*Hv - (ene2 - 2.0*ene*ene)*v);
        if constexpr (std::is_same<Scalar,double>::value){
            grad *= 2.0;
        }
        grad += norm_grad * v;
    }

//    tools::log->trace("log10 var: {:<24.18f} Energy: {:<24.18f} |Grad|: {:<24.18f} |Grad|_inf: {:<24.18f} SqNorm: {:<24.18f} Norm: {:<24.18f} Norm_func: {:<24.18f} |Norm_grad *v|: {:<24.18f} fx: {:<24.18f}",
//                      std::log10(std::abs(var)/length),
//                      std::real(ene + energy_reduced) / length,
//                      grad.norm(),
//                      grad.cwiseAbs().maxCoeff(),
//                      vv,
//                      norm,
//                      norm_func,
//                      (norm_grad * v).norm(),
//                      fx[0]);


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
        throw std::runtime_error("Direct functor failed at counter = " + std::to_string(counter) );
    }



    counter++;
    tools::common::profile::t_op->toc();
    return true;
}




template<typename Scalar>
void ceres_direct_functor<Scalar>::get_H2v (const VectorType &v)const{

    tools::common::profile::t_vH2->tic();
    size_t log2chiL  = std::log2(dsizes[1]);
    size_t log2chiR  = std::log2(dsizes[2]);
    size_t log2spin  = std::log2(dsizes[0]);
    Eigen::Tensor<Scalar,3> vH2;

    if (log2spin > log2chiL + log2chiR){
        if (log2chiL > log2chiR){
            if (print_path) tools::log->trace("get_H2v path: log2spin > log2chiL + log2chiR  and  log2chiL > log2chiR ");
            Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), dsizes).shuffle(Textra::array3{1,0,2});
            H2v_tensor.device(omp.dev) =
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
            H2v_tensor.device(omp.dev) =
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
        H2v_tensor.device(omp.dev) =
                theta
                        .contract(env2L, Textra::idx({0}, {0}))
                        .contract(mpo  , Textra::idx({0,3}, {2,0}))
                        .contract(mpo  , Textra::idx({4,2}, {2,0}))
                        .contract(env2R, Textra::idx({0,2,3}, {0,2,3}))
                        .shuffle(Textra::array3{1, 0, 2});
    }

    tools::common::profile::t_vH2->toc();
}



template<typename Scalar>
void ceres_direct_functor<Scalar>::get_Hv (const VectorType &v)const{
    tools::common::profile::t_vH->tic();
    size_t log2chiL  = std::log2(dsizes[1]);
    size_t log2chiR  = std::log2(dsizes[2]);
//            size_t log2spin  = std::log2(multiComponents.dsizes[0]);
    if (log2chiL > log2chiR){
        if (print_path) tools::log->trace("get_Hv path: log2chiL > log2chiR ");

        Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), dsizes).shuffle(Textra::array3{1,0,2});
        Hv_tensor.device(omp.dev) =
                theta
                        .contract(envL, Textra::idx({0}, {0}))
                        .contract(mpo , Textra::idx({0,3}, {2,0}))
                        .contract(envR, Textra::idx({0,2}, {0, 2}))
                        .shuffle(Textra::array3{1, 0, 2});
    }else{
        if (print_path) tools::log->trace("get_Hv path: log2chiL <= log2chiR ");

        Eigen::Tensor<Scalar,3> theta = Eigen::TensorMap<const Eigen::Tensor<const Scalar,3>>(v.derived().data(), dsizes).shuffle(Textra::array3{2,0,1});
        Hv_tensor.device(omp.dev) =
                theta
                        .contract(envR, Textra::idx({0}, {0}))
                        .contract(mpo , Textra::idx({0,3}, {2,1}))
                        .contract(envL, Textra::idx({0,2}, {0,2}))
                        .shuffle(Textra::array3{1, 2, 0});
    }

    tools::common::profile::t_vH->toc();
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



template class tools::finite::opt::internal::ceres_direct_functor<double>;
template class tools::finite::opt::internal::ceres_direct_functor<std::complex<double>>;



