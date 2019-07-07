//
// Created by david on 2019-03-18.
//
#include <string>
#include <iomanip>
#include <sstream>
#include <spdlog/spdlog.h>
#include <state/tools/finite/opt.h>
#include <state/class_finite_state.h>
#include <state/class_environment.h>
#include <model/class_model_base.h>
#include <math/class_eigsolver.h>
#include <math/arpack_extra/matrix_product_hamiltonian.h>
#include <simulation/nmspc_settings.h>

Eigen::Tensor<class_finite_state::Scalar,3>
tools::finite::opt::find_excited_state(const class_finite_state &state, const class_simulation_status &sim_status, OptMode optMode, OptSpace optSpace, OptType optType){
    tools::log->trace("Finding optimal excited state");
    using namespace opt::internals;
    using namespace Textra;
    std::stringstream problem_report;
    auto dims = state.active_dimensions();
    auto size = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());
    problem_report
            << "Starting optimization"
            << std::setprecision(10)
            << "\t mode [ "     << optMode << " ]"
            << "\t space [ "    << optSpace << " ]"
            << "\t type [ "     << optType << " ]"
            << "\t position [ " << state.get_position() << " ]"
            << "\t shape "      << "[ " << dims[0] << " " << dims[1]<< " " << dims[2] << " ] = [ " << size << " ]" << std::flush;
    tools::log->debug(problem_report.str());


    switch (optSpace){
        case OptSpace::SUBSPACE:    return internals::subspace_optimization(state, sim_status, optType, optMode);
        case OptSpace::DIRECT:      return internals::direct_optimization  (state, sim_status, optType);
    }
}

Eigen::Tensor<std::complex<double>,4> tools::finite::opt::find_ground_state(
        const class_finite_state &state, std::string ritz){
    return internals::ground_state_optimization(state,ritz);

}





void tools::finite::opt::internals::reset_timers(){
    t_opt-> reset();
    t_eig-> reset();
    t_ham-> reset();
    t_tot-> reset();
    t_vH2v->reset();
    t_vHv ->reset();
    t_vH2 ->reset();
    t_vH  ->reset();
    t_op  ->reset();
}

template<typename Scalar>
tools::finite::opt::internals::MultiComponents<Scalar>::MultiComponents(const class_finite_state & state){
    tools::log->trace("Generating multi components");
    if constexpr (std::is_same<Scalar,double>::value){
        mpo                          = state.get_multimpo().real();
//        mpo2                         = mpo.contract (mpo, Textra::idx({3},{2}));
//        auto [envL_cplx,envR_cplx]   = state.get_multienv();
//        auto [env2L_cplx,env2R_cplx] = state.get_multienv2();
        auto & envL_cplx  = state.get_ENVL(state.active_sites.front());
        auto & envR_cplx  = state.get_ENVR(state.active_sites.back());
        auto & env2L_cplx = state.get_ENV2L(state.active_sites.front());
        auto & env2R_cplx = state.get_ENV2R(state.active_sites.back());

        envL  = envL_cplx.block.real();         envR  = envR_cplx.block.real();
        env2L = env2L_cplx.block.real();        env2R = env2R_cplx.block.real();
    }

    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        mpo                          = state.get_multimpo();
//        mpo2                         = mpo.contract (mpo, Textra::idx({3},{2}));
//        auto [envL_cplx,envR_cplx]   = state.get_multienv();
//        auto [env2L_cplx,env2R_cplx] = state.get_multienv2();
        auto & envL_cplx  = state.get_ENVL(state.active_sites.front());
        auto & envR_cplx  = state.get_ENVR(state.active_sites.back());
        auto & env2L_cplx = state.get_ENV2L(state.active_sites.front());
        auto & env2R_cplx = state.get_ENV2R(state.active_sites.back());


        envL  = envL_cplx.block;         envR  = envR_cplx.block;
        env2L = env2L_cplx.block;        env2R = env2R_cplx.block;
    }
    dsizes        = state.active_dimensions();
}

template struct tools::finite::opt::internals::MultiComponents<double>;
template struct tools::finite::opt::internals::MultiComponents<std::complex<double>>;


double tools::finite::opt::internals::windowed_func_abs(double x,double window){
    if (std::abs(x) >= window){
        return std::abs(x)-window;
    }else{
        return 0;
    }
}
double tools::finite::opt::internals::windowed_grad_abs(double x,double window){
    if (std::abs(x) >= window){
        return sgn(x);
    }else{
        return 0.0;
    }
}



double tools::finite::opt::internals::windowed_func_pow(double x,double window){
    if (std::abs(x) >= window){
        return x*x - window*window;
    }else{
        return 0.0;
    }
}
double tools::finite::opt::internals::windowed_grad_pow(double x,double window){
    if (std::abs(x) >= window){
        return 2.0*x;
    }else{
        return 0.0;
    }
}


