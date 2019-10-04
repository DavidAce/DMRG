//
// Created by david on 2019-03-18.
//
#include <string>
#include <iomanip>
#include <sstream>
#include <spdlog/spdlog.h>
#include <tools/finite/opt.h>
#include <state/class_finite_state.h>
#include <state/class_environment.h>
#include <model/class_model_base.h>
#include <math/class_eigsolver.h>
#include <math/arpack_extra/matrix_product_hamiltonian.h>
#include <simulation/nmspc_settings.h>
#include <glog/logging.h>

Eigen::Tensor<class_finite_state::Scalar,3>
tools::finite::opt::find_excited_state(const class_finite_state &state, const class_simulation_status &sim_status, OptMode optMode, OptSpace optSpace, OptType optType){
    tools::log->trace("Optimizing excited state");
    using namespace opt::internals;
    using namespace Textra;
    static bool googleLogginghasInitialized = false;
    if(not googleLogginghasInitialized){
        google::InitGoogleLogging(tools::log->name().c_str());
        googleLogginghasInitialized = true;
    }

    std::stringstream problem_report;
    auto dims = state.active_dimensions();
    problem_report << fmt::format("Starting optimization: ")
                   << fmt::format("mode [ {} ] "             , optMode.str())
                   << fmt::format("space [ {} ] "            , optSpace.str())
                   << fmt::format("type [ {} ] "             , optType.str())
                   << fmt::format("position [ {} ] "         , state.get_position())
                   << fmt::format("shape [ {} {} {} ] = {} " , dims[0],dims[1],dims[2], state.active_problem_size());
    tools::log->debug(problem_report.str());

    switch (optSpace.option){
        case opt::SPACE::SUBSPACE:    return internals::ceres_subspace_optimization(state,sim_status, optType, optMode);
        case opt::SPACE::DIRECT:      return internals::ceres_direct_optimization(state,sim_status, optType);
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



//
//template<typename Scalar>
//tools::finite::opt::internals::MultiComponents<Scalar>::MultiComponents(const class_finite_state & state){
//    tools::log->trace("Generating multi components");
//    if constexpr (std::is_same<Scalar,double>::value){
//        mpo                          = state.get_multimpo().real();
//        auto & envL_cplx  = state.get_ENVL(state.active_sites.front());
//        auto & envR_cplx  = state.get_ENVR(state.active_sites.back());
//        auto & env2L_cplx = state.get_ENV2L(state.active_sites.front());
//        auto & env2R_cplx = state.get_ENV2R(state.active_sites.back());
//
//        envL  = envL_cplx.block.real();         envR  = envR_cplx.block.real();
//        env2L = env2L_cplx.block.real();        env2R = env2R_cplx.block.real();
//    }
//
//    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
//        mpo                          = state.get_multimpo();
//        auto & envL_cplx  = state.get_ENVL(state.active_sites.front());
//        auto & envR_cplx  = state.get_ENVR(state.active_sites.back());
//        auto & env2L_cplx = state.get_ENV2L(state.active_sites.front());
//        auto & env2R_cplx = state.get_ENV2R(state.active_sites.back());
//        envL  = envL_cplx.block;         envR  = envR_cplx.block;
//        env2L = env2L_cplx.block;        env2R = env2R_cplx.block;
//    }
//
//
//
//    dsizes        = state.active_dimensions();
//    tools::log->trace("Finished building multicomponents");
//}
//
//template struct tools::finite::opt::internals::MultiComponents<double>;
//template struct tools::finite::opt::internals::MultiComponents<std::complex<double>>;


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



std::pair<double,double> tools::finite::opt::internals::windowed_func_grad(double x,double window){
    double func = 0;
    double grad = 0;
    if (std::abs(x) >= window){
        func = x*x - window*window;
        grad = 2*x;
    }
    return std::make_pair(func,grad);
}

