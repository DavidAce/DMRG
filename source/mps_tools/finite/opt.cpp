//
// Created by david on 2019-03-18.
//
#include <string>
#include <iomanip>
#include <sstream>
#include <spdlog/spdlog.h>
#include <mps_tools/finite/opt.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_state/class_environment.h>
#include <model/class_hamiltonian_base.h>
#include <general/class_eigsolver.h>
#include <general/arpack_extra/matrix_product_hamiltonian.h>


std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
mpstools::finite::opt::find_optimal_excited_state(const class_finite_chain_state & state, const class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace,OptType optType){
    mpstools::log->trace("Finding optimal excited state");
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
            << "\t shape "      << dims << " = [ " << size << " ]" << std::flush;
    mpstools::log->debug(problem_report.str());



    switch (optSpace){
        case OptSpace::FULL:        return internals::subspace_optimization(state, sim_state , optMode, optSpace, optType);
        case OptSpace::PARTIAL:     return internals::subspace_optimization(state, sim_state , optMode, optSpace, optType);
        case OptSpace::DIRECT:      return internals::direct_optimization  (state, sim_state, optType);
    }
}

std::tuple<Eigen::Tensor<std::complex<double>,4>, double> mpstools::finite::opt::find_optimal_ground_state(const class_finite_chain_state & state, const class_simulation_state & sim_state, std::string ritz){
    return internals::ground_state_optimization(state,sim_state,ritz);

}



void mpstools::finite::opt::internals::reset_timers(){
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
mpstools::finite::opt::internals::superblock_components<Scalar>::superblock_components(const class_superblock &superblock){
    if constexpr (std::is_same<Scalar,double>::value){
        HA_MPO        = superblock.HA->MPO().real();
        HB_MPO        = superblock.HB->MPO().real();
        Lblock        = superblock.Lblock->block.real();
        Rblock        = superblock.Rblock->block.real();
        Lblock2       = superblock.Lblock2->block.real();
        Rblock2       = superblock.Rblock2->block.real();
    }

    if constexpr (std::is_same<Scalar,std::complex<double>>::value){
        HA_MPO        = superblock.HA->MPO();
        HB_MPO        = superblock.HB->MPO();
        Lblock        = superblock.Lblock->block;
        Rblock        = superblock.Rblock->block;
        Lblock2       = superblock.Lblock2->block;
        Rblock2       = superblock.Rblock2->block;
    }
    HAHB          = HA_MPO.contract (HB_MPO, Textra::idx({1},{0}));
    HAHA          = HA_MPO.contract (HA_MPO, Textra::idx({3},{2}));
    HBHB          = HB_MPO.contract (HB_MPO, Textra::idx({3},{2}));
    Lblock2HAHA   = Lblock2.contract(HAHA, Textra::idx({2,3},{0,3})).shuffle(Textra::array6{0,3,2,4,5,1});
    Rblock2HBHB   = Rblock2.contract(HBHB, Textra::idx({2,3},{1,4})).shuffle(Textra::array6{0,3,2,4,5,1});
    HAHB2         = HAHB.contract(HAHB, Textra::idx({2,5},{1,4}));
    dsizes        = superblock.dimensions();

}

template struct mpstools::finite::opt::internals::superblock_components<double>;
template struct mpstools::finite::opt::internals::superblock_components<std::complex<double>>;


double mpstools::finite::opt::internals::windowed_func_abs(double x,double window){
    if (std::abs(x) >= window){
        return std::abs(x)-window;
    }else{
        return 0;
    }
}
double mpstools::finite::opt::internals::windowed_grad_abs(double x,double window){
    if (std::abs(x) >= window){
        return sgn(x);
    }else{
        return 0.0;
    }
}



double mpstools::finite::opt::internals::windowed_func_pow(double x,double window){
    if (std::abs(x) >= window){
        return x*x - window*window;
    }else{
        return 0.0;
    }
}
double mpstools::finite::opt::internals::windowed_grad_pow(double x,double window){
    if (std::abs(x) >= window){
        return 2.0*x;
    }else{
        return 0.0;
    }
}


