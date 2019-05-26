//
// Created by david on 2019-03-18.
//
#include <string>
#include <iomanip>
#include <spdlog/spdlog.h>
#include <mps_routines/mps_tools/finite/opt.h>
#include <mps_routines/class_superblock.h>
#include <general/class_tic_toc.h>
#include <LBFGS.h>




std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
MPS_Tools::Finite::Opt::find_optimal_excited_state(const class_superblock & superblock, class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace){
    MPS_Tools::log->trace("Finding optimal excited state");
    internals::initialize_timers();
    internals::initialize_params();
    switch (optSpace){
        case OptSpace::DIRECT:
            return internals::direct_optimization(superblock);
        case OptSpace::GUIDED:
            return internals::guided_optimization(superblock,sim_state);
        case OptSpace::FULL:
            return internals::subspace_optimization(superblock, sim_state , optMode, optSpace);
        case OptSpace::PARTIAL:
            return internals::subspace_optimization(superblock, sim_state , optMode, optSpace);

    }

}

namespace MPS_Tools::Finite::Opt::internals{
    std::shared_ptr<class_tic_toc> t_opt = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_eig = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_ham = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_tot = std::make_shared<class_tic_toc>();

    std::shared_ptr<class_tic_toc> t_vH2v = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_vHv  = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_vH2  = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_vH   = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_op   = std::make_shared<class_tic_toc>();



    std::shared_ptr<LBFGSpp::LBFGSParam<double>> params   = std::make_shared<LBFGSpp::LBFGSParam<double>>();
}


void MPS_Tools::Finite::Opt::internals::initialize_timers(){
    t_opt->reset();
    t_eig->reset();
    t_ham->reset();
    t_tot->reset();

    t_vH2v->reset();
    t_vHv ->reset();
    t_vH2 ->reset();
    t_vH  ->reset();
    t_op  ->reset();

    t_opt->set_properties(true,5,"");
    t_eig->set_properties(true,5,"");
    t_ham->set_properties(true,5,"");
    t_tot->set_properties(true,5,"");

    t_vH2v->set_properties(true,5,"");
    t_vHv ->set_properties(true,5,"");
    t_vH2 ->set_properties(true,5,"");
    t_vH  ->set_properties(true,5,"");
    t_op  ->set_properties(true,5,"");

}


void MPS_Tools::Finite::Opt::internals::initialize_params(){
    using namespace LBFGSpp;
    // READ HERE http://pages.mtu.edu/~msgocken/ma5630spring2003/lectures/lines/lines/node3.html
    // I think c1 corresponds to ftol, and c2 corresponds to wolfe
    params->max_iterations = 2000;
    params->max_linesearch = 60; // Default is 20. 5 is really bad, 80 seems better.
    params->m              = 8;     // Default is 6
    params->past           = 1;     // Or perhaps it was this one that helped.
    params->epsilon        = 1e-5;  // Default is 1e-5.
    params->delta          = 1e-8;  // Default is 0. Trying this one instead of ftol.
    params->ftol           = 1e-4;  // Default is 1e-4. this really helped at threshold 1e-8. Perhaps it should be low. Ok..it didn't
    params->wolfe          = 0.9;   // Default is 0.9
    params->min_step       = 1e-40;
    params->max_step       = 1e+40;
    params->linesearch     = LINE_SEARCH_ALGORITHM::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
}







