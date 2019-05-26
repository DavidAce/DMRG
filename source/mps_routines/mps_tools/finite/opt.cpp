//
// Created by david on 2019-03-18.
//
#include <string>
#include <iomanip>
#include <spdlog/spdlog.h>
#include <mps_routines/mps_tools/finite/opt.h>
#include <mps_routines/class_superblock.h>
#include <general/class_tic_toc.h>





std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
MPS_Tools::Finite::Opt::find_optimal_excited_state(const class_superblock & superblock, class_simulation_state & sim_state, OptMode optMode, OptSpace optSpace){
    MPS_Tools::log->trace("Finding optimal excited state");
    internals::initialize_timers();
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

}


void MPS_Tools::Finite::Opt::internals::initialize_timers(){
//    t_opt = std::make_shared<class_tic_toc>();
//    t_eig = std::make_shared<class_tic_toc>();
//    t_ham = std::make_shared<class_tic_toc>();
//    t_tot = std::make_shared<class_tic_toc>();

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










