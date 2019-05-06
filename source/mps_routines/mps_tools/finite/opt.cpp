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
MPS_Tools::Finite::Opt::find_optimal_excited_state(const class_superblock & superblock, double energy_shift, OptMode optMode, OptSpace optSpace){
    MPS_Tools::log->trace("Finding optimal excited state");
    internals::initialize_timers();
    if (optSpace == OptSpace::DIRECT){
        return internals::direct_optimization(superblock);
    }else {
        return internals::subspace_optimization(superblock, energy_shift , optMode, optSpace);
    }
}

namespace MPS_Tools::Finite::Opt::internals{
    std::shared_ptr<class_tic_toc> t_opt = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_eig = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_ham = std::make_shared<class_tic_toc>();
    std::shared_ptr<class_tic_toc> t_tot = std::make_shared<class_tic_toc>();
}


void MPS_Tools::Finite::Opt::internals::initialize_timers(){
    t_opt = std::make_shared<class_tic_toc>();
    t_eig = std::make_shared<class_tic_toc>();
    t_ham = std::make_shared<class_tic_toc>();
    t_tot = std::make_shared<class_tic_toc>();
    t_opt->set_properties(true,5,"");
    t_eig->set_properties(true,5,"");
    t_ham->set_properties(true,5,"");
    t_tot->set_properties(true,5,"");
}










