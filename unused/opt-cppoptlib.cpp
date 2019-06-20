//
// Created by david on 2019-06-02.
//

#include <spdlog/spdlog.h>
#include <cppoptlib/solver/lbfgssolver.h>
#include <cppoptlib/meta.h>
#include <mps_tools/finite/opt.h>
#include <general/class_tic_toc.h>
#include <model/class_hamiltonian_base.h>
#include <mps_state/class_environment.h>
#include <mps_state/class_superblock.h>





std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
mpstools::finite::opt::internals::cppoptlib_optimization(const class_superblock & superblock, const class_simulation_state &sim_state){
    mpstools::log->trace("Optimizing in cppoptlib mode");
    using Scalar = std::complex<double>;
    t_opt->tic();
    double chain_length    = superblock.get_length();
    double energy_new,variance_new,overlap_new;

    auto theta = superblock.get_theta();

    std::vector<reports::direct_opt_tuple> opt_log;

    using cScalar = std::complex<double>;
    int iter_0 = 0;
    double energy_0   = mpstools::common::measure::energy_mpo(superblock,theta);
    double variance_0 = mpstools::common::measure::energy_variance_mpo(superblock,theta,energy_0);
    t_opt->toc();
    opt_log.emplace_back("Start (best overlap)",theta.size(), energy_0/chain_length, std::log10(variance_0/chain_length), 1.0, iter_0 ,0,t_opt->get_last_time_interval());
    t_opt->tic();
    using namespace mpstools::finite::opt::internals;
    auto crit = cppoptlib::Criteria<cScalar>::defaults(); // Create a Criteria class to set the solver's stop conditions
    crit.iterations = 400;
    crit.fDelta = 0;
    crit.xDelta = 0;
    crit.gradNorm = 1e-2;
    cppoptlib_functor<cScalar> functor (superblock,sim_state);
    cppoptlib::LbfgsSolver<cppoptlib_functor<cScalar>> solver; // Create solver and function object
    solver.setStopCriteria(crit);
    // x will be overwritten to be the best point found
    mpstools::log->trace("Running LBFGS");
    Eigen::VectorXcd xstart = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size());
    solver.minimize(functor, xstart);
    std::cout << "Solver status: " << solver.status() << std::endl;
    t_opt->toc();
    xstart.normalize();
    auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size());

    energy_new   = functor.get_energy() / chain_length;
    variance_new = functor.get_variance()/chain_length;
    overlap_new  = (theta_old.adjoint() * xstart).cwiseAbs().sum();
    opt_log.emplace_back("cppoptlib",theta.size(), energy_new, std::log10(variance_new), overlap_new, functor.get_iter(),functor.get_count(), t_opt->get_last_time_interval());
    mpstools::log->trace("Finished cppoptlib");

    reports::print_report(opt_log);

    reports::print_report(std::make_tuple(
            mpstools::finite::opt::internals::t_vH2v->get_measured_time(),
            mpstools::finite::opt::internals::t_vHv->get_measured_time(),
            mpstools::finite::opt::internals::t_vH2->get_measured_time(),
            mpstools::finite::opt::internals::t_vH->get_measured_time(),
            mpstools::finite::opt::internals::t_op->get_measured_time()
    ));

    return  std::make_tuple(Textra::Matrix_to_Tensor(xstart, superblock.dimensions()), energy_new);
}



