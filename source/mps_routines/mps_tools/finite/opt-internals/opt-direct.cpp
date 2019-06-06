//
// Created by david on 2019-03-18.
//

#include <spdlog/spdlog.h>
#include <LBFGS.h>
#include <mps_routines/mps_tools/finite/opt.h>
#include <general/class_tic_toc.h>
#include <mps_routines/class_superblock.h>
#include <model/class_hamiltonian_base.h>
#include <mps_routines/class_environment.h>





std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
MPS_Tools::Finite::Opt::internals::direct_optimization(const class_superblock & superblock, const class_simulation_state &sim_state){
    MPS_Tools::log->trace("Optimizing in DIRECT mode");
    using Scalar = std::complex<double>;
    t_opt->tic();
    double chain_length    = superblock.get_length();
    double energy_new,variance_new,overlap_new;

    auto theta = superblock.get_theta();
    Eigen::VectorXd xstart = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();

    std::vector<reports::direct_opt_tuple> opt_log;
    {
        int iter_0 = 0;
        double energy_0   = MPS_Tools::Common::Measure::energy_mpo(superblock,theta);
        double variance_0 = MPS_Tools::Common::Measure::energy_variance_mpo(superblock,theta,energy_0);
        t_opt->toc();
        opt_log.emplace_back("Start (best overlap)",theta.size(), energy_0/chain_length, std::log10(variance_0/chain_length), 1.0, iter_0 ,0,t_opt->get_last_time_interval());
        t_opt->tic();
        using namespace LBFGSpp;
        MPS_Tools::Finite::Opt::internals::direct_functor functor (superblock,sim_state);
        LBFGSpp::LBFGSSolver<double> solver(*params); // Create solver and function object

        // x will be overwritten to be the best point found
        double fx;
        MPS_Tools::log->trace("Running LBFGS");
        int niter = solver.minimize(functor, xstart, fx);
        int counter = functor.get_count();
        t_opt->toc();
        xstart.normalize();
        auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();

        energy_new   = functor.get_energy() / chain_length;
        variance_new = functor.get_variance()/chain_length;
        overlap_new  = (theta_old.adjoint() * xstart).cwiseAbs().sum();
        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt->get_last_time_interval());
        MPS_Tools::log->trace("Finished LBFGS");
    }
    reports::print_report(opt_log);

    reports::print_report(std::make_tuple(
            MPS_Tools::Finite::Opt::internals::t_vH2v->get_measured_time(),
            MPS_Tools::Finite::Opt::internals::t_vHv->get_measured_time(),
            MPS_Tools::Finite::Opt::internals::t_vH2->get_measured_time(),
            MPS_Tools::Finite::Opt::internals::t_vH->get_measured_time(),
            MPS_Tools::Finite::Opt::internals::t_op->get_measured_time()
    ));

    return  std::make_tuple(Textra::Matrix_to_Tensor(xstart.cast<Scalar>(), superblock.dimensions()), energy_new);
}






