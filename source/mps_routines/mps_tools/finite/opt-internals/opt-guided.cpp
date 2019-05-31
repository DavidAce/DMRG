//
// Created by david on 2019-05-25.
//


#include <spdlog/spdlog.h>
#include <LBFGS.h>
#include <mps_routines/mps_tools/finite/opt.h>
#include <general/class_tic_toc.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <algorithms/class_simulation_state.h>



std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
MPS_Tools::Finite::Opt::internals::guided_optimization(const class_superblock & superblock, class_simulation_state &sim_state){
    MPS_Tools::log->trace("Optimizing in GUIDED mode");
    using Scalar = std::complex<double>;
    auto theta = superblock.get_theta();


    t_opt->tic();
    double chain_length    = superblock.get_length();
    double energy_new,variance_new,overlap_new;
    Eigen::VectorXd theta_new = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();
    Eigen::VectorXd xstart    = theta_new;
    xstart.conservativeResize(theta.size()+1);
    xstart.tail(1).setConstant(1.0);

    std::vector<reports::direct_opt_tuple> opt_log;
    {
        int iter_0 = 0;
        double energy_0   = MPS_Tools::Common::Measure::energy_mpo(superblock,theta);
        double variance_0 = MPS_Tools::Common::Measure::energy_variance_mpo(superblock,theta,energy_0);
        t_opt->toc();
        std::cout << "t_opt: "<< t_opt->get_last_time_interval() << std::endl;
        opt_log.emplace_back("Initial",theta.size(), energy_0/chain_length, std::log10(variance_0/chain_length), 1.0, iter_0 ,0,t_opt->get_last_time_interval());
        std::cout << "t_opt: "<< t_opt->get_last_time_interval() << std::endl;
        t_opt->tic();
        using namespace LBFGSpp;
        MPS_Tools::Finite::Opt::internals::guided_functor functor (superblock, sim_state);

        // Create solver and function object
        LBFGSpp::LBFGSSolver<double> solver(*params);
        // x will be overwritten to be the best point found
        MPS_Tools::log->trace("Running LBFGS");
        double fx;
        int niter = solver.minimize(functor, xstart, fx);
        int counter = functor.get_count();
        t_opt->toc();
        std::cout << "t_opt: "<< t_opt->get_last_time_interval() << std::endl;

        theta_new = xstart.head(theta.size());
        theta_new.normalize();
        auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();

        energy_new   = functor.get_energy() / chain_length;
        variance_new = functor.get_variance()/chain_length;
        overlap_new  = theta_old.dot(theta_new);
        std::cout << "t_opt: "<< t_opt->get_last_time_interval() << std::endl;
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

    return  std::make_tuple(Textra::Matrix_to_Tensor(theta_new.cast<Scalar>(), superblock.dimensions()), energy_new);
}


