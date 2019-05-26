//
// Created by david on 2019-05-25.
//


#include <spdlog/spdlog.h>
#include <LBFGS.h>
#include <iomanip>
#include <mps_routines/mps_tools/finite/opt.h>
#include <general/class_tic_toc.h>
#include <mps_routines/class_superblock.h>
//#include <model/class_hamiltonian_base.h>
#include <mps_routines/class_environment.h>
#include <algorithms/class_simulation_state.h>



std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
MPS_Tools::Finite::Opt::internals::guided_optimization(const class_superblock & superblock, class_simulation_state &sim_state){
    MPS_Tools::log->trace("Optimizing in GUIDED mode");
    using Scalar = std::complex<double>;
    class_tic_toc t_lbfgs(true,5,"lbfgs");

    auto theta = superblock.get_theta();

    std::stringstream problem_report;
    problem_report
            << "Starting optimization \n"
            << std::setprecision(10)
            << "      mode        : "    << OptSpace::GUIDED << '\n'
            << "      position    : "    << superblock.get_position() << '\n'
            << "      chi         : "    << superblock.get_chi() << '\n'
            << "      shape       : "    << theta.size() << " x " << theta.size() << '\n' << '\n' << std::flush;
    MPS_Tools::log->debug(problem_report.str());




    t_opt->tic();
    double chain_length    = superblock.get_length();
    double energy_new,variance_new,overlap_new;
    Eigen::VectorXd xstart = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();
    xstart.conservativeResize(theta.size() + 2); // Add lambda1 and lambda2, lagrange multipliers
    xstart.tail(2).setConstant(1.0);
    Eigen::VectorXd theta_new;





    std::vector<std::tuple<std::string,int,double,Scalar,double,int,int,double>> opt_log;
    {
        int iter_0 = 0;
        double energy_0   = MPS_Tools::Common::Measure::energy_mpo(superblock,theta);
        double variance_0 = MPS_Tools::Common::Measure::energy_variance_mpo(superblock,theta,energy_0);
        t_opt->toc();
        opt_log.emplace_back("Start (best overlap)",theta.size(), energy_0/chain_length, std::log10(variance_0/chain_length), 1.0, iter_0 ,0,t_opt->get_last_time_interval());
        t_opt->tic();
        using namespace LBFGSpp;
        MPS_Tools::Finite::Opt::internals::guided_functor functor (superblock, sim_state);

        // Create solver and function object
        LBFGSpp::LBFGSSolver<double> solver(get_lbfgs_params());
        // x will be overwritten to be the best point found
        double fx;
        MPS_Tools::log->trace("Running LBFGS");
        int niter = solver.minimize(functor, xstart, fx);
        int counter = functor.get_count();
        t_opt->toc();
//        xstart.normalize();
        theta_new = Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,1>>(xstart.data(),theta.size()).real().normalized();
        auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();

        energy_new   = functor.get_energy() / chain_length;
        variance_new = functor.get_variance()/chain_length;
        overlap_new  = theta_old.dot(theta_new);
//        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt->get_last_time_interval(), 0);
        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt->get_last_time_interval());
        MPS_Tools::log->trace("Time in function = {:.3f}", t_opt->get_measured_time()*1000);
        MPS_Tools::log->trace("Finished LBFGS");

    }
    std::stringstream report;
    report    << std::setprecision(16) << '\n'
              <<"    "<< std::setw(24) << std::left << "Algorithm"
              <<"    "<< std::setw(8)  << std::left << "size"
              <<"    "<< std::setw(24) << std::left << "energy"
              <<"    "<< std::setw(44) << std::left << "variance"
              <<"    "<< std::setw(24) << std::left << "overlap"
              <<"    "<< std::setw(8)  << std::left << "iter"
              <<"    "<< std::setw(8)  << std::left << "counter"
              <<"    "<< std::setw(20) << std::left << "Elapsed time [ms]"
              <<"    "<< std::setw(20) << std::left << "Time per count [ms]"
              << '\n';
    for(auto &log : opt_log){
        report   << std::setprecision(16)
                 << "    " << std::setw(24) << std::left << std::fixed << std::get<0>(log)
                 << "    " << std::setw(8)  << std::left << std::fixed << std::get<1>(log)
                 << "    " << std::setw(24) << std::left << std::fixed << std::get<2>(log)
                 << "    " << std::setw(44) << std::left << std::fixed << std::get<3>(log)
                 << "    " << std::setw(24) << std::left << std::fixed << std::get<4>(log)
                 << "    " << std::setw(8)  << std::left << std::fixed << std::get<5>(log) << std::setprecision(3)
                 << "    " << std::setw(8)  << std::left << std::fixed << std::get<6>(log) << std::setprecision(3)
                 << "    " << std::setw(20) << std::left << std::fixed << std::get<7>(log)*1000
                 << "    " << std::setw(20) << std::left << std::fixed << std::get<7>(log)*1000 / (double)std::get<6>(log)
                 << '\n';
    }
    report << '\n';
    MPS_Tools::log->debug(report.str());


    std::stringstream lbfgs_report;
    lbfgs_report
            << std::setprecision(3) << '\n'
            << "    " << std::setw(24) << std::left << "LBFGS Time report"
            << "    " << std::setw(12) << std::left << "vH2v  [ms]"
            << "    " << std::setw(12) << std::left << "vHv  [ms]"
            << "    " << std::setw(12) << std::left << "vH2  [ms]"
            << "    " << std::setw(12) << std::left << "vH  [ms]"
            << "    " << std::setw(12) << std::left << "tot  [ms]"
            << "    " << std::setw(12) << std::left << "op  [ms]"
            << '\n';
    lbfgs_report
            << std::setprecision(3)
            << "    " << std::setw(24) << std::left << " "
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * MPS_Tools::Finite::Opt::internals::t_vH2v->get_measured_time()
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * MPS_Tools::Finite::Opt::internals::t_vHv->get_measured_time()
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * MPS_Tools::Finite::Opt::internals::t_vH2->get_measured_time()
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * MPS_Tools::Finite::Opt::internals::t_vH->get_measured_time()
            << "    " << std::setw(12) << std::left << std::fixed << 1000 *
                                                                     (  MPS_Tools::Finite::Opt::internals::t_vH2v->get_measured_time()
                                                                        + MPS_Tools::Finite::Opt::internals::t_vHv->get_measured_time()
                                                                        + MPS_Tools::Finite::Opt::internals::t_vH2->get_measured_time()
                                                                        + MPS_Tools::Finite::Opt::internals::t_vH->get_measured_time()
                                                                     )
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * MPS_Tools::Finite::Opt::internals::t_op->get_measured_time()
            << '\n';

    lbfgs_report << '\n';
    MPS_Tools::log->debug(lbfgs_report.str());

    // Check if energy is outside of the window or not
//    if(energy_new < sim_state.energy_lbound or energy_new > sim_state.energy_ubound) {
//        MPS_Tools::log->debug("DISCARDED");
//        return  std::make_tuple(theta, energy_new);
//    }
//    else{
//    }

    return  std::make_tuple(Textra::Matrix_to_Tensor(theta_new.cast<Scalar>(), superblock.dimensions()), energy_new);


}


