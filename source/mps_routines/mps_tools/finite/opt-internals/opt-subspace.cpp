//
// Created by david on 2019-03-18.
//
#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include <algorithms/class_simulation_state.h>
#include <general/class_eigsolver.h>
#include <general/arpack_extra/matrix_product_stl.h>
#include <general/arpack_extra/matrix_product_sparse.h>
#include <spdlog/spdlog.h>
#include <mps_routines/mps_tools/finite/opt.h>
#include <mps_routines/class_superblock.h>
#include <LBFGS.h>



std::vector<int> MPS_Tools::Finite::Opt::internals::generate_size_list(const int shape){
    std::vector<int> nev_list;
    if (shape <= 512){
        nev_list.push_back(-1);
        return nev_list;
    }

    int min_nev = 1;
    min_nev =  shape > 1024 ? std::min(8,shape/16) : min_nev;
    int max_nev = std::max(1,shape/16);
    max_nev = std::min(max_nev,256);
    int tmp_nev = min_nev;
    while (tmp_nev <= max_nev){
        nev_list.push_back(tmp_nev);
        tmp_nev *= 4;
    }

    if (shape <= settings::precision::MaxSizeFullDiag or settings::precision::MaxSizeFullDiag <= 0){ // Only do this for small enough matrices
        nev_list.push_back(-1); // "-1" means doing a full diagonalization with lapack instead of arpack.
    }
    return nev_list;
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
MPS_Tools::Finite::Opt::internals::find_subspace(const class_superblock & superblock, class_simulation_state & sim_state, OptMode & optMode, OptSpace &optSpace){
    MPS_Tools::log->trace("Finding subspace");

    using namespace eigutils::eigSetting;
    t_ham->tic();
    Eigen::MatrixXd H_local = superblock.get_H_local_matrix_real();
    t_ham->toc();
    t_eig->tic();

    // You need to copy the data into StlMatrixProduct, because the PartialPivLU will overwrite the data in H_local otherwise.
    StlMatrixProduct<double> hamiltonian_sparse (H_local.data(),H_local.rows(),Form::SYMMETRIC,Side::R, true);
    t_eig->toc();
    Eigen::VectorXd overlaps;
    std::vector<std::tuple<int,double,double,double,double,double,double>> result_log;
    const auto theta = superblock.get_theta();
    Eigen::Map<const Eigen::VectorXcd> theta_old (theta.data(),theta.size());
    int chain_length = superblock.get_length();
    double t_lu         = 0;

    long best_state_idx, worst_state_idx;
    double max_overlap, min_overlap, sq_sum_overlap;
    double subspace_quality;
//    double offset;

    double prec                       = settings::precision::VarConvergenceThreshold;
    double max_overlap_threshold      = 1 - prec; //1.0/std::sqrt(2); //Slightly less than 1/sqrt(2), in case that the choice is between cat states.
    double subspace_quality_threshold = prec;
    double sparcity    = (double)(H_local.array().cwiseAbs() > 1e-12).count()/(double)H_local.size();
    std::stringstream problem_report;
    problem_report
            << "Starting eigensolver \n"
            << std::setprecision(10)
            << "      mode        : "    << optMode << '\n'
            << "      space       : "    << optSpace << '\n'
            << "      position    : "    << superblock.get_position() << '\n'
            << "      chi         : "    << superblock.get_chi() << '\n'
            << "      shape       : "    << theta.size() << " x " << theta.size() << '\n'
            << "      sparcity    : "    << sparcity << '\n' << '\n' << std::flush;
    MPS_Tools::log->debug(problem_report.str());
    class_eigsolver solver;
    std::string reason = "none";
    bool has_solution = false;
    for (auto nev : generate_size_list(theta.size())){
        t_eig->tic();
        if (nev <= 0 and optMode == OptMode::OVERLAP and has_solution){reason = "good enough for overlap mode"; break;}
        if (nev > 0 and  optSpace == OptSpace::PARTIAL){
            hamiltonian_sparse.set_shift(sim_state.energy_now*chain_length);
            hamiltonian_sparse.FactorOP();
            t_lu = hamiltonian_sparse.t_factorOp.get_last_time_interval();
            hamiltonian_sparse.t_factorOp.reset();
            solver.eigs_stl(hamiltonian_sparse,nev,-1, sim_state.energy_now*chain_length,Form::SYMMETRIC,Ritz::LM,Side::R, true,false);
        }
        else if (nev <= 0 or optSpace == OptSpace::FULL){
            optSpace = OptSpace::FULL;
            nev = theta.size();
            solver.eig<Type::REAL, Form::SYMMETRIC>(H_local,true,false);
        }
        t_eig->toc();
//        auto eigvals           = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solver.solution.meta.cols);
//        auto eigvecs           = Eigen::Map<const Eigen::MatrixXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
        Eigen::VectorXd eigvals           = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solver.solution.meta.cols);
        Eigen::MatrixXd eigvecs           = Eigen::Map<const Eigen::MatrixXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
        Eigen::VectorXcd theta_cplx       = Eigen::Map<const Eigen::VectorXcd>(theta.data(),theta.size());

        overlaps         = (theta_old.adjoint() * eigvecs).cwiseAbs();
        max_overlap      = overlaps.maxCoeff(&best_state_idx);
        min_overlap      = overlaps.minCoeff(&worst_state_idx);
        sq_sum_overlap   = overlaps.cwiseAbs2().sum();
        subspace_quality = 1.0 - sq_sum_overlap;
        has_solution     = true;
//        offset           = sim_state.energy_target - eigvals(best_state_idx)/chain_length;
        result_log.emplace_back(nev, max_overlap,min_overlap,sq_sum_overlap,std::log10(subspace_quality),t_eig->get_last_time_interval(),t_lu);
        if(max_overlap    > 1.0 + 1e-10) throw std::runtime_error("max_overlap larger than one : " + std::to_string(max_overlap));
        if(sq_sum_overlap > 1.0 + 1e-10) throw std::runtime_error("eps larger than one : " + std::to_string(sq_sum_overlap));
        if(min_overlap    < 0.0)         throw std::runtime_error("min_overlap smaller than zero: " + std::to_string(min_overlap));


        if(optSpace == OptSpace::FULL)                    {reason = "full diag"; break;}
        if(max_overlap >= max_overlap_threshold )         {reason = "overlap is good"; break;}
        if(subspace_quality < subspace_quality_threshold) {reason = "subspace quality is good"; break;}

    }
    H_local.resize(0,0);
    MPS_Tools::log->debug("Finished eigensolver -- condition: {}",reason);

    auto eigvals           = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvals<Form::SYMMETRIC>().data(),solver.solution.meta.cols);
    auto eigvecs           = Eigen::Map<const Eigen::MatrixXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
    std::stringstream solver_report;
    solver_report << '\n'
                  << std::setw(12) << std::right << "n eigvecs"
                  << std::setw(24) << std::right << "max overlap"
                  << std::setw(24) << std::right << "min overlap"
                  << std::setw(34) << std::right << "eps = Σ_i |<state_i|old>|^2"
                  << std::setw(32) << std::right << "quality = log10(1 - eps)"
                  << std::setw(18) << std::right << "Eig Time[ms]"
                  << std::setw(18) << std::right << "LU  Time[ms]"
                  << '\n';
    for(auto &log : result_log){
        solver_report
                << std::setprecision(16)
                << std::setw(12) << std::right << std::get<0>(log)
                << std::setw(24) << std::right << std::get<1>(log)
                << std::setw(24) << std::right << std::get<2>(log)
                << std::setw(34) << std::right << std::get<3>(log)
                << std::setw(32) << std::right << std::get<4>(log) << std::setprecision(3)
                << std::setw(18) << std::right << std::get<5>(log)*1000
                << std::setw(18) << std::right << std::get<6>(log)*1000
                << '\n';
    }
    solver_report << '\n' << std::flush;
    MPS_Tools::log->debug(solver_report.str());

    if (optMode == OptMode::OVERLAP){
        return std::make_tuple(eigvecs.col(best_state_idx),eigvals.row(best_state_idx));
    }else{
        return std::make_tuple(eigvecs,eigvals);
    }
}


std::tuple<Eigen::Tensor<std::complex<double>,4>, double>
MPS_Tools::Finite::Opt::internals::subspace_optimization(const class_superblock & superblock, class_simulation_state & sim_state, OptMode &optMode, OptSpace &optSpace){
    MPS_Tools::log->trace("Optimizing in SUBSPACE mode");
    using Scalar = std::complex<double>;
    auto [eigvecs,eigvals]  = find_subspace(superblock,sim_state, optMode,optSpace);

    if (optMode == OptMode::OVERLAP and eigvecs.cols() == 1 and eigvals.rows() == 1){
        return std::make_tuple(
                Textra::Matrix_to_Tensor(eigvecs.cast<Scalar>(), superblock.dimensions()),
                eigvals(0)/superblock.get_length()
        );
    }

    using namespace eigutils::eigSetting;
    double chain_length    = superblock.get_length();
    auto theta     = superblock.get_theta();
    auto theta_old = Eigen::Map<const Eigen::VectorXcd> (theta.data(),theta.size());
    Eigen::VectorXcd theta_new;

    double overlap_new  = 0;
    double energy_new,variance_new;
    //Should really use xstart as the projection towards the previous theta, not best overlapping!
    // Note that alpha_i = <theta_old | theta_new_i> is not supposed to be squared! The overlap
    // Between xstart and theta_old should be
    Eigen::VectorXd xstart = (theta_old.adjoint() * eigvecs).normalized().real();
    std::vector<std::tuple<std::string,int,double,double,double,int,int,double>> opt_log;

    {
        t_opt->tic();
        double energy_0   = MPS_Tools::Common::Measure::energy_mpo(superblock,theta);
        double variance_0 = MPS_Tools::Common::Measure::energy_variance_mpo(superblock,theta,energy_0);
        t_opt->toc();
        Eigen::VectorXcd theta_0 = eigvecs * xstart;
        int iter_0 = 0;
        double overlap_0 = (theta_old.adjoint() * theta_0).cwiseAbs().sum();

        opt_log.emplace_back("Start (best overlap)",theta.size(), energy_0/chain_length, std::log10(variance_0/chain_length), overlap_0, iter_0,0, t_opt->get_last_time_interval());
        using namespace LBFGSpp;
        MPS_Tools::Finite::Opt::internals::subspace_functor
                functor (
                superblock,
                eigvecs,
                eigvals);

        // Create solver and function object
        LBFGSpp::LBFGSSolver<double> solver_3(get_lbfgs_params());
        // x will be overwritten to be the best point found
        double fx;
        t_opt->tic();
        MPS_Tools::log->trace("Running LBFGS");
        int niter = solver_3.minimize(functor, xstart, fx);
        int counter = functor.get_count();
        t_opt->toc();
        xstart.normalize();
        theta_new    = eigvecs * xstart;
        energy_new   = functor.get_energy() / chain_length;
        variance_new = functor.get_variance()/chain_length;
        overlap_new = (theta_old.adjoint() * theta_new).cwiseAbs().sum();
//        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_lbfgs.get_last_time_interval(), 0);
        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt->get_last_time_interval());
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
    MPS_Tools::log->debug(report.str());

    return std::make_tuple(
            Textra::Matrix_to_Tensor(theta_new, superblock.dimensions()),
            energy_new
    );

}







