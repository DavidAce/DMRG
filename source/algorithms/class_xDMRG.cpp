//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <IO/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_finite_chain.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/nmspc_eigsolver_props.h>
#include <general/class_eigsolver_arpack.h>
#include <general/class_eigsolver.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <LBFGS.h>
#include <algorithms/class_xDMRG_functor.h>
#include <algorithms/class_xDMRG_full_functor.h>

#include "class_xDMRG.h"

using namespace std;
using namespace Textra;

class_xDMRG::class_xDMRG(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_algorithm_base(std::move(hdf5_), "xDMRG",SimulationType::xDMRG) {
    initialize_constants();
    table_xdmrg         = std::make_unique<class_hdf5_table<class_table_dmrg>>(hdf5, sim_name,sim_name);
    table_xdmrg_chain   = std::make_unique<class_hdf5_table<class_table_finite_chain>>(hdf5, sim_name,sim_name + "_chain");
    env_storage         = std::make_shared<class_finite_chain>(num_sites, superblock, hdf5,sim_type,sim_name);
    measurement         = std::make_shared<class_measurement>(superblock, env_storage, sim_type);
    initialize_state(settings::model::initial_state);
    min_saturation_length = 1 * (int)(1.0 * num_sites);
    max_saturation_length = 1 * (int)(2.5 * num_sites);
}



void class_xDMRG::run() {
    if (!settings::xdmrg::on) { return; }
    rn::seed((unsigned long)seed);
    ccout(2) << "\nStarting " << sim_name << " simulation" << std::endl;
    t_tot.tic();
    initialize_chain();
    set_random_fields_in_chain_mpo();
    std::cout << "Parameters on the finite chain: \n" <<std::endl;
    env_storage->print_hamiltonians();
    find_energy_range();
    while(true) {
        single_xDMRG_step();
        env_storage_overwrite_local_ALL();
        store_table_entry_to_file();
        store_chain_entry_to_file();
        store_profiling_to_file_total();
        store_state_to_file();

        check_convergence();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if (iteration >= min_sweeps and env_storage->position_is_the_middle_any_direction() and
           (iteration >= max_sweeps or simulation_has_converged or simulation_has_to_stop))
        {
            break;
        }

        update_bond_dimension(min_saturation_length);
        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        iteration = env_storage->get_sweeps();
        step++;
        ccout(3) << "STATUS: Finished single xDMRG step\n";
    }
    t_tot.toc();
    print_status_full();
    measurement->compute_all_observables_from_finite_chain();
//  Write the wavefunction (this is only defined for short enough chain ( L < 14 say)
    if(settings::xdmrg::store_wavefn){
        hdf5->write_dataset(Textra::to_RowMajor(measurement->mps_chain), sim_name + "/chain/wavefunction");
    }
    print_profiling();
    set_file_OK();
}

void class_xDMRG::single_xDMRG_step() {
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();
    t_opt.tic();
    ccout(3) << "STATUS: Starting single xDMRG step\n";
    Eigen::Tensor<Scalar,4> theta = superblock->MPS->get_theta();
    xDMRG_Mode mode = xDMRG_Mode::KEEP_BEST_OVERLAP;
    mode = chi_temp <= 16 ? mode : xDMRG_Mode::PARTIAL_EIG_OPT;
    mode = chi_temp <= 32 ? mode : xDMRG_Mode::DIRECT_OPT;
//    mode = chi_temp >= 16 and chi_temp <= 64 ? xDMRG_Mode::PARTIAL_EIG_OPT : mode;
//    mode = chi_temp > 64 ? xDMRG_Mode::DIRECT_OPT : mode;
//    mode = theta.size() >= 4096 ? xDMRG_Mode::DIRECT_OPT : mode;

    theta = find_state_with_greatest_overlap_part_diag3(theta,mode);
    t_opt.toc();
    ccout(3) << "STATUS: Truncating theta \n";

    t_svd.tic();
    superblock->truncate_MPS(theta, chi_temp, settings::precision::SVDThreshold);
    t_svd.toc();

    measurement->set_not_measured();
    t_sim.toc();

}


void class_xDMRG::check_convergence(){
//    if(not env_storage->position_is_the_middle()){return;}
//    if(iteration < 5){return;}
    t_sim.tic();
    t_con.tic();
    check_convergence_variance_mpo();
//    if (iteration < min_sweeps){variance_mpo_saturated_for = 0;}
    ccout(2) << "Variance has saturated for " << variance_mpo_saturated_for << " steps \n";
    if(variance_mpo_has_converged)
    {
        ccout(1) << "Simulation has converged\n";
        simulation_has_converged = true;
    }

    else if (variance_mpo_has_saturated and
             bond_dimension_has_reached_max and
             //        env_storage->position_is_the_middle()
             variance_mpo_saturated_for > max_saturation_length
            )
    {
        ccout(1) << "Simulation has to stop\n";
        simulation_has_to_stop = true;
    }


    t_con.toc();
    t_sim.toc();
}

std::vector<int> class_xDMRG::generate_size_list(const int shape){
    std::vector<int> size_list;
    int min_size = 1;
    min_size =  shape > 1024 ? std::min(8,shape/16) : min_size;
    int max_size = std::max(1,shape/16);
//    max_size = std::min(max_size,256);
    int tmp_size = min_size;
    while (tmp_size <= max_size){
        size_list.push_back(tmp_size);
        tmp_size *= 4;
    }

    if (shape <= settings::precision::MaxSizeFullDiag or settings::precision::MaxSizeFullDiag <= 0){ // Only do this for small enough matrices
        size_list.push_back(-1); // "-1" means doing a full diagonalization with lapack instead of arpack.
    }
    return size_list;
}

void class_xDMRG::sort_and_filter_eigenstates(Eigen::VectorXcd &eigvals,
                                 Eigen::MatrixXcd &eigvecs,
                                 Eigen::VectorXd  &overlaps,
                                 int &nev,
                                 double overlap_cutoff)
{
    if (nev == 1){return;}
    auto sorted_idx = make_sorted_index(overlaps);
    int rows = eigvecs.rows();
    int cols = eigvecs.cols();
    Eigen::VectorXcd sorted_eigvals(cols);
    Eigen::MatrixXcd sorted_eigvecs(rows,cols);
    Eigen::VectorXd  sorted_overlaps(cols);

    int i = 0;
    for (auto s : sorted_idx){
        if (overlaps(s) > overlap_cutoff) {
            sorted_eigvals(i)       = eigvals(s);
            sorted_eigvecs.col(i)   = eigvecs.col(s);
            sorted_overlaps(i)      = overlaps(s);
            i++;
        }
    }
    nev      = i;
    sorted_eigvals.conservativeResize(nev);
    sorted_eigvecs.conservativeResize(Eigen::NoChange, nev);
    sorted_overlaps.conservativeResize(nev);
    eigvals  = sorted_eigvals;
    eigvecs  = sorted_eigvecs;
    overlaps = sorted_overlaps;
}


Eigen::Tensor<class_xDMRG::Scalar,4> class_xDMRG::find_state_with_greatest_overlap_part_diag3(Eigen::Tensor<Scalar,4> &theta, xDMRG_Mode mode) {
//    Textra::VectorType <Scalar> theta_new;
    Textra::VectorType <Scalar> theta_res;
    double energy_new   = 0;
    double variance_new = 0;
    double overlap_new  = 0;
    long shape = theta.size();
    double chain_length    = env_storage->get_length();

    if (mode == xDMRG_Mode::DIRECT_OPT){
        theta_res = direct_optimization(energy_new, variance_new, theta);
    }else{

        //    double start_time = t_tot.get_age();
        using namespace eigutils::eigSetting;


        t_ham.tic();
        ccout(3) << "STATUS: Starting construction of H_local \n";
        Eigen::MatrixXd H_local = superblock->get_H_local_matrix_real();
        ccout(3) << "STATUS: Finished construction H_local \n";
        t_ham.toc();
        t_eig.tic();
        ccout(3) << "STATUS: Instantiating StlMatrixProduct \n";
        // You need to copy the data into StlMatrixProduct, because the PartialPivLU will overwrite the data in H_local otherwise.
        StlMatrixProduct<double> hamiltonian_sparse (H_local.data(),H_local.rows(),Form::SYMMETRIC,Side::R, true);
        t_eig.toc();
        Eigen::VectorXd overlaps;
        std::vector<std::tuple<int,double,double,double,double,double,double,double>> result_log;
        Eigen::Map<Textra::VectorType<Scalar>> theta_old (theta.data(),shape);


        double t_lu         = 0;

        long best_state_idx, worst_state_idx;
        double max_overlap, min_overlap, sq_sum_overlap;
        double subspace_quality;
        double offset;

        double prec                       = 1e-10;//settings::precision::VarConvergenceThreshold;
    //    double overlap_cutoff             = prec;
        double max_overlap_threshold      = 1 - prec; //1.0/std::sqrt(2); //Slightly less than 1/sqrt(2), in case that the choice is between cat states.
    //    double max_overlap_threshold      = 1 - 1e-10; //1.0/std::sqrt(2); //Slightly less than 1/sqrt(2), in case that the choice is between cat states.
        double subspace_quality_threshold = prec;
    //    double subspace_quality_threshold = 1e-10;
        double sparcity    = (double)(H_local.array().cwiseAbs() > 1e-15).count()/(double)H_local.size();
        ccout(2) << "Starting eigensolver \n"
                 << std::setprecision(10)
                 << "      position    : "    << env_storage->get_position() << '\n'
                 << "      chi         : "    << chi_temp << '\n'
                 << "      shape       : "    << shape    << " x " << shape << '\n'
                 << "      sparcity    : "    << sparcity << '\n'
                 << "      Wall time   : "    << t_tot.get_age() << '\n' << '\n' << std::flush;

        class_eigsolver solver;
        std::string reason = "none";
        for (auto nev : generate_size_list(shape)){
            t_eig.tic();
            double start_time =  t_tot.get_age();
            if (nev > 0 and (  mode == xDMRG_Mode::KEEP_BEST_OVERLAP or mode == xDMRG_Mode::PARTIAL_EIG_OPT)){
                hamiltonian_sparse.set_shift(energy_now*chain_length);
                hamiltonian_sparse.FactorOP();
                t_lu = hamiltonian_sparse.t_factorOp.get_last_time_interval();
                hamiltonian_sparse.t_factorOp.reset();
                solver.eigs_stl(hamiltonian_sparse,nev,-1, energy_now*chain_length,Form::SYMMETRIC,Ritz::LM,Side::R, true,false);
            }else if (nev <= 0 and mode == xDMRG_Mode::KEEP_BEST_OVERLAP){
                nev = solver.solution.meta.cols;
                reason = "Keeping best overlap";
                break;
            }else{
                nev = shape;
                solver.eig<Type::REAL, Form::SYMMETRIC>(H_local,true,true);
            }
            t_eig.toc();
            auto eigvals           = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solver.solution.meta.cols);
            auto eigvecs           = Eigen::Map<const Eigen::MatrixXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
            overlaps         = (theta_old.adjoint() * eigvecs).cwiseAbs();
    //        sort_and_filter_eigenstates(eigvals, eigvecs, overlaps, nev, overlap_cutoff);
            max_overlap      = overlaps.maxCoeff(&best_state_idx);
            min_overlap      = overlaps.minCoeff(&worst_state_idx);
            sq_sum_overlap   = overlaps.cwiseAbs2().sum();
            subspace_quality = 1.0 - sq_sum_overlap;
            offset           = energy_target - eigvals(best_state_idx)/chain_length;
            result_log.emplace_back(nev, max_overlap,min_overlap,sq_sum_overlap,std::log10(subspace_quality),t_eig.get_last_time_interval(),t_lu,start_time);

            if(max_overlap >= max_overlap_threshold ){reason = "overlap"; break;}
            if(subspace_quality < subspace_quality_threshold){reason = "subspace quality"; break;}

        }
        H_local.resize(0,0);
        ccout(2) << "Finished eigensolver -- condition: " << reason << '\n';

        auto eigvals           = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvals<Form::SYMMETRIC>().data(),solver.solution.meta.cols);
        auto eigvecs           = Eigen::Map<const Eigen::MatrixXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
        int nev = solver.solution.meta.cols;

        ccout(2)  << setw(12) << right << "n eigvecs"
                  << setw(24) << right << "max overlap"
                  << setw(24) << right << "min overlap"
                  << setw(34) << right << "eps = Σ_i |<state_i|old>|^2"
                  << setw(32) << right << "quality = log10(1 - eps)"
                  << setw(18) << right << "Eig Time[ms]"
                  << setw(18) << right << "LU  Time[ms]"
                  << setw(18) << right << "Wall Time [s]"
                  << '\n';
        for(auto &log : result_log){
            ccout(2)  << std::setprecision(16)
                      << setw(12) << right << std::get<0>(log)
                      << setw(24) << right << std::get<1>(log)
                      << setw(24) << right << std::get<2>(log)
                      << setw(34) << right << std::get<3>(log)
                      << setw(32) << right << std::get<4>(log) << std::setprecision(3)
                      << setw(18) << right << std::get<5>(log)*1000
                      << setw(18) << right << std::get<6>(log)*1000
                      << setw(18) << right << std::get<7>(log)
                      << '\n';
        }
        ccout(2) << '\n' << std::flush;

        std::cout << setprecision(16);
        if(nev >= 2 and max_overlap < max_overlap_threshold and mode != xDMRG_Mode::KEEP_BEST_OVERLAP){
            theta_res = subspace_optimization(energy_new,variance_new,nev,eigvecs.data(),eigvals.data(),theta);
        }else {
            theta_res       = eigvecs.col(best_state_idx);
            energy_new      = eigvals(best_state_idx)/chain_length;
            variance_new    = measurement->compute_energy_variance_mpo(theta_res.data(),theta.dimensions(),eigvals(best_state_idx));
        }
    }
    energy_now = energy_new;
//    measurement->set_variance(std::abs(variance_new));

    double energy_ubound = energy_target + 0.1*(energy_max-energy_min);
    double energy_lbound = energy_target - 0.1*(energy_max-energy_min);
    if(energy_now < energy_lbound or energy_now > energy_ubound){
        std::cout << "WARNING: Partial diagonlization -- Energy far from mid-spectrum: " << "Energy =  " << energy_now << " | target = " << energy_target << std::endl;
    }
//    return theta;
    return Textra::Matrix_to_Tensor(theta_res, theta.dimensions());
}

Eigen::Matrix<class_xDMRG::Scalar,Eigen::Dynamic,1> class_xDMRG::subspace_optimization(double &energy_new,double &variance_new,int nev, const double * eigvecs_ptr, const double *eigvals_ptr,const Eigen::Tensor<Scalar,4> &theta){
    using namespace eigutils::eigSetting;
    long shape = theta.size();
    double chain_length    = env_storage->get_length();
    auto eigvals = Eigen::Map<const Eigen::VectorXd> (eigvals_ptr,nev);
    auto eigvecs = Eigen::Map<const Eigen::MatrixXd> (eigvecs_ptr,shape,nev);
    auto theta_old = Eigen::Map<const Textra::VectorType<Scalar>> (theta.data(),shape);
    Textra::VectorType <Scalar> theta_new;
    Textra::VectorType <Scalar> theta_opt;
    Textra::VectorType <Scalar> theta_eigen;
    Textra::VectorType <Scalar> theta_res;

//    double energy_new   = 0;
//    Scalar variance_new = 0;
    double overlap_new  = 0;
    double start_time = t_tot.get_age();
    //Should really use xstart as the projection towards the previous theta, not best overlapping!
    // Note that alpha_i = <theta_old | theta_new_i> is not supposed to be squared! The overlap
    // Between xstart and theta_old should be
    Eigen::VectorXd xstart = (theta_old.adjoint() * eigvecs).normalized().real();

    class_tic_toc t_lbfgs(true,5,"lbfgs");
    std::vector<std::tuple<std::string,size_t,double,Scalar,double,size_t,size_t,double,double>> opt_log;

    {
        ccout(3) << "STATUS: Starting LBFGS \n";
        t_lbfgs.tic();
        ccout(3) << "STATUS: Starting construction of H_local_sq \n";
        Eigen::MatrixXcd H_local_sq = superblock->get_H_local_sq_matrix();
        ccout(3) << "STATUS: Finished construction of H_local_sq \n";
        t_lbfgs.toc();
        Textra::VectorType<Scalar> theta_0 = eigvecs * xstart;
        size_t iter_0 = 0;
        double energy_0 = xstart.cwiseAbs2().cwiseProduct(eigvals.real()).sum() / chain_length;
//            Scalar variance_0 = std::log10(measurement->compute_energy_variance_mpo(theta_0.data(),theta.dimensions(),energy_0));
        Scalar variance_0 = std::log10(
                ((theta_0.adjoint() * H_local_sq.selfadjointView<Eigen::Upper>() * theta_0).sum() -
                 energy_0 * energy_0 * chain_length * chain_length) / chain_length);
        double overlap_0 = (theta_old.adjoint() * theta_0).cwiseAbs().sum();

        opt_log.emplace_back("Start (best overlap)",theta.size(), energy_0, variance_0, overlap_0, iter_0,0, t_lbfgs.get_last_time_interval(), start_time);
        start_time = t_tot.get_age();
        using namespace LBFGSpp;
        class_xDMRG_functor functor
                (shape,nev,
                 H_local_sq,
                 eigvecs_ptr,
                 eigvals_ptr
                );
        double threshold = 1e-10;
        LBFGSpp::LBFGSParam<double> param;
        param.max_iterations = 2000;
        param.m              = 16;
        param.max_linesearch = 200; // Default is 20. 5 is really bad, 80 seems better.
        param.epsilon        = 1e-6;//1e-5; // Default is 1e-5.
        param.delta          = 1e-8; // Default is 0;
        param.ftol           = 1e-5; //this really helped at threshold 1e-8. Perhaps it should be low. Ok..it didn't
        param.past           = 1;   // Or perhaps it was this one that helped... nope


        // Create solver and function object
        LBFGSpp::LBFGSSolver<double> solver_3(param);
        // x will be overwritten to be the best point found
        double fx;
        t_lbfgs.tic();
        ccout(3) << "STATUS: Running LBFGS\n";
        int niter = solver_3.minimize(functor, xstart, fx);
        size_t counter = functor.get_count();
        t_lbfgs.toc();
        xstart.normalize();
        theta_new    = eigvecs * xstart;
        energy_new   = functor.get_energy() / chain_length;
        variance_new = functor.get_variance()/chain_length;
        overlap_new = (theta_old.adjoint() * theta_new).cwiseAbs().sum();
        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_lbfgs.get_last_time_interval(), t_tot.get_age());

        ccout(3) << "STATUS: Finished LBFGS \n";
    }

    ccout(2)  << std::setprecision(16)
              <<"    "<< setw(24) << left << "Algorithm"
              <<"    "<< setw(8)  << left << "size"
              <<"    "<< setw(24) << left << "energy"
              <<"    "<< setw(44) << left << "variance"
              <<"    "<< setw(24) << left << "overlap"
              <<"    "<< setw(8)  << left << "iter"
              <<"    "<< setw(8)  << left << "counter"
              <<"    "<< setw(20) << left << "Elapsed time [ms]"
              <<"    "<< setw(20) << left << "Time per count [ms]"
              <<"    "<< setw(20) << left << "Wall time [s]"
              << '\n';
    for(auto &log : opt_log){
        ccout(2) << std::setprecision(16)
                 << "    " <<setw(24) << left << std::get<0>(log)
                 << "    " <<setw(8)  << left << std::get<1>(log)
                 << "    " <<setw(24) << left << std::get<2>(log)
                 << "    " <<setw(44) << left << std::get<3>(log)
                 << "    " <<setw(24) << left << std::get<4>(log)
                 << "    " <<setw(8)  << left << std::get<5>(log) << std::setprecision(3)
                 << "    " <<setw(8)  << left << std::get<6>(log) << std::setprecision(3)
                 << "    " <<setw(20) << left << std::get<7>(log)*1000
                 << "    " <<setw(20) << left << std::get<7>(log)*1000 / (double)std::get<6>(log)
                 << "    " <<setw(20) << left << std::get<8>(log)
                 << '\n';
    }
    ccout(2) << '\n';

    return theta_new;

}

Eigen::Matrix<class_xDMRG::Scalar,Eigen::Dynamic,1> class_xDMRG::direct_optimization(double &energy_new,
                                                                                     double &variance_new,
                                                                                     const Eigen::Tensor<Scalar, 4> &theta){
    class_tic_toc t_lbfgs(true,5,"lbfgs");
    t_lbfgs.tic();

    long shape = theta.size();
    double chain_length    = env_storage->get_length();

    //Should really use xstart as the projection towards the previous theta, not best overlapping!
    // Note that alpha_i = <theta_old | theta_new_i> is not supposed to be squared! The overlap
    // Between xstart and theta_old should be
    Eigen::VectorXd xstart = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();

    std::vector<std::tuple<std::string,size_t,double,Scalar,double,size_t,size_t,double,double>> opt_log;
    {
        ccout(3) << "STATUS: Starting Direct LBFGS \n";
//        ccout(3) << "STATUS: Starting construction of H_local_sq \n";
//        Eigen::MatrixXcd H_local_sq = superblock->get_H_local_sq_matrix();
//        ccout(3) << "STATUS: Finished construction of H_local_sq \n";
//        Textra::VectorType<Scalar> theta_0 = eigvecs * xstart;
        size_t iter_0 = 0;
        double energy_0   = measurement->compute_energy_mpo(theta.data(),theta.dimensions());
        double variance_0 = std::log10(measurement->compute_energy_variance_mpo(theta.data(),theta.dimensions(),energy_0));

//        Scalar variance_0 = std::log10(
//                ((theta_0.adjoint() * H_local_sq.selfadjointView<Eigen::Upper>() * theta_0).sum() -
//                 energy_0 * energy_0 * chain_length * chain_length) / chain_length);
//        double overlap_0 = (theta_old.adjoint() * theta_0).cwiseAbs().sum();
        t_lbfgs.toc();
        opt_log.emplace_back("Start (best overlap)",theta.size(), energy_0/chain_length, variance_0, 1.0, iter_0 ,0,t_lbfgs.get_last_time_interval(), t_tot.get_age());
        t_lbfgs.tic();
        using namespace LBFGSpp;
        class_xDMRG_full_functor<Scalar> functor(
                superblock->HA->MPO,
                superblock->HB->MPO,
                superblock->Lblock->block,
                superblock->Rblock->block,
                superblock->Lblock2->block,
                superblock->Rblock2->block,
                theta.dimensions()
                );
        double threshold = 1e-5;
        LBFGSpp::LBFGSParam<double> param;
        param.max_iterations = 2000;
        param.max_linesearch = 200; // Default is 20. 5 is really bad, 80 seems better.
        param.m              = 16;
        param.epsilon        = 1e-6; // Default is 1e-5.
        param.delta          = 1e-8; // Default is 0; //Trying this one instead of ftol.
        param.ftol           = 1e-5; //this really helped at threshold 1e-8. Perhaps it should be low. Ok..it didn't
        param.past           = 1; // Or perhaps it was this one that helped.
        param.wolfe          = 1e-1;

        // Create solver and function object
        LBFGSpp::LBFGSSolver<double> solver_3(param);
        // x will be overwritten to be the best point found
        double fx;
        ccout(3) << "STATUS: Running LBFGS\n";
        int niter = solver_3.minimize(functor, xstart, fx);
        size_t counter = functor.get_count();
        t_lbfgs.toc();
        xstart.normalize();
        auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();

        energy_new   = functor.get_energy() / chain_length;
        variance_new = functor.get_variance()/chain_length;
        double overlap_new = (theta_old.adjoint() * xstart).cwiseAbs().sum();
        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_lbfgs.get_last_time_interval(), t_tot.get_age());
        std::cout << "Time in function: " << std::setprecision(3) << functor.t_lbfgs.get_measured_time()*1000 << std::endl;
        ccout(3) << "STATUS: Finished LBFGS \n";
    }

    ccout(2)  << std::setprecision(16)
              <<"    "<< setw(24) << left << "Algorithm"
              <<"    "<< setw(8)  << left << "size"
              <<"    "<< setw(24) << left << "energy"
              <<"    "<< setw(44) << left << "variance"
              <<"    "<< setw(24) << left << "overlap"
              <<"    "<< setw(8)  << left << "iter"
              <<"    "<< setw(8)  << left << "counter"
              <<"    "<< setw(20) << left << "Elapsed time [ms]"
              <<"    "<< setw(20) << left << "Time per count [ms]"
              <<"    "<< setw(20) << left << "Wall time [s]"
              << '\n';
    for(auto &log : opt_log){
        ccout(2) << std::setprecision(16)
                 << "    " <<setw(24) << left << std::get<0>(log)
                 << "    " <<setw(8)  << left << std::get<1>(log)
                 << "    " <<setw(24) << left << std::get<2>(log)
                 << "    " <<setw(44) << left << std::get<3>(log)
                 << "    " <<setw(24) << left << std::get<4>(log)
                 << "    " <<setw(8)  << left << std::get<5>(log) << std::setprecision(3)
                 << "    " <<setw(8)  << left << std::get<6>(log) << std::setprecision(3)
                 << "    " <<setw(20) << left << std::get<7>(log)*1000
                 << "    " <<setw(20) << left << std::get<7>(log)*1000 / (double)std::get<6>(log)
                 << "    " <<setw(20) << left << std::get<8>(log)
                 << '\n';
    }
    ccout(2) << '\n';

    return xstart;
}

void class_xDMRG::initialize_chain() {
    while(true){
        env_storage_insert();
        if (superblock->environment_size + 2ul < (unsigned long) num_sites) {
            enlarge_environment();
            swap();
        } else {
            break;
        }
    }
}

void class_xDMRG::reset_chain_mps_to_random_product_state(std::string parity) {
    std::cout << "Resetting to random product state" << std::endl;
    assert(env_storage->get_length() == num_sites);

    iteration = env_storage->reset_sweeps();

    while(true) {
        // Random product state
        long chiA = superblock->MPS->chiA();
        long chiB = superblock->MPS->chiB();
        long d    = superblock->HA->get_spin_dimension();
        Eigen::Tensor<Scalar,4> theta;
        Eigen::MatrixXcd vecs1(d*chiA,d*chiB);
        Eigen::MatrixXcd vecs2(d*chiA,d*chiB);
        if (parity == "none"){
            theta = Textra::Matrix_to_Tensor(Eigen::MatrixXcd::Random(d*chiA,d*chiB),d,chiA,d,chiB);
        }else{
            if (parity == "sx"){
                vecs1.col(0) = qm::spinOneHalf::sx_eigvecs[0];
                vecs1.col(1) = qm::spinOneHalf::sx_eigvecs[1];
            }else if (parity == "sz"){
                vecs1.col(0) = qm::spinOneHalf::sz_eigvecs[0];
                vecs1.col(1) = qm::spinOneHalf::sz_eigvecs[1];
            }else if (parity == "sy"){
                vecs1.col(0) = qm::spinOneHalf::sy_eigvecs[0];
                vecs1.col(1) = qm::spinOneHalf::sy_eigvecs[1];
            }else{
                std::cerr << "Invalid parity name" << std::endl;
                exit(1);
            }
            theta.resize(d,chiA,d,chiB);
            Eigen::array<long, 4> extent4{2,1,2,1};
            Eigen::array<long, 2> extent2{2,2};
            if(rn::uniform_double_1() < 0.5){
                theta.slice(Eigen::array<long,4>{0,0,0,0},extent4).reshape(extent2) = Textra::Matrix_to_Tensor(vecs1);

            }else{
                theta.slice(Eigen::array<long,4>{0,0,0,0},extent4).reshape(extent2) = Textra::Matrix_to_Tensor(vecs2);
            }

        }
        //Get a properly normalized initial state.
        superblock->truncate_MPS(theta, 1, settings::precision::SVDThreshold);
        env_storage_overwrite_local_ALL();
//        env_storage->print_storage();
        // It's important not to perform the last step.
//        std::cout << "Position :" << env_storage->get_position() << std::endl;
        if(iteration > 1) {break;}
        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        iteration = env_storage->get_sweeps();

    }
    iteration = env_storage->reset_sweeps();
    measurement->set_not_measured();
}

void class_xDMRG::set_random_fields_in_chain_mpo() {
    std::cout << "Setting random fields in chain" << std::endl;
    assert(env_storage->get_length() == num_sites);
    std::vector<std::vector<double>> all_params;
    for (auto &mpo : env_storage->ref_MPO_L()){
        mpo->randomize_hamiltonian();
        all_params.push_back(mpo->get_parameter_values());
    }
    for (auto &mpo : env_storage->ref_MPO_R()){
        mpo->randomize_hamiltonian();
        all_params.push_back(mpo->get_parameter_values());
    }

    for (auto &mpo : env_storage->ref_MPO_L()){
        mpo->set_non_local_parameters(all_params);
    }
    for (auto &mpo : env_storage->ref_MPO_R()){
        mpo->set_non_local_parameters(all_params);
    }


    superblock->HA = env_storage->get_MPO_L().back()->clone();
    superblock->HB = env_storage->get_MPO_R().front()->clone();
    iteration = env_storage->reset_sweeps();
}

void class_xDMRG::find_energy_range() {
    std::cout << "Finding energy range" << std::endl;
    assert(env_storage->get_length() == num_sites);
    int max_sweeps_during_f_range = 5;
    int chi_during_f_range   = 4;
    iteration = env_storage->reset_sweeps();


    // Find energy minimum
    while(true) {
        single_DMRG_step(chi_during_f_range, eigsolver_properties::Ritz::SR);
        env_storage_overwrite_local_ALL();         //Needs to occurr after update_MPS...
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if(iteration >= max_sweeps_during_f_range) {break;}
        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        iteration = env_storage->get_sweeps();
    }
//    env_storage->print_hamiltonians();
//    exit(1);
    compute_observables();
    energy_min = measurement->get_energy_mpo();
    iteration = env_storage->reset_sweeps();

    reset_chain_mps_to_random_product_state();
    // Find energy maximum
    while(true) {
        single_DMRG_step(chi_during_f_range, eigsolver_properties::Ritz::LR);
        env_storage_overwrite_local_ALL();         //Needs to occurr after update_MPS...
        print_status_update();
        // It's important not to perform the last step.
        // That last state would not get optimized
        if(iteration >= max_sweeps_during_f_range) {break;}
        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        iteration = env_storage->get_sweeps();
    }
    compute_observables();
    energy_max = measurement->get_energy_mpo();

    iteration = env_storage->reset_sweeps();
    energy_target = 0.5*(energy_max+energy_min);
    energy_now = superblock->E_optimal / env_storage->get_length();
    std::cout << setprecision(10) ;
    std::cout << "Energy minimum (per site) = " << energy_min << std::endl;
    std::cout << "Energy maximum (per site) = " << energy_max << std::endl;
    std::cout << "Energy target  (per site) = " << energy_target << std::endl;
    double energy_ubound = energy_target + 0.05*(energy_max-energy_min);
    double energy_lbound = energy_target - 0.05*(energy_max-energy_min);
    while(energy_now < energy_lbound or energy_now > energy_ubound){
        reset_chain_mps_to_random_product_state("none");
        compute_observables();
        energy_now = measurement->get_energy_mpo();
    }
    std::cout << "Energy initial (per site) = " << energy_now <<  std::endl;


}

void class_xDMRG::store_table_entry_to_file(bool force){
    if(not force) {
        if (Math::mod(iteration, store_freq) != 0) { return; }
        if (not env_storage->position_is_the_middle_any_direction()) { return; }
        if (store_freq == 0) { return; }
    }
    ccout(3) << "STATUS: Storing table_entry to file\n";
    compute_observables();
    t_sto.tic();
    table_xdmrg->append_record(
            iteration,
            measurement->get_chain_length(),
            env_storage->get_position(),
            measurement->get_chi(),
            chi_max,
            measurement->get_energy_mpo(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            energy_min,
            energy_max,
            energy_target,
            measurement->get_variance_mpo(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            measurement->get_entanglement_entropy(),
            measurement->get_truncation_error(),
            t_tot.get_age());
    t_sto.toc();
}

void class_xDMRG::store_chain_entry_to_file(bool force){
    if (not force) {
        if (Math::mod(iteration, store_freq) != 0) { return; }
        if (store_freq == 0) { return; }
        if (not(env_storage->get_direction() == 1 or env_storage->position_is_the_middle_any_direction())) { return; }
    }
    ccout(3) << "STATUS: Storing chain_entry to file\n";
    t_sto.tic();
    table_xdmrg_chain->append_record(
            iteration,
            measurement->get_chain_length(),
            env_storage->get_position(),
            measurement->get_chi(),
            energy_now,
            measurement->get_entanglement_entropy(),
            measurement->get_truncation_error()
    );
    t_sto.toc();
}

void class_xDMRG::store_state_to_file(bool force){
    if(not force){
        if (Math::mod(iteration, store_freq) != 0) {return;}
        if (not env_storage->position_is_the_middle_any_direction()) {return;}
        if (store_freq == 0){return;}
    }
    ccout(3) << "STATUS: Storing storing mps to file\n";
    t_sto.tic();
    env_storage->write_all_to_hdf5();
    if (settings::hdf5::resume_from_file){
        env_storage->write_full_mps_to_hdf5();
        env_storage->write_full_mpo_to_hdf5();
    }
    t_sto.toc();
}

void class_xDMRG::initialize_constants(){
    using namespace settings;
    num_sites   = xdmrg::num_sites;
    max_sweeps   = xdmrg::max_sweeps;
    chi_max      = xdmrg::chi_max   ;
    chi_grow     = xdmrg::chi_grow  ;
    print_freq   = xdmrg::print_freq;
    store_freq   = xdmrg::store_freq;
    seed         = xdmrg::seed      ;
}

void class_xDMRG::print_profiling(){
    if (settings::profiling::on) {
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_ste.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        print_profiling_sim(t_sim);
        measurement->print_profiling(t_obs);
    }
}

void class_xDMRG::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\n Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_opt.print_time_w_percent(t_parent);
        t_eig.print_time_w_percent(t_parent);
        t_ham.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_con.print_time_w_percent(t_parent);
    }
}