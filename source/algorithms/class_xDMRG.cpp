//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <IO/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_finite_chain_sweeper.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/nmspc_eigsolver_props.h>
//#include <general/class_eigsolver_elemental.h>
#include <general/class_eigsolver_armadillo.h>
#include <general/class_eigsolver_arpack.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/class_optimization_cppoptlib.h>
#include <general/class_bfgs_optimization.h>
#include <general/class_lm_optimization.h>
#include <LBFGS.h>
#include <general/class_lbfgspp_optimization.h>
#include <algorithms/class_xDMRG_functor.h>

#include "class_xDMRG.h"

using namespace std;
using namespace Textra;

class_xDMRG::class_xDMRG(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_algorithm_base(std::move(hdf5_), "xDMRG",SimulationType::xDMRG) {
    initialize_constants();
    table_xdmrg         = std::make_unique<class_hdf5_table<class_table_dmrg>>(hdf5, sim_name,sim_name);
    table_xdmrg_chain   = std::make_unique<class_hdf5_table<class_table_finite_chain>>(hdf5, sim_name,sim_name + "_chain");
    env_storage         = std::make_shared<class_finite_chain_sweeper>(max_length, superblock, hdf5,sim_type,sim_name);
    measurement         = std::make_shared<class_measurement>(superblock, env_storage, sim_type);
    initialize_state(settings::model::initial_state);
}



void class_xDMRG::run() {
    if (!settings::xdmrg::on) { return; }
    rn::seed((unsigned long)seed);
    ccout(0) << "\nStarting " << sim_name << " simulation" << std::endl;
    t_tot.tic();
    initialize_chain();
    set_random_fields_in_chain_mpo();
    std::cout << "Parameters on the finite chain: \n" <<std::endl;
    env_storage->print_hamiltonians();
    int full_iters = 5;
    find_energy_range();
//    hdf5->write_dataset(Textra::to_RowMajor(measurement->mps_chain), sim_name + "/chain/wavefunction_initial");

    while(true) {
        single_xDMRG_step();
        env_storage_overwrite_local_ALL();
        store_table_entry_to_file();
        store_chain_entry_to_file();
        store_profiling_to_file_total();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if(env_storage->position_is_the_middle()) {
            check_convergence_all();
            if (iteration >= full_iters
            and (iteration >= max_sweeps
                 or simulation_has_converged
                 or simulation_has_to_stop))
            {
                break;
            }
        }

        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        iteration = env_storage->get_sweeps();
    }
    t_tot.toc();
    print_status_full();
    env_storage->write_chain_to_file();
    measurement->compute_all_observables_from_finite_chain();
    // Write the wavefunction (this is only defined for short enough chain ( L < 14 say)
    if(settings::xdmrg::store_wavefn){
        hdf5->write_dataset(Textra::to_RowMajor(measurement->mps_chain), sim_name + "/chain/wavefunction");
    }
    store_profiling_to_file_total(true);
    print_profiling();

}

void class_xDMRG::single_xDMRG_step() {
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();
    t_opt.tic();
//    bool compare = false;
    Eigen::Tensor<Scalar,4> theta = superblock->MPS->get_theta();
//    xDMRG_Mode mode = chi_temp <= 12 ? xDMRG_Mode::FULL : xDMRG_Mode::PARTIAL;
    xDMRG_Mode mode = xDMRG_Mode::PARTIAL;

    switch (mode){
        case xDMRG_Mode::FULL:
            theta = find_state_with_greatest_overlap_full_diag(theta);
            t_opt.toc();
            break;
        case xDMRG_Mode::PARTIAL:
            theta = find_state_with_greatest_overlap_part_diag2(theta);
            t_opt.toc();
            break;
    }
    t_svd.tic();
    superblock->truncate_MPS(theta, chi_temp, settings::precision::SVDThreshold);
    t_svd.toc();

    measurement->set_not_measured();
    t_sim.toc();




//            if (false) {
//                std::unique_ptr<class_mps_2site> MPS_backup = std::make_unique<class_mps_2site>(*superblock->MPS);
////                superblock->truncate_MPS(theta, 256, settings::precision::SVDThreshold);
//                measurement->set_not_measured();
//                measurement->compute_all_observables_from_superblock(theta);
//                compare = true;
//                std::cout << "Iteration: " << iteration << "    position: " << env_storage->get_position() << std::endl;
//                std::cout << "      Variance before: " << setw(12) << std::setprecision(10)
//                          << std::log10(measurement->get_variance_mpo()) << std::endl;
//                if (iteration == 5 and env_storage->get_position() == 7) {
//                    std::cout << "Singular values before : \n" << superblock->MPS->LC << std::endl;
//                }
//                superblock->MPS->set_mps(
//                        MPS_backup->MPS_A->get_L(),
//                        MPS_backup->MPS_A->get_G(),
//                        MPS_backup->LC,
//                        MPS_backup->MPS_B->get_G(),
//                        MPS_backup->MPS_B->get_L()
//                        );
//
//            }
////            measurement->set_not_measured();
////            measurement->compute_all_observables_from_superblock();
//            if (compare) {
//                std::cout << "      Variance after : " << setw(12) << std::setprecision(10)
//                          << std::log10(measurement->get_variance_mpo()) << std::endl;
//                if (iteration == 5 and env_storage->get_position() == 7) {
//                    std::cout << "Singular values after  : \n" << superblock->MPS->LC << std::endl;
//                }
////                std::cout << "HA = \n"
////                          << superblock->HA->get_parameter_values()[0] << ' '
////                          << superblock->HA->get_parameter_values()[1] << ' '
////                          << superblock->HA->get_parameter_values()[2] << ' '
////                          << std::endl;
////
////
////                std::cout << "HB = \n"
////                          << superblock->HB->get_parameter_values()[0] << ' '
////                          << superblock->HB->get_parameter_values()[1] << ' '
////                          << superblock->HB->get_parameter_values()[2] << ' '
////                          << std::endl;
//
//            }
//
}

Eigen::Tensor<class_xDMRG::Scalar,4> class_xDMRG::find_state_with_greatest_overlap_full_diag(Eigen::Tensor<Scalar, 4> &theta){
//    Eigen::Tensor<Scalar,1> theta = superblock->MPS->get_theta().reshape(superblock->shape1);
    long shape = theta.size();
    double L   = env_storage->get_length();
    t_ham.tic();
    auto H_local = superblock->get_H_local_rank2();
    assert(H_local.dimension(0) == shape);
    assert(H_local.dimension(1) == shape);
    t_ham.toc();
    t_ham.toc();
    //Sparcity seems to be about 5-15%. Better to do dense.
    class_eigsolver_armadillo solver;
    Eigen::VectorXd  eigvals(shape);
    Eigen::MatrixXcd eigvecs(shape,shape);
    t_eig.tic();
    solver.eig_sym(shape, H_local.data(), eigvals.data(), eigvecs.data());
    t_eig.toc();
    Eigen::Map<Textra::VectorType<Scalar>> theta_vector (theta.data(),shape);
    Textra::VectorType<double> overlaps = (theta_vector.conjugate().transpose() * eigvecs).cwiseAbs();
    long best_state;
    double max_overlap = overlaps.cwiseAbs().maxCoeff(&best_state);
    double offset = energy_target - eigvals(best_state)/L;
    if (false){
//    if (std::abs(offset) > 0.01 * (energy_max - energy_min)){
        (eigvals.array() - energy_target*L).cwiseAbs().minCoeff(&best_state);
        max_overlap = overlaps(best_state);
        std::cout << "Finding state with closest energy:  " << std::setprecision(10) << max_overlap << " offset: " << setw(12) << offset << " (target = " << energy_target << ")" << " position: " << env_storage->get_position()<< std::endl;
    }else{
//        std::cout << "Finding state with closest overlap: " << std::setprecision(10) << max_overlap << " offset: " << setw(12) << offset << " (target = " << energy_target << ")" << " position: " << env_storage->get_position()<< std::endl;
    }
//    std::cerr << "Armadillo Found overlap: " << setprecision(14) << setw(18) << max_overlap << std::endl;

    Textra::MatrixType<Scalar> state = eigvecs.col(best_state);
    energy_now = eigvals(best_state)/L;
    return Textra::Matrix_to_Tensor(state, theta.dimensions());
}

Eigen::Tensor<class_xDMRG::Scalar,4> class_xDMRG::find_state_with_greatest_overlap_part_diag(Eigen::Tensor<Scalar,4> &theta) {
    long shape = theta.size();
    double L    = env_storage->get_length();
    t_ham.tic();
    auto H_local = superblock->get_H_local_rank2();
    assert(H_local.dimension(0) == shape);
    assert(H_local.dimension(1) == shape);
    t_ham.toc();

//    long chiA = superblock->MPS->chiA();
//    long chiB = superblock->MPS->chiB();
//    long d    = superblock->HA->get_spin_dimension();
//    std::array<long,4> shape4_theta = {d,chiA,d,chiB};
//    std::array<long,4> shape4_mpo   = superblock->HA->MPO.dimensions();;

    class_eigsolver_arpack<std::complex<double>, Form::GENERAL> arpack_solver;
    Eigen::VectorXcd eigvals;
    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd overlaps;
    std::vector<std::tuple<int,double,double>> results;
    Eigen::Map<Textra::VectorType<Scalar>> theta_vector (theta.data(),shape);
    Textra::VectorType <Scalar> theta_res = theta_vector;
    int nev = std::min((int)(shape/4), 1);
    int ncv = std::min((int)(shape/2), nev*2);
    long best_state_idx;
    double max_overlap;
    double offset;
    double overlap_threshold = 1.0/std::sqrt(2); //Slightly less than 1/sqrt(2), in case that the choice is between cat states.
    for (int iter = 0; iter  <  4; iter++){
        t_eig.tic();
//        arpack_solver.eig_shift_invert(H_local_ptr.data(), shape, nev,ncv,energy_now*L, Ritz::LM, true,true, theta_res.data());
        arpack_solver.eig_shift_invert(H_local.data(), shape, nev,ncv,energy_now*L, Ritz::LM, true, false, theta_res.data());
        t_eig.toc();
        eigvals         = Eigen::Map<const Eigen::VectorXcd> (arpack_solver.get_eigvals().data(),arpack_solver.GetNevFound());
        eigvecs         = Eigen::Map<const Eigen::MatrixXcd> (arpack_solver.get_eigvecs().data(),arpack_solver.Rows(),arpack_solver.Cols());
        overlaps        = (theta_vector.conjugate().transpose() * eigvecs).cwiseAbs();
        max_overlap     = overlaps.maxCoeff(&best_state_idx);
        offset          = energy_target - eigvals(best_state_idx).real()/L;
        results.emplace_back(nev, max_overlap,t_tot.get_age());

        if(max_overlap >= overlap_threshold ){break;}
//        overlap_threshold  *= 0.8;
        nev = std::min((int)(shape/4), nev*4);
        ncv = std::min((int)(shape/2), nev*2);
        if(nev >= shape/4 ){break;}


    }

    if(max_overlap <  overlap_threshold){
        std::cerr << "WARNING: Partial diagonlization -- Overlap very small: " << std::setprecision(10) << max_overlap << "      position: " << env_storage->get_position() << std::endl;
        std::cerr << "          Overlaps:"  << std::endl;
        for (auto &res : results){
            std::cerr << "          nev: " << setw(10) << std::get<0>(res) << "   overlap: " << setprecision(14) << setw(18) << std::get<1>(res) << "       Wall time: " << std::get<2>(res) << std::endl;
        }
        std::cerr << "          Running full diagonalization. Wall time: " << t_tot.get_age() << std::endl;
        class_eigsolver_armadillo solver;
        Eigen::VectorXd eigvals_real(shape);
        eigvecs.resize(shape,shape);
        t_eig.tic();
        solver.eig_sym(shape, H_local.data(), eigvals_real.data(), eigvecs.data());
        eigvals = eigvals_real;
        t_eig.toc();
        overlaps        = (theta_vector.conjugate().transpose() * eigvecs).cwiseAbs();
        max_overlap     = overlaps.maxCoeff(&best_state_idx);
        offset          = energy_target - eigvals(best_state_idx).real()/L;
        results.emplace_back(std::make_tuple(nev, max_overlap,t_tot.get_age()));
        std::cerr << "          Full diag complete.           Wall time: " << t_tot.get_age() << std::endl;
        std::cerr << "          Found overlap: " << setprecision(14) << setw(18) << max_overlap << std::endl;
    }
//    std::cerr << "Arpack Found overlap: " << setprecision(14) << setw(18) << max_overlap << std::endl;

    energy_now      = eigvals(best_state_idx).real()/L;
    theta_res       = eigvecs.col(best_state_idx);

    double energy_ubound = energy_target + 0.1*(energy_max-energy_min);
    double energy_lbound = energy_target - 0.1*(energy_max-energy_min);
    if(energy_now < energy_lbound or energy_now > energy_ubound){
        std::cerr << "WARNING: Partial diagonlization -- Energy far from mid-spectrum: " << "Energy =  " << energy_now << " | target = " << energy_target << std::endl;
    }


//    energy_target = energy_now;
//    if (env_storage->position_is_the_middle()){energy_target = energy_now;}
    return Textra::Matrix_to_Tensor(theta_res, theta.dimensions());
}

std::vector<int> generate_size_list(const int shape){
    std::vector<int> size_list;
    int min_size = 1;
    int max_size = std::max(1,shape/32);
    int tmp_size = min_size;
    while (tmp_size <= max_size){
        size_list.push_back(tmp_size);
        tmp_size *= 8;
    }
    return size_list;
}


Eigen::Tensor<class_xDMRG::Scalar,4> class_xDMRG::find_state_with_greatest_overlap_part_diag2(Eigen::Tensor<Scalar,4> &theta) {
    long shape = theta.size();
    double L    = env_storage->get_length();
    t_ham.tic();
    auto H_local = superblock->get_H_local_rank2();
    auto H_local_sq = superblock->get_H_local_sq_rank2();
    assert(H_local.dimension(0) == shape);
    assert(H_local.dimension(1) == shape);
    t_ham.toc();

    class_eigsolver_arpack<std::complex<double>, Form::GENERAL> arpack_solver;
    Eigen::VectorXcd eigvals;
    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd overlaps;
    std::vector<std::tuple<int,double,double,double,double,double>> result_log;
    Eigen::Map<Textra::VectorType<Scalar>> theta_old (theta.data(),shape);
    Textra::VectorType <Scalar> theta_new;
    Textra::VectorType <Scalar> theta_opt;
    Textra::VectorType <Scalar> theta_eigen;
    Textra::VectorType <Scalar> theta_res;
    double energy_old   = 0;
    double energy_new   = 0;
    Scalar variance_old = 0;
    Scalar variance_new = 0;
    double overlap_new   = 0;

//    int nev = std::min((int)(shape/4), 1);
    int ncv = std::min((int)(shape/2), 4);
    long best_state_idx, worst_state_idx;
    double max_overlap, min_overlap, sq_sum_overlap;
    double subspace_quality;
    double offset;
    double overlap_threshold          = 1 - 1e-4; //1.0/std::sqrt(2); //Slightly less than 1/sqrt(2), in case that the choice is between cat states.
    double subspace_quality_threshold = 1e-2;
    for (auto &nev : generate_size_list(shape)){
        t_eig.tic();
        arpack_solver.eig_shift_invert(H_local.data(), shape, nev,ncv,energy_now*L, Ritz::LM, true, false, theta_res.data());
        t_eig.toc();
        eigvals         = Eigen::Map<const Eigen::VectorXcd> (arpack_solver.get_eigvals().data(),arpack_solver.GetNevFound());
        eigvecs         = Eigen::Map<const Eigen::MatrixXcd> (arpack_solver.get_eigvecs().data(),arpack_solver.Rows(),arpack_solver.Cols());
        overlaps         = (theta_old.adjoint() * eigvecs).cwiseAbs();
        max_overlap      = overlaps.maxCoeff(&best_state_idx);
        min_overlap      = overlaps.minCoeff(&worst_state_idx);
        sq_sum_overlap   = overlaps.cwiseAbs2().sum();
        subspace_quality = 1.0 - sq_sum_overlap;
        offset           = energy_target - eigvals(best_state_idx).real()/L;
        result_log.emplace_back(nev, max_overlap,min_overlap,sq_sum_overlap,std::log10(subspace_quality),t_eig.get_last_time_interval());

        if(max_overlap >= overlap_threshold ){break;}
        if(subspace_quality < subspace_quality_threshold){break;}
        ncv = std::min((int)(shape/2), nev*16);
    }


    if (subspace_quality >= subspace_quality_threshold and max_overlap < overlap_threshold){
//        std::cerr << "          Running full diagonalization. Wall time: " << t_tot.get_age() << std::endl;
        class_eigsolver_armadillo solver;
        Eigen::VectorXd eigvals_real(shape);
        eigvecs.resize(shape,shape);
        t_eig.tic();
        solver.eig_sym(shape, H_local.data(), eigvals_real.data(), eigvecs.data());
        eigvals = eigvals_real;
        t_eig.toc();
        overlaps         = (theta_old.adjoint() * eigvecs).cwiseAbs();
        max_overlap      = overlaps.maxCoeff(&best_state_idx);
        min_overlap      = overlaps.minCoeff(&worst_state_idx);
        sq_sum_overlap   = overlaps.cwiseAbs2().sum();
        subspace_quality = 1.0 - sq_sum_overlap;
        offset           = energy_target - eigvals(best_state_idx).real()/L;
        result_log.emplace_back(shape, max_overlap,min_overlap,sq_sum_overlap,std::log10(subspace_quality),t_eig.get_last_time_interval());
//        std::cerr << "          Full diag complete.           Wall time: " << t_tot.get_age() << std::endl;
    }

    int nev = eigvals.size();


    std::cout << setprecision(16);
    if(nev >= 2){
        double sparcity = Textra::Tensor2_to_Matrix(H_local).real().nonZeros() / (double)H_local.size();
        std::cout << "WARNING: Partial diagonlization -- overlap too small. Starting subspace optimization \n"
                     << std::setprecision(10)
                     << "      max overlap : "    << max_overlap << std::endl
                     << "      position    : "    << env_storage->get_position() << std::endl
                     << "      chi         : "    << chi_temp << std::endl
                     << "      shape       : "    << shape    << " x " << shape << std::endl
                     << "      sparcity    : "    << sparcity << std::endl
                     << std::endl<< std::flush;


        std::cout << std::setprecision(16)
                  << "      " <<setw(12) << left << "n eigvecs"
                  << "      " <<setw(24) << left << "max overlap"
                  << "      " <<setw(24) << left << "min overlap"
                  << "      " <<setw(32) << left << "eps = Σ_i |<state_i|old>|^2"
                  << "      " <<setw(32) << left << "quality = log10(1 - eps)"
                  << "      " <<setw(24) << left << "Wall Time"
                  << std::endl;
        for(auto &log : result_log){
            std::cout << std::setprecision(16)
                      << "      " << setw(12) << left << std::get<0>(log)
                      << "      " << setw(24) << left << std::get<1>(log)
                      << "      " << setw(24) << left << std::get<2>(log)
                      << "      " << setw(32) << left << std::get<3>(log)
                      << "      " << setw(32) << left << std::get<4>(log)
                      << "      " << setw(24) << left << std::get<5>(log)
                      << std::endl;
        }
        std::cout << std::endl << std::flush;

        //Should really use xstart as the projection towards the previous theta, not best overlapping!
        Eigen::VectorXd xstart = (theta_old.adjoint() * eigvecs).cwiseAbs2().normalized();


        class_tic_toc t_lbfgs(true,5,"lbfgs");
        std::vector<std::tuple<std::string,double,Scalar,double,size_t,double>> opt_log;
        {
            t_lbfgs.tic();
            Textra::VectorType<Scalar> theta_0 = eigvecs * xstart;
            size_t iter_0 = 0;
            double energy_0 = xstart.cwiseAbs2().cwiseProduct(eigvals.real()).sum() / L;
            Scalar variance_0 = std::log10(
                    ((theta_0.adjoint() * Textra::Tensor2_to_Matrix(H_local_sq) * theta_0).sum() -
                     energy_0 * energy_0 * L * L) / L);
            double overlap_0 = (theta_old.adjoint() * theta_0).cwiseAbs().sum();
            t_lbfgs.toc();
            opt_log.emplace_back("Start (best overlap)", energy_0, variance_0, overlap_0, iter_0, t_lbfgs.get_last_time_interval());
        }


        {
            using namespace LBFGSpp;
            class_xDMRG_functor functor
                    (shape,nev,
                     H_local.data(),
                     H_local_sq.data(),
                     eigvecs.data(),
                     eigvals.data()
                    );
            double threshold = 1e-3;
            LBFGSpp::LBFGSParam<double> param;
            param.max_iterations = 100;
            param.epsilon        = threshold;
            param.delta          = threshold;
            param.past           = 1;

            // Create solver and function object
            LBFGSpp::LBFGSSolver<double> solver_3(param);
            // x will be overwritten to be the best point found
            double fx;
            t_lbfgs.tic();
            int niter = solver_3.minimize(functor, xstart, fx);
            t_lbfgs.toc();
            xstart.normalize();
            theta_new    = eigvecs * xstart;
            energy_new   = xstart.cwiseAbs2().cwiseProduct(eigvals.real()).sum() / L;
            variance_new = std::log10(
                    ((theta_new.adjoint() * Textra::Tensor2_to_Matrix(H_local_sq) * theta_new).sum() -
                     energy_new * energy_new * L * L) / L);
            overlap_new = (theta_old.adjoint() * theta_new).cwiseAbs().sum();
            opt_log.emplace_back("LBFGS++", energy_new, variance_new, overlap_new, niter, t_lbfgs.get_last_time_interval());
        }






        std::cout << std::setprecision(16)
                  <<"       "<< setw(24) << left << "Algorithm"
                  <<"       "<< setw(24) << left << "energy"
                  <<"       "<< setw(44) << left << "variance"
                  <<"       "<< setw(24) << left << "overlap"
                  <<"       "<< setw(8)  << left << "iter"
                  <<"       "<< setw(20) << left << "time"
                  << std::endl;
        for(auto &log : opt_log){
            std::cout << std::setprecision(16)
                      << "      " <<setw(24) << left << std::get<0>(log)
                      << "      " <<setw(24) << left << std::get<1>(log)
                      << "      " <<setw(44) << left << std::get<2>(log)
                      << "      " <<setw(24) << left << std::get<3>(log)
                      << "      " <<setw(8)  << left << std::get<4>(log)
                      << "      " <<setw(20) << left << std::get<5>(log)
                      << std::endl;
        }
        std::cout << std::endl;


    }else {
        theta_new       = eigvecs.col(best_state_idx);
        energy_new      = eigvals(best_state_idx).real()/L;
    }
    theta_res  = theta_new;
    energy_now = energy_new;

    double energy_ubound = energy_target + 0.1*(energy_max-energy_min);
    double energy_lbound = energy_target - 0.1*(energy_max-energy_min);
    if(energy_now < energy_lbound or energy_now > energy_ubound){
        std::cout << "WARNING: Partial diagonlization -- Energy far from mid-spectrum: " << "Energy =  " << energy_now << " | target = " << energy_target << std::endl;
    }
    return Textra::Matrix_to_Tensor(theta_res, theta.dimensions());
}





void class_xDMRG::initialize_chain() {
    while(true){
        env_storage_insert();
        if (superblock->environment_size + 2ul < (unsigned long) max_length) {
            enlarge_environment();
            swap();
        } else {
            break;
        }
    }
}

void class_xDMRG::reset_chain_mps_to_random_product_state(std::string parity) {
    std::cout << "Resetting to random product state" << std::endl;
    assert(env_storage->get_length() == max_length);

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
    assert(env_storage->get_length() == max_length);
    for (auto &mpo : env_storage->ref_MPO_L()){
        mpo->randomize_hamiltonian();
    }
    for (auto &mpo : env_storage->ref_MPO_R()){
        mpo->randomize_hamiltonian();
    }

    superblock->HA = env_storage->get_MPO_L().back()->clone();
    superblock->HB = env_storage->get_MPO_R().front()->clone();
//    std::cout << "HA initial = \n"
//              << superblock->HA->get_parameter_values()[0] << ' '
//              << superblock->HA->get_parameter_values()[1] << ' '
//              << superblock->HA->get_parameter_values()[2] << ' '
//            << "\n"
//            << env_storage->get_MPO_L().back()->get_parameter_values()[0]<< ' '
//            << env_storage->get_MPO_L().back()->get_parameter_values()[1]<< ' '
//            << env_storage->get_MPO_L().back()->get_parameter_values()[2]<< ' '
//            << std::endl;
//
//
//    std::cout << "HB initial = \n"
//              << superblock->HB->get_parameter_values()[0] << ' '
//              << superblock->HB->get_parameter_values()[1] << ' '
//              << superblock->HB->get_parameter_values()[2] << ' '
//            << "\n"
//            << env_storage->get_MPO_R().front()->get_parameter_values()[0] << ' '
//            << env_storage->get_MPO_R().front()->get_parameter_values()[1] << ' '
//            << env_storage->get_MPO_R().front()->get_parameter_values()[2] << ' '
//            << std::endl;

    iteration = env_storage->reset_sweeps();

//
//
//    while(true) {
//        superblock->HA->randomize_hamiltonian();
//        superblock->HB->randomize_hamiltonian();
//        env_storage_overwrite_local_MPO();
//        // It's important not to perform the last step.
//        if(iteration > 1) {break;}
//        enlarge_environment(env_storage->get_direction());
//        env_storage_move();
//        iteration = env_storage->get_sweeps();
//    }
//    iteration = env_storage->reset_sweeps();
}

void class_xDMRG::find_energy_range() {
    std::cout << "Finding energy range" << std::endl;
    assert(env_storage->get_length() == max_length);
    int max_sweeps_during_f_range = 5;
    int chi_during_f_range   = 4;
    iteration = env_storage->reset_sweeps();


    // Find energy minimum


    while(true) {
        single_DMRG_step(chi_during_f_range, Ritz::SR);
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
        single_DMRG_step(chi_during_f_range, Ritz::LR);
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

void class_xDMRG::store_table_entry_to_file(){
    if (Math::mod(iteration, store_freq) != 0) {return;}
    if (not env_storage->position_is_the_middle()) {return;}
    if (store_freq == 0){return;}
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

void class_xDMRG::store_chain_entry_to_file(){
    if (Math::mod(iteration, store_freq) != 0) {return;}
    if (store_freq == 0){return;}
    if (not (env_storage->get_direction() == 1 or env_storage->position_is_the_right_edge())){return;}
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

void class_xDMRG::check_convergence_all(){
    if(not env_storage->position_is_the_middle()){return;}
//    if(iteration < 5){return;}
    t_con.tic();
    check_convergence_variance_mpo();
//    check_convergence_entanglement();
    update_bond_dimension();
    if(variance_mpo_has_converged)
    {
        simulation_has_converged = true;
    }

    if (variance_mpo_has_saturated and bond_dimension_has_reached_max){
        simulation_has_to_stop = true;
    }


    t_con.toc();
}

void class_xDMRG::initialize_constants(){
    using namespace settings;
    max_length   = xdmrg::max_length;
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