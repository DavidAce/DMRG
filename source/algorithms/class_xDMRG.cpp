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
//
//        if(iteration < full_iters){
//            if(iteration == 0 and env_storage->position_is_the_middle()){std::cout << "--FULL DIAGONALIZATION MODE--" << std::endl;}
//            single_xDMRG_step(xDMRG_Mode::FULL);
//        }else{
//            if(iteration == full_iters and env_storage->position_is_the_middle()){std::cout << "--PARTIAL DIAGONALIZATION MODE--" << std::endl;}
//            single_xDMRG_step(xDMRG_Mode::PARTIAL);
//        }
//        env_storage->print_hamiltonians();
//        exit(1);
        env_storage_overwrite_local_ALL();
        store_table_entry_to_file();
        store_chain_entry_to_file();
        store_profiling_to_file();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if(env_storage->position_is_the_middle()) {
            check_convergence_all();
            if (iteration >= full_iters and (iteration >= max_sweeps or simulation_has_converged)) {
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
    if(measurement->get_chain_length() < 15){
        hdf5->write_dataset(Textra::to_RowMajor(measurement->mps_chain), sim_name + "/chain/wavefunction");
    }

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
    xDMRG_Mode mode = chi_temp <= 12 ? xDMRG_Mode::FULL : xDMRG_Mode::PARTIAL;

    switch (mode){
        case xDMRG_Mode::FULL:
            theta = find_state_with_greatest_overlap_full_diag(theta);
            t_opt.toc();
            break;
        case xDMRG_Mode::PARTIAL:
            theta = find_state_with_greatest_overlap_part_diag(theta);
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
    Eigen::Tensor<Scalar,2> H_local  =
            superblock->Lblock->block
                    .contract(superblock->HA->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
                    .contract(superblock->HB->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
                    .contract(superblock->Rblock->block, Textra::idx({4},{2}))
                    .shuffle(Textra::array8{3,1,5,7,2,0,4,6})
                    .reshape(Textra::array2{shape, shape});
    assert(H_local.dimension(0) == shape);
    assert(H_local.dimension(1) == shape);
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
    std::cerr << "Armadillo Found overlap: " << setprecision(14) << setw(18) << max_overlap << std::endl;

    Textra::MatrixType<Scalar> state = eigvecs.col(best_state);
    energy_now = eigvals(best_state)/L;
    return Textra::Matrix_to_Tensor(state, theta.dimensions());
}

Eigen::Tensor<class_xDMRG::Scalar,4> class_xDMRG::find_state_with_greatest_overlap_part_diag(Eigen::Tensor<Scalar,4> &theta) {
    long shape = theta.size();
    double L    = env_storage->get_length();
    t_ham.tic();
    Eigen::Tensor<Scalar,2> H_local  =
            superblock->Lblock->block
                    .contract(superblock->HA->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
                    .contract(superblock->HB->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
                    .contract(superblock->Rblock->block, Textra::idx({4},{2}))
                    .shuffle(Textra::array8{3,1,5,7,2,0,4,6})
                    .reshape(Textra::array2{shape, shape});
    assert(H_local.dimension(0) == shape);
    assert(H_local.dimension(1) == shape);
    t_ham.toc();

    long chiA = superblock->MPS->chiA();
    long chiB = superblock->MPS->chiB();
    long d    = superblock->HA->get_spin_dimension();
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
    for (int iter = 0; iter  <  5; iter++){
        t_eig.tic();
//        arpack_solver.eig_shift_invert2(superblock->Lblock->block.data(),
//                                        superblock->Rblock->block.data(),
//                                        superblock->HA->MPO.data(),
//                                        superblock->HB->MPO.data(),
//                                        shape4_theta,
//                                        shape4_mpo,
//                                        nev, ncv,
//                                        energy_now*L,
//                                        Ritz::LM,
//                                        true, true,
//                                        theta.data());
        arpack_solver.eig_shift_invert(H_local.data(), shape, nev,ncv,energy_now*L, Ritz::LM, true,true, theta_res.data());
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
//    if( overlaps.maxCoeff() < 1e-2){
//        std::cerr << "Failed to find a great eigenstate! Choosing state closest in energy instead." << std::endl;
//        (eigvals.array() - energy_target*L).cwiseAbs().minCoeff(&best_state_idx);
////        theta_res       = eigvecs.col(best_state_idx);
//
//    }
    if(max_overlap <  overlap_threshold){
        std::cerr << "WARNING:  Partial diagonlization -- Overlap very small: " << std::setprecision(10) << max_overlap << "      position: " << env_storage->get_position() << std::endl;
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
    std::cerr << "Arpack Found overlap: " << setprecision(14) << setw(18) << max_overlap << std::endl;

    energy_now      = eigvals(best_state_idx).real()/L;
    theta_res       = eigvecs.col(best_state_idx);

    double energy_ubound = energy_target + 0.05*(energy_max-energy_min);
    double energy_lbound = energy_target - 0.05*(energy_max-energy_min);
    if(energy_now < energy_lbound or energy_now > energy_ubound){
        std::cerr << "WARNING: Partial diagonlization -- Energy far from mid-spectrum: " << "Energy =  " << energy_now << " | target = " << energy_target << std::endl;
    }


//    energy_target = energy_now;
//    if (env_storage->position_is_the_middle()){energy_target = energy_now;}
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
    double energy_ubound = energy_target + 0.01*(energy_max-energy_min);
    double energy_lbound = energy_target - 0.01*(energy_max-energy_min);
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
    if(variance_mpo_has_converged )
    {
        simulation_has_converged = true;
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