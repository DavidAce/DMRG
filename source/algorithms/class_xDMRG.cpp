//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <IO/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_finite_chain_sweeper.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/nmspc_eigsolver_props.h>
#include <general/class_elemental_eigsolver.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <general/nmspc_tensor_extra.h>
#include <general/class_arpack_eigsolver.h>
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
    find_energy_range();
    while(true) {
        single_xDMRG_step(chi_max_temp,energy_mid);
        env_storage_overwrite_local_ALL();
        store_table_entry_to_file();
        store_chain_entry_to_file();
        store_profiling_to_file();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if(iteration >= max_sweeps or simulation_has_converged) {break;}
        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        check_convergence_overall();
        iteration = env_storage->get_sweeps();
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
    measurement->compute_all_observables_from_finite_chain();
    env_storage->write_chain_to_file();
}

void class_xDMRG::single_xDMRG_step(long chi_max, double energy_target) {
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();
    t_opt.tic();
    superblock->MPS->theta = superblock->MPS->get_theta();
    if (iteration <= 5){
        chi_max = std::min(chi_max,12l);
        superblock->MPS->theta = find_state_with_greatest_overlap_full_diag(superblock->MPS->theta, energy_target);
    }else{
        superblock->MPS->theta = find_state_with_greatest_overlap_part_diag(superblock->MPS->theta, energy_target);
    }

    t_opt.toc();
    t_svd.tic();
    superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, chi_max, settings::precision::SVDThreshold);
    t_svd.toc();
    measurement->set_not_measured();
    t_sim.toc();
}


Eigen::Tensor<class_xDMRG::Scalar,4> class_xDMRG::find_state_with_greatest_overlap_full_diag(Eigen::Tensor<Scalar, 4> &theta, double energy_target){
//    Eigen::Tensor<Scalar,1> theta = superblock->MPS->get_theta().reshape(superblock->shape1);
    long shape = theta.size();
    double L   = env_storage->get_length();
    Eigen::Tensor<Scalar,2> H_local  =
            superblock->Lblock->block
            .contract(superblock->HA->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
            .contract(superblock->HB->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
            .contract(superblock->Rblock->block, Textra::idx({4},{2}))
            .shuffle(Textra::array8{3,1,5,7,2,0,4,6})
            .reshape(Textra::array2{shape, shape});
    assert(H_local.dimension(0) == shape);
    assert(H_local.dimension(1) == shape);

    //Sparcity seems to be about 5-15%. Better to do dense.

    class_elemental_eigsolver solver;
    Textra::VectorType<double> eigvals(shape);
    Textra::MatrixType<Scalar> eigvecs(shape,shape);
    solver.eig(shape, H_local.data(),eigvals.data(),eigvecs.data());
    Eigen::Map<Textra::VectorType<Scalar>> theta_vector (theta.data(),shape);
    Textra::VectorType<Scalar> overlaps = theta_vector.conjugate().transpose() * eigvecs;
    long best_state;
    if (overlaps.isApproxToConstant(0.0, 1e-3) or std::abs(energy_at_site - energy_target) > 1e-1){
        best_state = std::lower_bound(eigvals.data(), eigvals.data() + eigvals.size(), energy_target*L) - eigvals.data();
        best_state = std::min(best_state, overlaps.size()-1);
        best_state = std::max(best_state,0l);
//        std::cout << setprecision(10) << "diff: " << std::abs(energy_at_site - energy_target);
//        std::cout << "  best state: " << best_state << "   energy_target = " << energy_target*L << "energy chosen = " << eigvals(best_state)<<  std::endl;
//        std::cout << " eigvals \n " << eigvals << std::endl;

    }else{
        overlaps.cwiseAbs().maxCoeff(&best_state);
    }
    Textra::MatrixType<Scalar> state = eigvecs.col(best_state);
    energy_at_site = eigvals(best_state)/L;
    return Textra::Matrix_to_Tensor(state, theta.dimensions());
}



bool complex_less(const std::complex<double>& lhs, const std::complex<double> &rhs){
    return std::real(lhs) < std::real(rhs);
}



Eigen::Tensor<class_xDMRG::Scalar,4> class_xDMRG::find_state_with_greatest_overlap_part_diag(Eigen::Tensor<Scalar,4> &theta, double energy_target) {
    long shape = theta.size();
    double L    = env_storage->get_length();
    Eigen::Tensor<Scalar,2> H_local  =
            superblock->Lblock->block
                    .contract(superblock->HA->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
                    .contract(superblock->HB->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
                    .contract(superblock->Rblock->block, Textra::idx({4},{2}))
                    .shuffle(Textra::array8{3,1,5,7,2,0,4,6})
                    .reshape(Textra::array2{shape, shape});
    assert(H_local.dimension(0) == shape);
    assert(H_local.dimension(1) == shape);


    int nev = std::min((int)(shape/4), 20);
    int ncv = std::min((int)(shape/2), nev*4);
    class_arpack_eigsolver<std::complex<double>> arpack_solver;
    arpack_solver.eig_shift_invert(H_local.data(), shape, nev,ncv,energy_target, Ritz::LM, true,true);
    Eigen::Map<const Eigen::VectorXcd> eigvals(arpack_solver.ref_eigvals().data(),arpack_solver.GetNevFound());
    Eigen::Map<const Eigen::MatrixXcd> eigvecs(arpack_solver.ref_eigvecs().data(),arpack_solver.Rows(),arpack_solver.Cols());
    Eigen::Map<Textra::VectorType<Scalar>> theta_vector (theta.data(),shape);
    Textra::VectorType<Scalar> overlaps = theta_vector.conjugate().transpose() * eigvecs;
    long best_state;
    if (overlaps.isApproxToConstant(0.0, 1e-3)){
        best_state = std::lower_bound(eigvals.data(), eigvals.data() + eigvals.size(), energy_target*L, complex_less) - eigvals.data();
        best_state = std::min(best_state, overlaps.size()-1);
        best_state = std::max(best_state,0l);
    }else{
        overlaps.cwiseAbs().maxCoeff(&best_state);
    }
    Textra::MatrixType<Scalar> state = eigvecs.col(best_state);
    energy_at_site = eigvals(best_state).real()/L;
    return Textra::Matrix_to_Tensor(state, theta.dimensions());
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


void class_xDMRG::reset_chain_mps_to_random_product_state() {
    std::cout << "Resetting to random product state" << std::endl;
    assert(env_storage->get_length() == max_length);

    iteration = env_storage->reset_sweeps();
    while(true) {
        // Random product state
        long chiA = superblock->MPS->chiA();
        long chiB = superblock->MPS->chiB();
        long d    = superblock->HA->get_spin_dimension();
        superblock->MPS->theta = Textra::Matrix_to_Tensor(Eigen::MatrixXcd::Random(d*chiA,d*chiB),d,chiA,d,chiB);
        //Get a properly normalized initial state.
        superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, 1, settings::precision::SVDThreshold);
        env_storage_overwrite_local_ALL();
//        env_storage->print_storage();
        // It's important not to perform the last step.
        if(iteration > 1) {break;}
        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        iteration = env_storage->get_sweeps();

    }
    iteration = env_storage->reset_sweeps();
}

void class_xDMRG::set_random_fields_in_chain_mpo() {
    std::cout << "Setting random fields in chain" << std::endl;
    assert(env_storage->get_length() == max_length);
    iteration = env_storage->reset_sweeps();
    while(true) {
        superblock->HA->randomize_hamiltonian();
        superblock->HB->randomize_hamiltonian();
        env_storage_overwrite_local_MPO();

        // It's important not to perform the last step.
        if(iteration > 1) {break;}
        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        iteration = env_storage->get_sweeps();
    }
    iteration = env_storage->reset_sweeps();
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
    compute_observables();
    energy_min = measurement->get_energy_mpo();
    std::cout << "Energy minimum = " << energy_min << std::endl;
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
    std::cout << "Energy maximum = " << energy_max << std::endl;

    iteration = env_storage->reset_sweeps();
    energy_mid = 0.5*(energy_max+energy_min);
    energy_at_site = superblock->E_optimal / env_storage->get_length();
    std::cout << "Energy_at_site: " << energy_at_site << std::endl;
    std::cout << "Energy medium  = " << energy_mid << std::endl;
    reset_chain_mps_to_random_product_state();
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
            energy_at_site / measurement->get_chain_length(),
            measurement->get_entanglement_entropy(),
            measurement->get_truncation_error()
    );
    t_sto.toc();
}

void class_xDMRG::check_convergence_overall(){
    if(not env_storage->position_is_the_middle()){return;}
    t_con.tic();
    check_convergence_entanglement(1e-6);
    check_convergence_variance_mpo();
    check_convergence_bond_dimension();
    if(entanglement_has_converged and
       variance_mpo_has_converged and
       bond_dimension_has_converged)
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
    r_strength   = xdmrg::r_strength;
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
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_con.print_time_w_percent(t_parent);
    }
}