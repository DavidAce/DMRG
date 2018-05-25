//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <IO/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_finite_chain_sweeper.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <general/nmspc_tensor_extra.h>
#include "class_xDMRG.h"

using namespace std;
using namespace Textra;

class_xDMRG::class_xDMRG(std::shared_ptr<class_hdf5_file> hdf5_)
        : class_base_algorithm(std::move(hdf5_), "xDMRG",SimulationType::xDMRG) {
    initialize_constants();
    table_xdmrg = std::make_unique<class_hdf5_table<class_table_dmrg>>(hdf5, sim_name,sim_name);
    env_storage = std::make_shared<class_finite_chain_sweeper>(max_length, superblock, hdf5 );
    measurement  = std::make_shared<class_measurement>(superblock, env_storage, sim_type);
    superblock->HA->set_random_field(r_strength);
    superblock->HB->set_random_field(r_strength);
    initialize_state(settings::model::initial_state);

}



void class_xDMRG::run() {
    if (!settings::xdmrg::on) { return; }
    ccout(0) << "\nStarting " << sim_name << " simulation" << std::endl;
    t_tot.tic();
//    chi_temp = chi_max;
    initialize_random_chain();
    while(true) {
        single_xDMRG_step(chi_temp);
        env_storage_overwrite_MPS();         //Needs to occurr after update_MPS...
        store_table_entry_to_file();
        store_profiling_to_file();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if(iteration >= max_sweeps) {break;}
        enlarge_environment(env_storage->get_direction());
        env_storage_move();
        update_chi();
        iteration = env_storage->get_sweeps();
    }
    t_tot.toc();
    print_status_full();
    print_profiling();
    env_storage->print_storage();
    env_storage->print_hamiltonians();
    measurement->compute_finite_norm();
    measurement->compute_finite_energy();
    measurement->compute_energy_variance_finite_chain();

}

auto class_xDMRG::find_greatest_overlap(Eigen::Tensor<Scalar,4> &theta){
//    Eigen::Tensor<Scalar,1> theta = superblock->MPS->get_theta().reshape(superblock->shape1);
    Eigen::Tensor<Scalar,2> H_local =
            superblock->Lblock->block
            .contract(superblock->HA->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
            .contract(superblock->HB->MPO      , Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
            .contract(superblock->Rblock->block, Textra::idx({4},{2}))
            .shuffle(Textra::array8{2,0,4,6,3,1,5,7})
            .reshape(Textra::array2{superblock->shape1[0], superblock->shape1[0]});
    Eigen::SparseMatrix<Scalar> H_sparse = Eigen::Map<Textra::MatrixType<Scalar>>(H_local.data(),H_local.dimension(0),H_local.dimension(1)).sparseView();
    Eigen::Map<Textra::VectorType<Scalar>> theta_vector (theta.data(),theta.size());
    //Sparcity seems to be about 5-15%.
    Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<Scalar>> es(H_sparse);
    Textra::VectorType<Scalar> overlaps = theta_vector.transpose() * es.eigenvectors();
    int best_state;
    overlaps.cwiseAbs().maxCoeff(&best_state);
    Textra::MatrixType<Scalar> state = es.eigenvectors().col(best_state);

    return Textra::Matrix_to_Tensor(state, superblock->shape4);
}


void class_xDMRG::single_xDMRG_step(long chi_max) {
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();
    superblock->set_current_dimensions();
    t_opt.tic();
    superblock->MPS->theta = superblock->MPS->get_theta();
    superblock->MPS->theta = find_greatest_overlap(superblock->MPS->theta);
    t_opt.toc();
    t_svd.tic();
    superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, chi_max, settings::precision::SVDThreshold);
    t_svd.toc();
    measurement->is_measured = false;
    t_sim.toc();
}



void class_xDMRG::initialize_random_chain() {
    rn::seed((unsigned long)seed);
    while (true) {
        long d    = superblock->MPS->GA.dimension(0);
        long chiB = superblock->MPS->GA.dimension(1);
        long chiA = superblock->MPS->GA.dimension(2);
        superblock->MPS->theta = Textra::Matrix_to_Tensor(Eigen::MatrixXd::Random(d*chiB,d*chiA).cast<Scalar>(),d,chiB,d,chiA);
        superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, chiA, settings::precision::SVDThreshold);

        superblock->HA->set_random_field(r_strength);
        superblock->HB->set_random_field(r_strength);
        print_status_update();
        env_storage_insert();
        if (superblock->environment_size + 2ul < (unsigned long) max_length) {
            enlarge_environment();
            swap();
        } else {
            break;
        }
    }
}



void class_xDMRG::store_table_entry_to_file(){
    if (Math::mod(iteration, store_freq) != 0) {return;}
    compute_observables();
    t_sto.tic();
    table_xdmrg->append_record(
            iteration,
            measurement->get_chain_length(),
            env_storage->get_position(),
            measurement->get_chi(),
            chi_max,
            measurement->get_energy1(),
            measurement->get_energy2(),
            measurement->get_energy3(),
            measurement->get_energy4(),
            measurement->get_energy5(),
            measurement->get_energy6(),
            measurement->get_variance1(),
            measurement->get_variance2(),
            measurement->get_variance3(),
            measurement->get_variance4(),
            measurement->get_variance5(),
            measurement->get_variance6(),
            measurement->get_entanglement_entropy(),
            measurement->get_truncation_error(),
            t_tot.get_age());

    t_sto.toc();
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

void class_xDMRG::update_chi(){
    t_chi.tic();
    if(entropy_has_converged()) {
        simulation_has_converged = chi_temp == chi_max;
        chi_temp = chi_grow ? std::min(chi_max, chi_temp + 4) : chi_max;
    }
    if(not chi_grow){
        chi_temp = chi_max;
    }
    t_chi.toc();
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
        t_chi.print_time_w_percent(t_parent);
    }
}