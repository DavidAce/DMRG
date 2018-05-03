//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_finite_chain_storage.h>
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
        : class_base_algorithm(std::move(hdf5_), "xDMRG", "xDMRG",SimulationType::xDMRG) {
    env_storage  = std::make_shared<class_finite_chain_storage>(max_length, superblock, hdf5);

}



void class_xDMRG::run() {
    if (!settings::xdmrg::on) { return; }
    ccout(0) << "\nStarting " << table_name << " simulation" << std::endl;
    t_tot.tic();
    chi_temp = chi_max;
    randomness_strength = 1;
    initialize_random_chain();
    while(sweeps < max_sweeps) {
        single_xDMRG_step(chi_temp);
        env_storage_overwrite_MPS();         //Needs to occurr after update_MPS...
        store_table_entry();
        print_status_update();
        iteration++;
        position = enlarge_environment(direction);
        position = env_storage_move();
        update_chi();
//        if(iteration > 10 ){exit(0);}
    }
    t_tot.toc();
    print_status_full();
    print_profiling();

}


auto class_xDMRG::find_greatest_overlap(){
    Eigen::Tensor<Scalar,1> theta = superblock->MPS->get_theta().reshape(superblock->shape1);
    Eigen::Tensor<Scalar,2> H_local =
            superblock->Lblock->block
            .contract(superblock->HA->MPO_zero_site_energy() ,         Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
            .contract(superblock->HB->MPO_zero_site_energy() ,         Textra::idx({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
            .contract(superblock->Rblock->block, Textra::idx({4},{2}))
            .shuffle(Textra::array8{2,0,4,6,3,1,5,7})
            .reshape(Textra::array2{superblock->shape1[0], superblock->shape1[0]});

    using SparseType = Eigen::SparseMatrix<Scalar>;
    Eigen::Map<Textra::MatrixType<Scalar>> H_matrix     (H_local.data(),H_local.dimension(0),H_local.dimension(1));
    Eigen::Map<Textra::VectorType<Scalar>> theta_matrix (theta.data(),theta.dimension(0));
    SparseType H_sparse = H_matrix.sparseView();

    Eigen::SelfAdjointEigenSolver<SparseType> es(H_sparse);
    Textra::VectorType<Scalar> overlaps = theta_matrix.transpose() * es.eigenvectors();

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
    superblock->MPS->theta = find_greatest_overlap();
    t_opt.toc();
    t_svd.tic();
    superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, chi_max, settings::precision::SVDThreshold);
    t_svd.toc();
    t_sim.toc();
}



int class_xDMRG::initialize_random_chain() {
    rn::seed(7);
    while (true) {
        auto r1 = rn::uniform_complex_1();
        auto r2 = rn::uniform_complex_1();
        superblock->MPS->GA(1, 0, 0) = r1.real();
        superblock->MPS->GA(0, 0, 0) = r1.imag();
        superblock->MPS->GB(1, 0, 0) = r2.real();
        superblock->MPS->GB(0, 0, 0) = r2.imag();
        superblock->H->HA = superblock->H->compute_H_MPO_custom_field(
                rn::uniform_double(-randomness_strength, randomness_strength), 0.0);
        superblock->H->HB = superblock->H->compute_H_MPO_custom_field(
                rn::uniform_double(-randomness_strength, randomness_strength), 0.0);
        position = env_storage_insert();
        if (superblock->chain_length < max_length) {
            enlarge_environment();
            position++;
        } else {
            break;
        }
    }
    return position;
}



int class_xDMRG::env_storage_insert() {
    t_ste.tic();
    int position = env_storage->insert();
    t_ste.toc();
    return position;
}

void class_xDMRG::env_storage_overwrite_MPS(){
    t_ste.tic();
    env_storage->overwrite_MPS();
    t_ste.toc();
}
int class_xDMRG::env_storage_move(){
    t_ste.tic();
    int position = env_storage->move(direction, sweeps);
    t_ste.toc();
    return position;
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