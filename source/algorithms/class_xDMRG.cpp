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
#include <mps_routines/class_environment_storage.h>
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
        : class_algorithm_base(std::move(hdf5_), "xDMRG", "xDMRG",SimulationType::xDMRG) {
    env_storage  = std::make_shared<class_environment_storage>(max_length, superblock, hdf5);

}



void class_xDMRG::run() {
    if (!settings::xdmrg::on) { return; }
    ccout(0) << "\nStarting " << table_name << " simulation" << std::endl;
    t_tot.tic();
    chi_temp = chi_max;
//    while(superblock->chain_length < max_length and !simulation_has_converged){
//        single_xDMRG_step(chi_temp);
//        print_status_update();
//        store_table_entry();
//
//        position = enlarge_environment();
//        iteration++;
//        update_chi();
//        swap();
//    }
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
    }
    t_tot.toc();
    print_status_full();
    print_profiling();

}

void class_xDMRG::single_xDMRG_step(long chi_max) {
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();
    superblock->update_bond_dimensions();
    t_eig.tic();    find_greatest_overlap();    t_eig.toc();
    t_svd.tic();    superblock->truncate         (chi_max, settings::precision::SVDThreshold);                   t_svd.toc();
    t_mps.tic();    superblock->update_MPS();                                                                    t_mps.toc();
    t_sim.toc();
}


void class_xDMRG::find_greatest_overlap(){
    Eigen::Tensor<Scalar,1> theta = superblock->MPS->get_theta().shuffle(array4{0,2,1,3}).reshape(superblock->shape1);
    //Best yet! The sparcity of the effective hamiltonian (Lblock MM Rblock) is about 58% nonzeros.
    Eigen::Tensor<Scalar,2> H_local =
            superblock->Lblock->block
            .contract(superblock->H->MM ,        Textra::idx<1>({2},{0}))//  idx<3>({1,2,3},{0,4,5}))
            .contract(superblock->Rblock->block, Textra::idx<1>({2},{2}))
            .shuffle(Textra::array8{2,0,3,6,4,1,5,7})
            .reshape(Textra::array2{superblock->shape1[0], superblock->shape1[0]});

    using SparseType = Eigen::SparseMatrix<Scalar>;
    Eigen::Map<Textra::MatrixType<Scalar>> H_matrix     (H_local.data(),H_local.dimension(0),H_local.dimension(1));
    Eigen::Map<Textra::VectorType<Scalar>> theta_matrix (theta.data(),theta.dimension(0));
    SparseType H_sparse = H_matrix.sparseView();

    Eigen::SelfAdjointEigenSolver<SparseType> es(H_sparse);
    Textra::VectorType<double> overlaps = theta_matrix.transpose() * es.eigenvectors();

    int best_state;
    overlaps.cwiseAbs().maxCoeff(&best_state);
    Textra::MatrixType<double> state = es.eigenvectors().col(best_state);

    superblock->ground_state = Textra::Matrix_to_Tensor<Scalar,2>(state, superblock->shape2);
    std::cout << setprecision(16);

//    std::cout << "eigenvalue: "<< es.eigenvalues().coeff(best_state)/superblock->chain_length << std::endl;// <<  es.eigenvalues().transpose()/superblock->chain_length << std::endl;
//    SparseType H_sparse = H_matrix.sparseView();

//    H_sparse = H_sparse.selfadjointView<Eigen::Upper>();
//    std::cout<<setprecision(8) << "H_matrix\n" <<H_matrix - H_matrix.transpose() << std::endl;
//    std::cout << "H_matrix: \n"<< H_matrix << std::endl;
//    std::cout << "theta_matrix: \n"<< theta_matrix.transpose() << std::endl;
//    std::cout << "eigenvectors: \n"<< es.eigenvectors() << std::endl;
//    std::cout << "overlaps:     \n"<< overlaps.transpose() << std::endl;
//    state *= overlaps.array().sign().coeff(best_state);

//    Textra::Tensor<Scalar,0> E =
//            H_local.contract(superblock->ground_state.reshape(superblock->shape1),             Textra::idx<1>({0},{0}))
//                   .contract(superblock->ground_state.reshape(superblock->shape1).conjugate(), Textra::idx<1>({0},{0}) );
////    std::cout << "chosen (" << best_state << "):     \n" << superblock->ground_state << std::endl;
//    exit(0);
//    std::cout << "E check: " << E(0)/superblock->chain_length << std::endl;
}

int class_xDMRG::initialize_random_chain() {
    rn::seed(5);

    while (true) {
        auto r1 = rn::uniform_complex_1();
        auto r2 = rn::uniform_complex_1();
        superblock->MPS->GA(1, 0, 0) = r1.real();
        superblock->MPS->GA(0, 0, 0) = r1.imag();
        superblock->MPS->GB(1, 0, 0) = r2.real();
        superblock->MPS->GB(0, 0, 0) = r2.imag();
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
        t_env.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        t_eig.print_time_w_percent(t_sim);
        t_svd.print_time_w_percent(t_sim);
        t_mps.print_time_w_percent(t_sim);
    }
}
