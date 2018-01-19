//
// Created by david on 2018-01-18.
//

#include <sim_parameters/n_sim_settings.h>
#include <IO/class_multidata_buffer.h>
#include "class_infinite_DMRG.h"
using namespace std;
using namespace Textra;

void class_infinite_DMRG::run() {
/*!
 * \fn void iDMRG(class_superblock &superblock, class_storage &S, int max_length)
 * \brief Infinite DMRG, grows the chain from 2 up to `max_idmrg::length` particles.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param storage A class that stores current MPS and environments at each iteration.
 * \param max_length Maximum chain length after which the algorithm stops.
 */

        if(!settings::idmrg::on){return;}
        ccout(0) << "Starting normal simulation using infinite-DMRG " << std::endl;
        using namespace settings::profiling;
        t_tot.set_properties(on, precision, "iDMRG Total Time           ");
        t_eig.set_properties(on, precision, "iDMRG Eigenvalue solver    ");
        t_svd.set_properties(on, precision, "iDMRG SVD Truncation       ");
        t_env.set_properties(on, precision, "iDMRG Enlarge environment  ");
        t_sto.set_properties(on, precision, "iDMRG Store MPS            ");
        t_mps.set_properties(on, precision, "iDMRG Update MPS           ");

        t_tot.tic();
        class_superblock superblock;
        class_observables observables (superblock, SimulationType::iDMRG);

        while(superblock.chain_length < max_length){
            class_multidata_buffer container(hdf5, "iDMRG/L", superblock.chain_length);
            single_DMRG_step(superblock, chi_max);
            container.push_back(observables);
            t_env.tic();    superblock.enlarge_environment();                       t_env.toc();
            superblock.swap_AB();

        }
        t_tot.toc();
        observables.print_status_full();
        t_eig.print_time_w_percent();
        t_svd.print_time_w_percent();
        t_env.print_time_w_percent();
        t_sto.print_time_w_percent();
        t_mps.print_time_w_percent();
        t_tot.print_time();
        cout << endl;
    }



void class_infinite_DMRG::run(class_superblock &superblock) {
/*!
 * \fn void iDMRG(class_superblock &superblock, class_storage &S, int max_length)
 * \brief Infinite DMRG, grows the chain from 2 up to `max_idmrg::length` particles.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 * \param storage A class that stores current MPS and environments at each iteration.
 * \param max_length Maximum chain length after which the algorithm stops.
 */

    if(!settings::idmrg::on){return;}
    ccout(0) << "Starting normal simulation using infinite-DMRG " << std::endl;
    using namespace settings::profiling;
    t_tot.set_properties(on, precision, "iDMRG Total Time           ");
    t_eig.set_properties(on, precision, "iDMRG Eigenvalue solver    ");
    t_svd.set_properties(on, precision, "iDMRG SVD Truncation       ");
    t_env.set_properties(on, precision, "iDMRG Enlarge environment  ");
    t_sto.set_properties(on, precision, "iDMRG Store MPS            ");
    t_mps.set_properties(on, precision, "iDMRG Update MPS           ");

    t_tot.tic();
    class_observables observables (superblock, SimulationType::iDMRG);
    class_multidata_buffer container(hdf5, "iDMRG/L", superblock.chain_length);

    while(superblock.chain_length < max_length){
        single_DMRG_step(superblock, chi_max);
        container.push_back(observables);
        t_env.tic();    superblock.enlarge_environment();                       t_env.toc();
        superblock.swap_AB();

    }
    t_tot.toc();
    observables.print_status_full();
    t_eig.print_time_w_percent();
    t_svd.print_time_w_percent();
    t_env.print_time_w_percent();
    t_sto.print_time_w_percent();
    t_mps.print_time_w_percent();
    t_tot.print_time();
    cout << endl;
}


