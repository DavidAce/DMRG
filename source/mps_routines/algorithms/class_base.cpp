//
// Created by david on 2018-01-18.
//

#include "class_base.h"
#include <IO/class_hdf5_file.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_measurement.h>
namespace s = settings;
using namespace std;
using namespace Textra;

class_base::class_base(std::shared_ptr<class_hdf5_file>         hdf5_,
                       std::shared_ptr<class_hdf5_table_buffer> table_buffer_,
                       std::shared_ptr<class_superblock>        superblock_,
                       std::shared_ptr<class_measurement>       observables_)
                :hdf5           (std::move(hdf5_)),
                table_buffer    (std::move(table_buffer_)),
                superblock      (std::move(superblock_)),
                observables     (std::move(observables_))
{
    simtype = observables->sim;
//    set_simtype_labels();
    set_profiling_labels();
};


class_base::class_base(){
//    set_simtype_labels();
    set_profiling_labels();
}


void class_base::single_DMRG_step(long chi_max){
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();
    superblock->update_bond_dimensions();
    t_eig.tic();    superblock->find_ground_state(s::precision::eigSteps, s::precision::eigThreshold);    t_eig.toc();
    t_svd.tic();    superblock->truncate         (chi_max,                s::precision::SVDThreshold);    t_svd.toc();
    t_mps.tic();    superblock->update_MPS();                                                             t_mps.toc();
    t_sim.toc();
}

void class_base::single_TEBD_step(long chi_max){
/*!
 * \fn single_iTEBD_step(class_superblock &superblock)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 */
    t_sim.tic();
    superblock->update_bond_dimensions();
    t_evo.tic();    superblock->time_evolve();                                       t_evo.toc();
    t_svd.tic();    superblock->truncate   (chi_max, s::precision::SVDThreshold);    t_svd.toc();
    t_mps.tic();    superblock->update_MPS();                                        t_mps.toc();
    t_sim.toc();
}

//void class_base::set_simtype_labels() {
//        simType2string = {
//            {SimulationType::iDMRG,     "iDMRG"},
//            {SimulationType::fDMRG,     "fDMRG"},
//            {SimulationType::FES_iDMRG, "FES_iDMRG"},
//            {SimulationType::iTEBD,     "iTEBD"},
//            {SimulationType::FES_iTEBD, "FES_iTEBD"},
//            {SimulationType::NONE,      "NONE"}
//    };
//}

void class_base::set_profiling_labels() {
    using namespace settings::profiling;
    t_tot.set_properties(on, precision,"+Total Time              ");
    t_sto.set_properties(on, precision,"\u21B3 Store Results          ");
    t_udt.set_properties(on, precision,"\u21B3 Update Timestep        ");
    t_env.set_properties(on, precision,"\u21B3 Update Environments    ");
    t_sim.set_properties(on, precision,"\u21B3+Simulation             ");
    t_evo.set_properties(on, precision," \u21B3 Time Evolution        ");
    t_eig.set_properties(on, precision," \u21B3 Eig. decomp.          ");
    t_svd.set_properties(on, precision," \u21B3 SVD Truncation        ");
    t_mps.set_properties(on, precision," \u21B3 Update MPS            ");
}

void class_base::enlarge_environment(){
    t_env.tic();
    superblock->enlarge_environment();
    t_env.toc();
}

void class_base::swap(){
    superblock->swap_AB();
}

