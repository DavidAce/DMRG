/*! \file */
#include <class_hdf5.h>
#include <class_tic_toc.h>
#include <n_tensor_extra.h>
#include <class_superblock.h>
#include <class_storage.h>
#include <class_algorithms.h>
#include <n_settings.h>
using namespace std;
using namespace Eigen;
using namespace Textra;


int main() {
    /*!
        \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
        \return an integer 0 upon exit success
    */


    //Initialize the algorithm collection, the superblock, storage and hdf5 disk storage classes
    class_algorithms algorithms;
    class_superblock superblock;
    class_storage storage;
    class_hdf5 hdf5(string("myfile.h5"), string("../output"), true);


    //Change some settings if you don't like the default values
    settings::max_idmrg_length  = 100;
    settings::max_fdmrg_sweeps  = 4;
    settings::max_itebd_steps   = 1000;
    settings::max_fes_steps     = 1000;
    settings::SVDThreshold      = 1e-8;
    settings::chi               = 10;
    settings::max_fes_chi       = 15;
    settings::verbosity         = 3;
    settings::graphics          = true;
    settings::profiling         = true;

    //Run the algorithms
    algorithms.iDMRG(superblock, storage, hdf5);
    algorithms.fDMRG(superblock, storage, hdf5);
    algorithms.iTEBD(superblock, hdf5);
    algorithms.FES(superblock,   hdf5);


    /*! \todo  Measure finite-entanglement scaling:
     *  - Grow the system once with iDMRG
     *  - Do fDMRG for some sweeps
     *  - Store mps and all observables to hdf5
     *  - repeat
     *
     *
     *  .... Or maybe use itebd as explained in the paper Jens H sent.
     */

    //Write results to file
    hdf5.write_to_hdf5(superblock.H.asMatrix,        "/Hamiltonian/H");
    hdf5.write_to_hdf5(superblock.H.asTensor4,       "/Hamiltonian/MPO");
    hdf5.write_to_hdf5(superblock.H.asTimeEvolution, "/Hamiltonian/time evolution");
    return 0;
}



