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

    //Change some settings if you don't like the default values
    settings::idmrg_max_length  = 100;
    settings::fdmrg_max_sweeps  = 4;
    settings::itebd_max_steps   = 1000;
    settings::fes_max_steps     = 1000;
    settings::SVDThreshold      = 1e-8;
    settings::chi               = 10;
    settings::fes_max_chi       = 15;
    settings::hdf5_save_to_file = true;
    settings::console_verbosity = 3;
    settings::console_graphics  = true;
    settings::profiling_on      = true;

    //Initialize the algorithm collection, the superblock, storage and hdf5 disk storage classes
    class_algorithms algorithms;
    class_superblock superblock;
    class_storage storage;
//    class_hdf5 hdf5(string("myfile.h5"), string("../output"), true);



    //Run the algorithms
    algorithms.iDMRG(superblock, storage);
    algorithms.fDMRG(superblock, storage);
    algorithms.iTEBD(superblock);
    algorithms.FES(superblock);


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
//    hdf5.write_to_file(superblock.H.asMatrix,        "/Hamiltonian/H");
//    hdf5.write_to_file(superblock.H.asTensor4,       "/Hamiltonian/MPO");
//    hdf5.write_to_file(superblock.H.asTimeEvolution, "/Hamiltonian/time evolution");
    return 0;
}



