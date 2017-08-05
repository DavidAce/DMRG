/*! \file */
#include <class_hdf5.h>
#include <class_tic_toc.h>
#include <n_tensor_extra.h>
#include <class_superblock.h>
#include <class_storage.h>
#include <class_algorithms.h>
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
    //Change some parameters if you don't like the default values
    algorithms.params.max_idmrg_length  = 100;
    algorithms.params.max_fdmrg_sweeps  = 4;
    algorithms.params.max_itebd_steps   = 100000;
    algorithms.params.delta_t           = 0.0001;
    algorithms.params.SVDThreshold      = 1e-10;
    algorithms.params.chi               = 10;
    algorithms.params.chi_max           = 30;
    algorithms.params.increasing_chi    = true;
    algorithms.params.verbosity         = 0;
    algorithms.params.graphics          = true;
    algorithms.params.profiling         = true;

    //Run the algorithms
    algorithms.iDMRG(superblock, storage);
    algorithms.fDMRG(superblock, storage);
    algorithms.iTEBD(superblock);

    hdf5.write_to_hdf5(superblock.H.asMatrix,        "Hamiltonian");
    hdf5.write_to_hdf5(superblock.H.asTensor4,       "Hamiltonian MPO");
    hdf5.write_to_hdf5(superblock.H.asTimeEvolution, "Hamiltonian time evolution");
    return 0;
}



