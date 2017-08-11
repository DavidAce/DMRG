/*! \file */
#include <class_superblock.h>
#include <class_storage.h>
#include <class_algorithms.h>
#include <n_settings.h>
using namespace std;
using namespace Eigen;
using namespace Textra;

/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/
int main() {

    std::cout << "Current git revision: " << GIT_REVISION << '\n';

    //Change some settings if you don't like the default values
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
    return 0;
}



