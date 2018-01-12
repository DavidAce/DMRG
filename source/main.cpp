/*! \file */
//#define EIGEN_USE_MKL_ALL
#include <sim_parameters/n_sim_settings.h>
#include <sim_parameters/n_model.h>
#include <mps_routines/class_algorithms.h>
#include <gitversion.h>
#include <IO/class_file_reader.h>


using namespace std;
using namespace Eigen;
using namespace Textra;



/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main() {

    // Print current Git status
    cout << "Branch: "          + GIT::BRANCH +
            " | Commit hash: "  + GIT::COMMIT_HASH +
            " | Revision: "     + GIT::REVISION << endl;

    //Overwrite default settings with the ones in the input file or your own.
    class_file_reader input("input/input.cfg");
    settings::fes_idmrg::chi_min = input.find_parameter<int>("fes_idmrg::chi_min");
    settings::fes_idmrg::chi_max = input.find_parameter<int>("fes_idmrg::chi_max");
    settings::hdf5::save_to_file = input.find_parameter<bool>("hdf5::save_to_file");
    settings::console::verbosity = input.find_parameter<int>("console::verbosity");
    settings::console::graphics  = false;
    settings::profiling::on      = true;
    //Change some model parameters if you don't like the default ones

    //Initialize the algorithm class
    class_algorithms algorithms;

    //Run the algorithms
//    algorithms.iDMRG();
//    algorithms.fDMRG();
//    algorithms.iTEBD();
//    algorithms.FES_iTEBD();
    algorithms.FES_iDMRG();
    return 0;
}



