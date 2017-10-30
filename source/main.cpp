/*! \file */
#include <n_settings.h>
#include <class_algorithms.h>
#include <gitversion.h>


using namespace std;
using namespace Eigen;
using namespace Textra;



/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main() {
    cout << "Branch: "          + GIT::BRANCH +
            " | Commit hash: "  + GIT::COMMIT_HASH +
            " | Revision: "     + GIT::REVISION << endl;

    cout << GIT::BRANCH << endl;
    //Change some settings if you don't like the default values
    settings::hdf5::save_to_file = true;
    settings::console::verbosity = 3;
    settings::console::graphics  = true;
    settings::profiling::on      = true;

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



