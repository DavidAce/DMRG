/*! \file */

#include <sim_parameters/nmspc_sim_settings.h>
#include <sim_parameters/nmspc_model.h>
#include <algorithms/class_algorithm_launcher.h>
#include <gitversion.h>
#include <IO/class_file_reader.h>


/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main(int argc, char* argv[]) {

    // Print current Git status
    std::cout << "Git Branch: " + GIT::BRANCH +
            " | Commit hash: "  + GIT::COMMIT_HASH +
            " | Revision: "     + GIT::REVISION << std::endl << std::flush;

    //If input file is given as command line argument, then use that.
    std::string inputfile = argc > 1 ? std::string(argv[0]) : std::string("input.cfg") ;
    class_file_reader indata(inputfile);
    settings::load_from_file(indata);

    //Initialize the algorithm class
    //This class stores simulationdata automatically to a file specified in the input file
    class_algorithm_launcher launcher;

    //Run the algorithms
    launcher.run_algorithms();


    return 0;
}



