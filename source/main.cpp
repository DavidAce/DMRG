/*! \file */
#include <sim_parameters/nmspc_sim_settings.h>
#include <sim_parameters/nmspc_model.h>
#include <algorithms/class_algorithm_launcher.h>
#include <gitversion.h>
#include <IO/class_file_reader.h>
//#include <mkl.h>
#include <omp.h>


/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main(int argc, char* argv[]) {

    // Print current Git status
    std::cout << "Git Branch: " + GIT::BRANCH +
            " | Commit hash: "  + GIT::COMMIT_HASH +
            " | Revision: "     + GIT::REVISION << std::endl << std::flush;
//    mkl_set_num_threads(8);
//    omp_set_num_threads(8);
//    Eigen::setNbThreads(8);
//      openblas_set_num_threads(8);

    std::cout << "Eigen Num Threads: " << Eigen::nbThreads( ) << std::endl;
    std::cout << "OMP   Max Threads: " << omp_get_max_threads() << std::endl;
//    std::cout << "MKL   Num Threads: " << mkl_get_num_threads() << std::endl;
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



