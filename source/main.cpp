/*! \file */


#include <sim_parameters/nmspc_sim_settings.h>
#include <sim_parameters/nmspc_model.h>
#include <algorithms/class_algorithm_launcher.h>
#include <gitversion.h>
#include <IO/class_file_reader.h>
//#include <El.hpp>

#ifdef OpenBLAS_AVAILABLE
#include <cblas.h>
#endif


#ifdef OpenMP_AVAILABLE
#include <omp.h>
#endif

#ifdef MKL_AVAILABLE
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl_service.h>
#include <mkl.h>
#endif


/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main(int argc, char* argv[]) {
    int num_threads = 4;
    #ifdef OpenBLAS_AVAILABLE
        openblas_set_num_threads(num_threads);
        std::cout << "Using OpenBLAS with " << openblas_get_num_threads() << " thread(s)" << std::endl;
    #endif

    #ifdef OpenMP_AVAILABLE
        omp_set_num_threads(num_threads);
        Eigen::setNbThreads(num_threads);
        std::cout << "Using OpenMP with " << omp_get_num_threads() << " thread(s)" << std::endl;
    #endif

    #ifdef MKL_AVAILABLE
        mkl_set_num_threads(num_threads);
        std::cout << "Using Intel MKL with " << num_threads << " thread(s)" << std::endl;
    #endif


    // Print current Git status
    std::cout << "Git Branch: " + GIT::BRANCH +
            " | Commit hash: "  + GIT::COMMIT_HASH +
            " | Revision: "     + GIT::REVISION << std::endl << std::flush;

    //Print all given parameters
    //Load input and output files from command line. If none were given use defaults.
    //Normally an output filename is given in the input file. But it can also be given from command line.
    std::string inputfile  = "input.cfg";
    std::string outputfile;
    bool outputfile_given  = false;
    for (int i=0; i < argc; i++){
        std::string arg_string = std::string(argv[i]);
        std::cout <<"Input argument [" << i << "] : " << arg_string << std::endl;
        if (arg_string.find(".cfg") != std::string::npos) {inputfile  = arg_string;}
        if (arg_string.find(".h5")  != std::string::npos) {outputfile = arg_string;outputfile_given=true;}
    }
    class_file_reader indata(inputfile);
    settings::load_from_file(indata);

    //If an output filename was given explicitly, overwrite the default , if it was given explicitly in command line.
    settings::hdf5::output_filename = outputfile_given ? outputfile : settings::hdf5::output_filename;

    //Initialize the algorithm class
    //This class stores simulation data_struct automatically to a file specified in the input file
    class_algorithm_launcher launcher;

    //Run the algorithms
    launcher.run_algorithms();


    return 0;
}



