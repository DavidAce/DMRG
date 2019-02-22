/*! \file */


#include <sim_parameters/nmspc_sim_settings.h>
#include <sim_parameters/nmspc_model.h>
#include <algorithms/class_algorithm_launcher.h>
#include <gitversion.h>
#include <io/class_settings_reader.h>
#include <io/class_hdf5_file.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
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

    spdlog::set_default_logger(spdlog::stdout_color_mt("DMRG"));
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v");
    spdlog::set_level(spdlog::level::trace);

    int openblas_num_threads = 1;
    #ifdef OpenBLAS_AVAILABLE
        openblas_set_num_threads(openblas_num_threads);
    std::cout << "OpenBLAS compiled with mode " << openblas_get_parallel()
              << " for target " << openblas_get_corename()
              << " with config " << openblas_get_config()
              << ". Running with " << openblas_get_num_threads() << " thread(s)" << std::endl;
    #endif

    #ifdef OpenMP_AVAILABLE
        Eigen::initParallel();
//        omp_set_num_threads(OpenMP_NUM_THREADS);
//        omp_set_dynamic(0);
//        Eigen::setNbThreads(OpenMP_NUM_THREADS);
//        Eigen::setNbThreads(0);
        spdlog::info("Using Eigen  with {} threads",Eigen::nbThreads( )  );
        spdlog::info("Using OpenMP with {} threads",omp_get_max_threads());
        #ifdef MKL_AVAILABLE
//            mkl_set_num_threads(OpenMP_NUM_THREADS);
            spdlog::info("Using Intel MKL with {} threads", mkl_get_max_threads());
        #endif
    #endif




    // Print current Git status
    spdlog::info("Git branch      : {}",GIT::BRANCH);
    spdlog::info("    commit hash : {}",GIT::COMMIT_HASH);
    spdlog::info("    revision    : {}",GIT::REVISION);

    //Print all given parameters
    //Load input and output files from command line. If none were given use defaults.
    //Normally an output filename is given in the input file. But it can also be given from command line.
    std::string inputfile  = "input.cfg";
    std::string outputfile = "output.h5";
    for (int i=0; i < argc; i++){
        std::string arg_string = std::string(argv[i]);
        spdlog::info("Input argument {} : {}",i,arg_string);
        if (arg_string.find(".cfg") != std::string::npos) {inputfile  = arg_string;}
        if (arg_string.find(".h5")  != std::string::npos) {outputfile = arg_string;}
    }
    class_settings_reader indata(inputfile);
    if(indata.found_file){
        settings::load_from_file(indata);
    }else{
        auto hdf5 = std::make_shared<class_hdf5_file> (outputfile,"",false,false,true);
        if (hdf5->fileMode == class_hdf5_file::FileMode::OPEN) {
            spdlog::info("The output file existed already: {}", hdf5->get_file_name());
            spdlog::info("Loading settings from existing file.");
            settings::load_from_hdf5(*hdf5);
        }else{
            spdlog::warn("Couldn't find an inputfile or previous outputfile to load settings. Running defaults.");
        }

    }


//    // Set console settings
//    if (settings::console::verbosity < 0 or settings::console::verbosity > 6){
//        std::cerr << "ERROR: Expected verbosity level integer in [0-6]. Got: " << settings::console::verbosity << std::endl;
//        exit(2);
//    }else{
//        spdlog::level::level_enum lvl = static_cast<spdlog::level::level_enum>(settings::console::verbosity);
//        spdlog::set_level(lvl);
//        spdlog::debug("Verbosity level: {}", spdlog::level::to_string_view(lvl));
//    }
//    if ( not settings::console::timestamp){
//        spdlog::set_pattern("[%n]%^[%=8l]%$ %v");
//    }




//
//
//
//    spdlog = spdlog::get("DMRG++");
//
//    spdlog->warn("TESTING MY LOGGER");
    //Initialize the algorithm class
    //This class stores simulation data automatically to a file specified in the input file
    class_algorithm_launcher launcher;

    //Run the algorithms
    launcher.run_algorithms();


    return 0;
}



