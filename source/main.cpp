/*! \file */


#include <simulation/nmspc_settings.h>
#include <general/nmspc_random_numbers.h>
#include <algorithms/class_algorithm_launcher.h>
#include <io/class_settings_reader.h>
#include <io/nmspc_logger.h>

#include <iostream>
#include <h5pp/h5pp.h>
#include <experimental/filesystem>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef OpenBLAS_AVAILABLE
#include <cblas.h>
#include <openblas_config.h>
#endif


#ifdef MKL_AVAILABLE
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl_service.h>
#include <mkl.h>
#endif

#include <thread>
//#include <unistd.h>
#include <getopt.h>
#include <gitversion.h>


void print_usage(){

std::cout<<
R"(
==========  DMRG++  ============
Usage                       : DMRG++ [-option <value>].
NOTE                        : Order of argument matters. In particular, set seeds AFTER the input file.
-h                          : Help. Shows this text.
-i <input file>             : Full or relative path of the input file that configures the simulation (default = input.cfg)
-l                          : Enables loading from the given output file (i.e., from a previous simulation)
-r <seed rng>               : Positive number that seeds the random number generator, used for model params (default = 1)
-s <seed state>             : Positive number that seed the initial state. Negative is unused (default -1)
-o <output file>            : Full or relative path to the output file (output)
-x                          : Do not append seed to the output filename.

)";

}



/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main(int argc, char* argv[]) {
    auto log = Logger::setLogger("DMRG",0);

    // print current Git status
    log->info("Git branch      : {}",GIT::BRANCH);
    log->info("    commit hash : {}",GIT::COMMIT_HASH);
    log->info("    revision    : {}",GIT::REVISION);



    bool append_seed = true;
    bool load_previous = false;
    while(true){
        char opt = getopt(argc, argv, "hi:lr:s:o:x");
        if (opt == EOF) break;
        if(optarg == nullptr) log->info("Parsing input argument: -{}",opt);
        else                  log->info("Parsing input argument: -{} {}",opt,optarg);

        switch(opt){
            case 'i': {
                settings::input::input_file = std::string(optarg);
                class_settings_reader indata(settings::input::input_file);
                if(indata.found_file){
                    settings::load_from_file(indata);
                }else{
                    log->critical("Could not find input file: {}", settings::input::input_file);
                    exit(1);
                }
                continue;
            }
            case 'l': load_previous = true; continue;
            case 'r': {
                int seed_init = (int) std::strtol(optarg,nullptr,10);
                if(seed_init >= 0){
                    log->info("Replacing model::seed_model {} -> {}", settings::model::seed_model, seed_init);
                    settings::model::seed_model = seed_init;

                }
                continue;
            }
            case 's': {
                int seed_state = (int) std::strtol(optarg,nullptr,10);
                if(seed_state >= 0) {
                    log->info("Replacing model::seed_state {} -> {}", settings::model::seed_state,seed_state);
                    settings::model::seed_state = seed_state;
                }
                continue;
            }
            case 'o': settings::output::output_filename = std::string(optarg); continue;
            case 'x': append_seed = false; continue;
            case ':': log->error("Option -{} needs a value", opt); break;
            case 'h':
            case '?':
            default: print_usage(); exit(0);
            case -1: break;
        }
        break;
    }



    if(load_previous){
        try{
            auto h5ppFile = std::make_shared<h5pp::File> (settings::output::output_filename, h5pp::AccessMode::READONLY, h5pp::CreateMode::OPEN);
            log->info("Loading settings from existing file [{}]", h5ppFile->getFilePath());
            settings::load_from_hdf5(*h5ppFile);
        }catch(std::exception &ex){
            log->info("Couldn't load from output file: {}", settings::output::output_filename, ex.what() );
            exit(0);
        }
    }

    if (not load_previous and append_seed and settings::model::seed_model >= 0 ){
        //Append the seed_model to the output filename
        namespace fs = std::experimental::filesystem;
        fs::path oldFileName = settings::output::output_filename;
        fs::path newFileName = settings::output::output_filename;
        newFileName.replace_filename(oldFileName.stem().string() + "_" + std::to_string(settings::model::seed_model) + oldFileName.extension().string() );
        settings::output::output_filename = newFileName.string();
        log->info("Appended model::seed_model to output filename: [{}] --> [{}]",oldFileName.string(), newFileName.string());
    }
    if (not load_previous and append_seed and settings::model::seed_state >= 0){
        //Append the seed_state to the output filename
        namespace fs = std::experimental::filesystem;
        fs::path oldFileName = settings::output::output_filename;
        fs::path newFileName = settings::output::output_filename;
        newFileName.replace_filename(oldFileName.stem().string() + "_" + std::to_string(settings::model::seed_state) + oldFileName.extension().string() );
        settings::output::output_filename = newFileName.string();
        log->info("Appended model::seed_state to output filename: [{}] --> [{}]",oldFileName.string(), newFileName.string());
    }


    //Set the number of threads to be used


    #ifdef _OPENMP
        if(settings::threading::num_threads_omp   <= 0) { settings::threading::num_threads_omp   = std::thread::hardware_concurrency(); }
        if(settings::threading::num_threads_eigen <= 0) { settings::threading::num_threads_eigen = std::thread::hardware_concurrency(); }
        if(settings::threading::num_threads_blas  <= 0) { settings::threading::num_threads_blas  = std::thread::hardware_concurrency(); }

        omp_set_num_threads(settings::threading::num_threads_omp);
        Eigen::setNbThreads(settings::threading::num_threads_eigen);
        log->info("Using Eigen  with {} threads",Eigen::nbThreads());
        log->info("Using OpenMP with {} threads",omp_get_max_threads());

        #ifdef OpenBLAS_AVAILABLE
                openblas_set_num_threads(settings::threading::num_threads_blas);
                std::cout << OPENBLAS_VERSION
                          << " compiled with parallel mode " << openblas_get_parallel()
                          << " for target " << openblas_get_corename()
                          << " with config " << openblas_get_config()
                          << " with multithread threshold " << OPENBLAS_GEMM_MULTITHREAD_THRESHOLD
                          << ". Running with " << openblas_get_num_threads() << " thread(s)" << std::endl;
        #endif

        #ifdef MKL_AVAILABLE
            if(settings::threading::num_threads_blas   <= 0){ settings::threading::num_threads_blas   = std::thread::hardware_concurrency(); }
            mkl_set_num_threads(settings::threading::num_threads_blas);
            log->info("Using Intel MKL with {} threads", mkl_get_max_threads());
        #endif

    #endif










    // Seed only this once (This also takes care of srand used by Eigen
    rn::seed(settings::model::seed_model);


    //Initialize the algorithm class
    //This class stores simulation data automatically to a file specified in the input file
    class_algorithm_launcher launcher;

    //Run the algorithms
    launcher.run_algorithms();

    return 0;
}



