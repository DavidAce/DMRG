/*! \file */


#include <simulation/nmspc_settings.h>
#include <math/nmspc_random.h>
#include <algorithms/class_algorithm_launcher.h>
#include <io/class_settings_reader.h>
#include <io/nmspc_logger.h>
#include <io/nmspc_filesystem.h>
#include <iostream>
#include <h5pp/h5pp.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <general/nmspc_omp.h>


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
-r <seed>                   : Positive number that seeds the random number generator (default = 1)
-s <state number>           : Positive number whose bitfield sets the initial product state. Negative is unused (default -1)
-t <num threads>            : Number of OpenMP threads
-o <output file>            : Full or relative path to the output file (output)
-x                          : Do not append seed to the output filename.

)";

}


#include <signal.h>

void signal_callback_handler(int signum) {
    switch(signum){
        case SIGTERM:  {std::cout  << "Caught SIGTERM" <<std::endl; break;}
        case SIGKILL:  {std::cout  << "Caught SIGKILL" <<std::endl; break;}
        case SIGINT :  {std::cout  << "Caught SIGINT"  <<std::endl; break;}
        case SIGHUP :  {std::cout  << "Caught SIGHUP"  <<std::endl; break;}
        case SIGQUIT : {std::cout  << "Caught SIGQUIT" <<std::endl; break;}
    }
    std::cout << "Exiting" << std::endl << std::flush;
    std::quick_exit(signum);
}




/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main(int argc, char* argv[]) {
    //Register termination codes and what to do in those cases
    signal(SIGTERM , signal_callback_handler);
    signal(SIGINT  , signal_callback_handler);
    signal(SIGKILL , signal_callback_handler);
    signal(SIGHUP  , signal_callback_handler);
    signal(SIGQUIT , signal_callback_handler);

    auto log = Logger::setLogger("DMRG",0);
    using namespace tools;
    // print current Git status
    log->info("Git branch      : {}",GIT::BRANCH);
    log->info("    commit hash : {}",GIT::COMMIT_HASH);
    log->info("    revision    : {}",GIT::REVISION);



    bool append_seed = true;
    bool load_previous = false;
    while(true){
        char opt = getopt(argc, argv, "hi:lr:s:t:o:x");
        if (opt == EOF) break;
        if(optarg == nullptr) log->info("Parsing input argument: -{}",opt);
        else                  log->info("Parsing input argument: -{} {}",opt,optarg);

        switch(opt){
            case 'i': {
                settings::input::input_filename = std::string(optarg);
                class_settings_reader indata(settings::input::input_filename);
                if(indata.found_file){
                    settings::load_from_file(indata);
                }else{
                    log->critical("Could not find input file: {}", settings::input::input_filename);
                    exit(1);
                }
                continue;
            }
            case 'l': load_previous = true; continue;
            case 'r': {
                long seed = std::strtol(optarg,nullptr,10);
                if(seed >= 0){
                    log->info("Replacing model::seed_model {} -> {}", settings::model::seed, seed);
                    settings::model::seed = seed;
                }
                continue;
            }
            case 's': {
                long state_number =  std::strtol(optarg,nullptr,10);
                if(state_number >= 0) {
                    log->info("Replacing model::seed_state {} -> {}", settings::model::state_number,state_number);
                    settings::model::state_number = state_number;
                }
                continue;
            }
            case 't': {
                int num_threads = std::strtol(optarg,nullptr,10);
                if(num_threads > 0) {
                    log->info("Setting OpenMP threads to {}", num_threads);
                    settings::threading::num_threads_eigen = num_threads;
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

    if (not load_previous and append_seed and settings::model::seed >= 0 ){
        //Append the seed_model to the output filename
        fs::path oldFileName = settings::output::output_filename;
        fs::path newFileName = settings::output::output_filename;
        newFileName.replace_filename(oldFileName.stem().string() + "_" + std::to_string(settings::model::seed) + oldFileName.extension().string() );
        settings::output::output_filename = newFileName.string();
        log->info("Appended model::seed_model to output filename: [{}] --> [{}]",oldFileName.string(), newFileName.string());
    }
    if (not load_previous and append_seed and settings::model::state_number >= 0){
        //Append the seed_state to the output filename
        fs::path oldFileName = settings::output::output_filename;
        fs::path newFileName = settings::output::output_filename;
        newFileName.replace_filename(oldFileName.stem().string() + "_" + std::to_string(settings::model::state_number) + oldFileName.extension().string() );
        settings::output::output_filename = newFileName.string();
        log->info("Appended model::seed_state to output filename: [{}] --> [{}]",oldFileName.string(), newFileName.string());
    }


    //Set the number of threads to be used


    #ifdef _OPENMP
        if(settings::threading::num_threads_omp   <= 0) { settings::threading::num_threads_omp   = (int)std::thread::hardware_concurrency(); }
        if(settings::threading::num_threads_eigen <= 0) { settings::threading::num_threads_eigen = (int)std::thread::hardware_concurrency(); }
        if(settings::threading::num_threads_blas  <= 0) { settings::threading::num_threads_blas  = (int)std::thread::hardware_concurrency(); }

        omp_set_num_threads(settings::threading::num_threads_omp);
        Eigen::setNbThreads(settings::threading::num_threads_eigen);
        Eigen::initParallel();
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
            mkl_set_num_threads(settings::threading::num_threads_blas);
            log->info("Using Intel MKL with {} threads", mkl_get_max_threads());
        #endif

    #endif










    // Seed with random::device initially (This also takes care of srand used by Eigen)
    // Note that we re-seed if a seed_model or seed_state is given also.
    // This is to make reproducible simulations
    rn::seed(settings::model::seed);


    //Initialize the algorithm class
    //This class stores simulation data automatically to a file specified in the input file
    class_algorithm_launcher launcher;

    //Run the algorithms
    launcher.run_algorithms();

    return 0;
}



