/*! \file */

#include <algorithms/class_algorithm_launcher.h>
#include <h5pp/h5pp.h>
#include <io/class_config_reader.h>
#include <io/nmspc_filesystem.h>
#include <io/nmspc_logger.h>
#include <iostream>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <general/nmspc_tensor_omp.h>

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
-c <existing .cfg filename> : Full or relative path of the config file that sets up the simulation (default = input.cfg)
-l <existing hdf5 filename> : Load from the given file (i.e., from a previous simulation).
-r <seed>                   : Positive number that seeds the random number generator (default = 1)
-s <state number>           : Positive number whose bitfield sets the initial product state. Negative is unused (default -1)
-t <num threads>            : Number of OpenMP threads
-o <output filename>        : Full or relative path to the output file (output). The seed number will be appended to this filename unless -x is passed.
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

std::string filename_append_number(const std::string_view filename,const long number){
    if(number < 0) return std::string(filename);
    auto log = Logger::setLogger("DMRG",0);
    //Append the seed_model to the output filename
    h5pp::fs::path oldFileName = filename;
    h5pp::fs::path newFileName = filename;
    newFileName.replace_filename(oldFileName.stem().string() + "_" + std::to_string(number) + oldFileName.extension().string() );
    log->info("Appended number [{}] to filename: [{}]",number, newFileName.string());
    return newFileName.string();
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
    std::string config_filename;
    std::string h5_load_filename;
    while(true){
        char opt = getopt(argc, argv, "hc:l:r:s:t:o:x");
        if (opt == EOF) break;
        if(optarg == nullptr) log->info("Parsing input argument: -{}",opt);
        else                  log->info("Parsing input argument: -{} {}",opt,optarg);

        switch(opt){
            case 'c': config_filename = std::string(optarg); continue;
            case 'l': h5_load_filename = std::string(optarg); continue;
            case 'r': {
                long seed = std::strtol(optarg,nullptr,10);
                if(seed >= 0){
                    log->info("Replacing model::seed {} -> {}", settings::model::seed, seed);
                    settings::model::seed = seed;
                }
                continue;
            }
            case 's': {
                long state_number =  std::strtol(optarg,nullptr,10);
                if(state_number >= 0) {
                    log->info("Replacing model::state_number {} -> {}", settings::model::state_number,state_number);
                    settings::model::state_number = state_number;
                }
                continue;
            }
            case 't': {
                int num_threads = std::strtol(optarg,nullptr,10);
                if(num_threads > 0) {
                    log->info("Setting OpenMP threads to {}", num_threads);
                    settings::threading::num_threads = num_threads;
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


    /*! There are some ways this could go. This decision tree goes in order from A:
        A1: hdf5 file not given: go to B
        A2: hdf5 file given: Load hdf5 file config -> go to B

        B1: config file not given: go to C if coming from A2, else use the default
        B2: config file given: load given config -> go to C

        C: generate output filename -> go to D

        D1: output file exists: check option in settings::output::create_mode:
            - case OPEN:        Load config and simulation state -> continue simulation
            - case RENAME:      Rename output file -> go to D2
            - case TRUNCATE:    Delete output file -> go to D2

        D2: output file does not exist: Check if coming from A2
            - Yes: load simulation state from given hdf5 file -> continue simulation
            - No: start new simulation

        NOTE that D is done separately inside of each algorithm (iDMRG, iTEBD... etc)
    */

    // A: Try loading config from given hdf5 file
    if(not h5_load_filename.empty()){
        settings::input::h5_load_filename = h5_load_filename;
        // A2: If we were given a previous simulation file it must be loaded first
        h5pp::File h5file (settings::input::h5_load_filename, h5pp::AccessMode::READONLY,h5pp::CreateMode::OPEN);
        settings::load_config_from_hdf5(h5file);
    }

    //B: Try loading given config file.
    //   Note that there is a default "input/input.cfg" if none was given
    if(not config_filename.empty()){
        settings::input::config_filename = config_filename;
        class_config_reader indata(settings::input::config_filename);
        if(indata.found_file){
            // B2
            settings::load_config_from_cfg(indata);
        }else{
            // Invalid file given
            throw std::runtime_error(fmt::format("Could not find config file: {}", settings::input::config_filename));
        }
    }//else use default cfg


    // C: Generate the correct output filename based on given seeds
    if (append_seed) {
        settings::output::output_filename = filename_append_number(settings::output::output_filename, settings::model::seed);
        settings::output::output_filename = filename_append_number(settings::output::output_filename, settings::model::state_number);
    }

    // D
    if(h5pp::fs::exists(settings::output::output_filename)){
        switch(settings::output::create_mode){
            case h5pp::CreateMode::OPEN : {

            }
            case h5pp::CreateMode::RENAME : {

            }
            case h5pp::CreateMode::TRUNCATE : {

            }
        }

    }else{

    }

    //Set the number of threads to be used

    #ifdef _OPENMP
        if(settings::threading::num_threads <= 0) { settings::threading::num_threads = (int)std::thread::hardware_concurrency(); }

        omp_set_num_threads(settings::threading::num_threads);
        Eigen::setNbThreads(settings::threading::num_threads);
        log->info("Using Eigen  with {} threads",Eigen::nbThreads());
        log->info("Using OpenMP with {} threads",omp_get_max_threads());

        #ifdef OpenBLAS_AVAILABLE
                openblas_set_num_threads(settings::threading::num_threads);
                std::cout << OPENBLAS_VERSION
                          << " compiled with parallel mode " << openblas_get_parallel()
                          << " for target " << openblas_get_corename()
                          << " with config " << openblas_get_config()
                          << " with multithread threshold " << OPENBLAS_GEMM_MULTITHREAD_THRESHOLD
                          << ". Running with " << openblas_get_num_threads() << " thread(s)" << std::endl;
        #endif

        #ifdef MKL_AVAILABLE
            mkl_set_num_threads(settings::threading::num_threads);
            log->info("Using Intel MKL with {} threads", mkl_get_max_threads());
        #endif

    #endif


    // Seed with random::device initially (This also takes care of srand used by Eigen)
    // Note that we re-seed if a seed_model or seed_state is given also.
    // This is to make reproducible simulations
    rn::seed(settings::model::seed);


    //Initialize the algorithm class
    //This class stores simulation data automatically to a file specified in the config file
    class_algorithm_launcher launcher;

    //Run the algorithms
    launcher.run_algorithms();

    return 0;
}



