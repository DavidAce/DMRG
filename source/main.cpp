/*! \file */

#include <algorithms/AlgorithmLauncher.h>
#include <config/settings.h>
#include <h5pp/h5pp.h>
#include <io/filesystem.h>
#include <math/rnd.h>
#include <math/tenx.h>
#include <tools/common/log.h>

#if defined(OPENBLAS_AVAILABLE)
    #include <cblas.h>
    #include <openblas_config.h>
#endif

#if defined(MKL_AVAILABLE)
    #define MKL_Complex8  std::complex<float>
    #define MKL_Complex16 std::complex<double>
//    #include <mkl.h>
    #include <mkl_service.h>
#endif

#include <config/loader.h>
#include <debug/stacktrace.h>
#include <getopt.h>
#include <gitversion.h>
#include <thread>

#if __has_include(<unistd.h>)
    #include <unistd.h>
#endif

void print_usage() {
    std::printf(
        R"(
==========  DMRG++  ============
Usage                       : DMRG++ [-option <value>].
-h                          : Help. Shows this text.
-b <positive integer>       : Integer whose bitfield sets the initial product state. Negative is unused (default -1)
-c <.cfg or .h5 filename>   : Full or relative path to a config file or hdf5 file from a previous simulation (which has a config file) (default = input.cfg)
-n <stl threads>            : Number of C++11 threads (Used by Eigen::Tensor)
-o <output filename base>   : Full or relative path to the output file (output). The seed number will be appended to this filename unless -x is passed.
-i <.cfg or .h5 filename>   : Full or relative path to a config file or hdf5 file from a previous simulation (which has a config file) (default = input.cfg)
-s <seed>                   : Positive number that seeds the random number generator (default = 1)
-t <omp threads>            : Number of OpenMP threads
-v                          : Enables trace-level verbosity
-x                          : Do not append seed to the output filename.

)");
}

void clean_up() {
    if(not settings::output::use_temp_dir) return;
    if(fs::exists(settings::output::tmp::hdf5_temp_path)) {
        try {
            tools::log->info("Cleaning up temporary file: [{}]", settings::output::tmp::hdf5_temp_path);
            h5pp::hdf5::moveFile(settings::output::tmp::hdf5_temp_path, settings::output::tmp::hdf5_final_path, h5pp::FilePermission::REPLACE);
        } catch(const std::exception &err) { tools::log->info("Cleaning not needed: {}", err.what()); }
    }
    H5garbage_collect();
    H5Eprint(H5E_DEFAULT, stderr);
}

std::string filename_append_number(const std::string &filename, const long number) {
    if(number < 0) return filename;
    // Append the seed_model to the output filename
    h5pp::fs::path oldFileName = filename;
    h5pp::fs::path newFileName = filename;
    if(oldFileName.stem().string().find(std::to_string(number)) != std::string::npos) return filename;
    newFileName.replace_filename(fmt::format("{}_{}{}", oldFileName.stem().string(), number, oldFileName.extension().string()));
    tools::log->info("Appended number [{}] to filename: [{}]", number, newFileName.string());
    return newFileName.string();
}

/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main(int argc, char *argv[]) {
    // Register termination codes and what to do in those cases
    debug::register_callbacks();

    // Make sure to move the file back from temp location
    std::atexit(clean_up);
    std::at_quick_exit(clean_up);

    tools::log = tools::Logger::setLogger("DMRG++ main", 0, true);
    using namespace tools;

#if __has_include(<unistd.h>)
    char name[HOST_NAME_MAX];
    auto err = gethostname(name, HOST_NAME_MAX);
    if(err == 0) tools::log->info("Hostname        : {}", name);
#endif

    // print current Git status
    tools::log->info("Git branch      : {}", GIT::BRANCH);
    tools::log->info("    commit hash : {}", GIT::COMMIT_HASH);
    tools::log->info("    revision    : {}", GIT::REVISION);

    // Here we use getopt to parse CLI input
    // Note that CLI input always override config-file values
    // wherever they are found (config file, h5 file)

    bool        append_seed = true;
    std::string config;
    std::string output;
    long        verbosity   = -1;
    long        seed        = -1;
    long        bitfield    = -1;
    long        omp_threads = -1;
    long        stl_threads = -1;

    while(true) {
        char opt = static_cast<char>(getopt(argc, argv, "hb:c:i:n:o:s:t:vx"));
        if(opt == EOF) break;
        if(optarg == nullptr)
            tools::log->info("Parsing input argument: -{}", opt);
        else
            tools::log->info("Parsing input argument: -{} {}", opt, optarg);
        switch(opt) {
            case 'b': bitfield = std::strtol(optarg, nullptr, 10); continue;
            case 'c':
            case 'i': config = std::string(optarg); continue;
            case 'n': stl_threads = std::strtol(optarg, nullptr, 10); continue;
            case 'o': output = std::string(optarg); continue;
            case 's': seed = std::strtol(optarg, nullptr, 10); continue;
            case 't': omp_threads = std::strtol(optarg, nullptr, 10); continue;
            case 'v': verbosity = 0; continue;
            case 'x': append_seed = false; continue;
            case ':': tools::log->error("Option -{} needs a value", opt); break;
            case 'h':
            case '?':
            default: print_usage(); exit(0);
            case -1: break;
        }
        break;
    }
    /*
     There may be multiple config files to consider:
         1) Given from CLI (.config/.h5)
         2) Inside the output file "<output>_<seed>.h5" if it already exists.

     What should one do?
     Simplest solution: Always ignore case 2)!
     Taking 2) into account leads to very confusing policies.

     What are the implications?
        - If CLI passes a .config file we have to override some of its settings:
            - output

    */

    /*! It's important that we do things in this order:
        A1: config file not given:  use the default input/input.config
        A2: config file given with .config/h5 extension: load given config

        B: Override settings with parameters given through CLI

        C: generate output filename. If the seed is already on the filename, it is not appended again.
    */

    // B: Try loading given config file.
    //   Note that there is a default "input/input.config" if none was given
    if(not config.empty()) {
        Loader dmrg_config(config);
        if(dmrg_config.file_exists) {
            dmrg_config.load();
            settings::load(dmrg_config); // B2
        } else
            throw std::runtime_error(fmt::format("Could not find config file: {}", config)); // Invalid file given
        settings::input::config_filename = config;
    } // else use default config

    // B: Override settings
    if(seed >= 0) settings::input::seed = seed;
    if(bitfield >= 0) settings::input::bitfield = bitfield;
    if(stl_threads >= 0) settings::threading::stl_threads = static_cast<int>(stl_threads);
    if(omp_threads >= 0) settings::threading::omp_threads = static_cast<int>(omp_threads);
    if(not output.empty()) settings::output::output_filepath = output;
    if(verbosity >= 0) settings::console::verbosity = static_cast<size_t>(verbosity);
    tools::log = tools::Logger::setLogger("DMRG++ main", settings::console::verbosity, settings::console::timestamp);

    // C: Generate the correct output filename based on given seeds
    if(append_seed) {
        settings::output::output_filepath = filename_append_number(settings::output::output_filepath, settings::input::seed);
        settings::output::output_filepath = filename_append_number(settings::output::output_filepath, settings::input::bitfield);
    }

    // Seed with random::device initially (This also takes care of srand used by Eigen)
    // This is to make reproducible simulations
    rnd::seed(settings::input::seed);

// Set the number of threads to be used
#if defined(EIGEN_USE_THREADS)
    if(settings::threading::stl_threads <= 0) { settings::threading::stl_threads = (int) std::thread::hardware_concurrency(); }
    tenx::omp::setNumThreads(settings::threading::stl_threads);
    tools::log->info("Using Eigen Tensor with {} C++11 threads", tenx::omp::num_threads);
#else
    if(settings::threading::stl_threads > 1)
        tools::log->warn("EIGEN_USE_THREADS is not defined: "
                         "Failed to enable threading in Eigen::Tensor with stl_threads = {}",
                         settings::threading::stl_threads);
#endif

#if defined(_OPENMP)
    if(settings::threading::omp_threads <= 0) { settings::threading::omp_threads = (int) std::thread::hardware_concurrency(); }
    omp_set_num_threads(settings::threading::omp_threads); // Should only need this. Both Eigen (non-Tensor) and MKL listen to this
                                                           //    omp_set_max_active_levels(1);
    tools::log->info("Using OpenMP with {} threads and {} active levels", omp_get_max_threads(), omp_get_max_active_levels());
#endif
#if defined(OPENBLAS_AVAILABLE)
    auto        openblas_parallel_mode = openblas_get_parallel();
    std::string openblas_parallel_str;
    if(openblas_parallel_mode == 0) openblas_parallel_str = "seq";
    if(openblas_parallel_mode == 1) openblas_parallel_str = "threads";
    if(openblas_parallel_mode == 2) openblas_parallel_str = "openmp";
    if(openblas_parallel_mode == 1) openblas_set_num_threads(settings::threading::omp_threads); // Use the omp_threads level for blas-related threading
    tools::log->info("{} compiled with parallel mode [{}] for target {} with config {} | multithread threshold {} | running with {} threads", OPENBLAS_VERSION,
                     openblas_parallel_str, openblas_get_corename(), openblas_get_config(), OPENBLAS_GEMM_MULTITHREAD_THRESHOLD, openblas_get_num_threads());
#endif
#if defined(MKL_AVAILABLE)
    tools::log->info("Using Intel MKL with {} threads", mkl_get_max_threads());
#endif

    // Initialize the algorithm class
    // This class stores simulation data automatically to a file specified in the config file
    AlgorithmLauncher launcher;

    // Run the algorithms
    launcher.run_algorithms();

    return 0;
}
