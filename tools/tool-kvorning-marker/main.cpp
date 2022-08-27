/*! \file */

#include "algorithms/AlgorithmLauncher.h"
#include "config/settings.h"
#include "io/filesystem.h"
#include "math/rnd.h"
#include "math/tenx.h"
#include "tools/common/log.h"
#include <h5pp/h5pp.h>

#if defined(OPENBLAS_AVAILABLE)
    #include <openblas/cblas.h>
    #include <openblas/openblas_config.h>
#endif

#if defined(MKL_AVAILABLE)
    #include <mkl.h>
    #include <mkl_service.h>
#endif

#include "conf.h"
#include "config/parse.h"
#include "debug/info.h"
#include "debug/stacktrace.h"
#include "env/environment.h"
#include <thread>

void clean_up() {
    if(not settings::storage::use_temp_dir) return;
    if(fs::exists(settings::storage::tmp::hdf5_temp_path)) {
        try {
            tools::log->info("Cleaning up temporary file: [{}]", settings::storage::tmp::hdf5_temp_path);
            h5pp::hdf5::moveFile(settings::storage::tmp::hdf5_temp_path, settings::storage::tmp::hdf5_final_path, h5pp::FilePermission::REPLACE);
        } catch(const std::exception &err) { tools::log->info("Cleaning not needed: {}", err.what()); }
    }
    H5garbage_collect();
    H5Eprint(H5E_DEFAULT, stderr);
}

// MWE: https://godbolt.org/z/jddxod53d
int parse(int argc, char **argv) {
    using namespace settings;
    using namespace h5pp;
    using namespace spdlog;

    auto s2e_log     = mapStr2Enum<spdlog::level::level_enum>("trace", "debug", "info");
    auto s2e_logh5pp = mapStr2Enum<h5pp::LogLevel>("trace", "debug", "info");
    auto s2e_model   = ModelType_s2e;
    int  dummy       = 0;
    bool noseedname  = false;

    auto preload = [&argc, &argv, &s2e_log]() -> int {
        CLI::App pre;
        pre.get_formatter()->column_width(90);
        pre.option_defaults()->always_capture_default();
        pre.allow_extras(true);
        pre.set_help_flag("--help-preload", "Help for preloading configuration");
        /* clang-format off */
        pre.add_option("-c,--config"                       , input::config_filename , "Path to a .cfg or .h5 file from a previous simulation");
        pre.add_option("-v,--log,--verbosity,--loglevel"   , console::loglevel      , "Log level of DMRG++")->transform(CLI::CheckedTransformer(s2e_log, CLI::ignore_case))->type_name("ENUM");
        pre.add_option("--timestamp"                       , console::timestamp     , "Log timestamp");
        /* clang-format on */
        pre.parse(argc, argv);
        tools::log = tools::Logger::setLogger("DMRG++ config", settings::console::loglevel, settings::console::timestamp);
        tools::log->info("Preloading {}", input::config_filename);
        //  Try loading the given config file.
        //  Note that there is a default "input/input.config" if none was given
        Loader dmrg_config(settings::input::config_filename);
        if(dmrg_config.file_exists) {
            dmrg_config.load();
            settings::load(dmrg_config);
        } else if(pre.get_option("--config")->empty()) {
            tools::log->warn("The default config file does not exist: {}", input::config_filename);
        } else
            throw except::runtime_error("Could not find config file: {}", settings::input::config_filename); // Invalid file
        return 0;
    };
    preload();

    CLI::App app;
    app.description("DMRG++: An MPS-based algorithm to calculate 1D quantum-states");
    app.get_formatter()->column_width(90);
    app.option_defaults()->always_capture_default();
    app.allow_extras(false);
    /* clang-format off */
    app.add_flag("--help-preload"                      , "Print help related to preloading configuration");
    app.add_option("-c,--config"                       , input::config_filename         , "Path to a .cfg or .h5 file from a previous simulation");
    app.add_option("-m,--model"                        , model::model_type              , "Select the Hamiltonian")->transform(CLI::CheckedTransformer(s2e_model, CLI::ignore_case));
    app.add_option("-b,--bitfield"                     , input::bitfield                , "Integer whose bitfield sets the initial product state. Negative is unused");
    app.add_option("-n,--stlthreads"                   , threading::stl_threads         , "Number of C++11 threads (Used by Eigen::Tensor)");
    app.add_option("-o,--outfile"                      , storage::output_filepath       , "Path to the output file. The seed number gets appended by default (see -x)");
    app.add_option("-s,--seed"                         , input::seed                    , "Positive number seeds the random number generator");
    app.add_option("-t,--ompthreads"                   , threading::omp_threads         , "Number of OpenMP threads");
    app.add_flag  ("-x,--noseedname"                   , noseedname                     , "Do not append seed to the output filename");
    app.add_option("-z,--compression"                  , storage::compression_level     , "Compression level of h5pp")->check(CLI::Range(0,9));
    app.add_flag  ("-r,--resume"                                                        , "Resume simulation from last iteration");
    app.add_option("--resume-iter"                     , storage::file_resume_iter      , "Resume from iteration");
    app.add_option("--resume-name"                     , storage::file_resume_name      , "Resume from state matching this name");
    app.add_option("-v,--log,--verbosity,--loglevel"   , console::loglevel              , "Log level of DMRG++")->transform(CLI::CheckedTransformer(s2e_log, CLI::ignore_case))->type_name("ENUM");
    app.add_option("-V,--logh5pp"                      , console::logh5pp               , "Log level of h5pp")->transform(CLI::CheckedTransformer(s2e_logh5pp, CLI::ignore_case))->type_name("ENUM");
    app.add_option("--timestamp"                       , console::timestamp             , "Log timestamp");
    app.add_option("--dummyrange"                      , dummy                          , "Dummy")->check(CLI::Range(0,3));
    /* clang-format on */

    app.parse(argc, argv);

    if(app.count("--resume") > 0 or app.count("--resume-iter") > 0 or app.count("--resume-name") > 0) {
        tools::log->info("Resuming from iter {}", storage::file_resume_iter);
        settings::storage::file_collision_policy = FileCollisionPolicy::RESUME;
    }

    //    for(const auto &res : app.get_options()) fmt::print("{:<32} = {}\n", res->get_name(), res->results());

    // Generate the correct output filename based on given seeds
    if(not noseedname) {
        settings::storage::output_filepath = filename_append_number(settings::storage::output_filepath, settings::input::seed);
        settings::storage::output_filepath = filename_append_number(settings::storage::output_filepath, settings::input::bitfield);
    }

    return 0;
}

/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/
int main(int argc, char *argv[]) {
    settings::parse(argc, argv);

    tools::log = tools::Logger::setLogger("DMRG++ main", settings::console::loglevel, settings::console::timestamp);

    // print environment and git status
    tools::log->info("Hostname        : {}", debug::hostname());
    tools::log->info("Build hostname  : {}", env::build::hostname);
    tools::log->info("Git branch      : {}", env::git::branch);
    tools::log->info("    commit hash : {}", env::git::commit_hash);
    tools::log->info("    revision    : {}", env::git::revision);

    // Register termination codes and what to do in those cases
    debug::register_callbacks();

    // Make sure to move the file back from temp location
    std::atexit(clean_up);
    std::at_quick_exit(clean_up);

    // Seed with random::device initially (This also takes care of srand used by Eigen)
    // This is to make reproducible simulations
    rnd::seed(settings::input::seed);

// Set the number of threads to be used
#if defined(EIGEN_USE_THREADS)
    if(settings::threading::stl_threads <= 0) { settings::threading::stl_threads = (int) std::thread::hardware_concurrency(); }
    tenx::omp::setNumThreads(settings::threading::stl_threads);
    tools::log->info("Eigen3 Tensor | stl threads {}", tenx::omp::num_threads);
#else
    if(settings::threading::stl_threads > 1)
        tools::log->warn("EIGEN_USE_THREADS is not defined: "
                         "Failed to enable threading in Eigen::Tensor with stl_threads = {}",
                         settings::threading::stl_threads);
#endif

#if defined(_OPENMP)
    if(settings::threading::omp_threads <= 0) { settings::threading::omp_threads = (int) std::thread::hardware_concurrency(); }
    omp_set_num_threads(settings::threading::omp_threads); // Should only need this. Both Eigen (non-Tensor) and MKL listen to this
    tools::log->info("OpenMP | threads {} | max active levels {}", omp_get_max_threads(), omp_get_max_active_levels());
#endif

#if defined(OPENBLAS_AVAILABLE)
    auto        openblas_parallel_mode = openblas_get_parallel();
    std::string openblas_parallel_str;
    if(openblas_parallel_mode == 0) openblas_parallel_str = "seq";
    if(openblas_parallel_mode == 1) openblas_parallel_str = "threads";
    if(openblas_parallel_mode == 2) openblas_parallel_str = "openmp";
    if(openblas_parallel_mode == 1) openblas_set_num_threads(settings::threading::omp_threads); // Use the omp_threads level for blas-related threading
    tools::log->info("{} threads {} | parallel mode [{}:{}] | core type {} | config {} | multithread threshold {}", OPENBLAS_VERSION,
                     openblas_get_num_threads(), openblas_parallel_mode, openblas_parallel_str, openblas_get_corename(), openblas_get_config(),
                     OPENBLAS_GEMM_MULTITHREAD_THRESHOLD);
#endif
#if defined(MKL_AVAILABLE)
    tools::log->info("Intel MKL | threads {}", mkl_get_max_threads());
#endif

#if defined(EIGEN_USE_MKL_ALL)
    tools::log->info("Eigen3 | threads {} | EIGEN_USE_MKL_ALL", Eigen::nbThreads());
#elif defined(EIGEN_USE_BLAS)
    tools::log->info("Eigen3 | threads {} | EIGEN_USE_BLAS", Eigen::nbThreads());
#else
    tools::log->info("Eigen3 | threads {} | No BLAS backend", Eigen::nbThreads());
#endif

    // Initialize the algorithm class
    // This class stores simulation data automatically to a file specified in the config file
    AlgorithmLauncher launcher;

    // Run the algorithms
    launcher.run_algorithms();

    return 0;
}
