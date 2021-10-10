/*! \file */

#include <algorithms/AlgorithmLauncher.h>
#include <config/settings.h>
#include <h5pp/h5pp.h>
#include <io/filesystem.h>
#include <math/rnd.h>
#include <math/tenx.h>
#include <tools/common/log.h>

#if defined(OPENBLAS_AVAILABLE)
    #include <openblas/cblas.h>
    #include <openblas/openblas_config.h>
#endif

#if defined(MKL_AVAILABLE)
    #define MKL_Complex8  std::complex<float>
    #define MKL_Complex16 std::complex<double>
    #include <mkl.h>
    #include <mkl_service.h>
#endif
#include <CLI/CLI.hpp>
#include <config/loader.h>
#include <debug/stacktrace.h>
#include <env/environment.h>
#include <thread>

void print_usage() {
    std::printf(
        R"(
==========  DMRG++  ============
Usage                       : DMRG++ [-option <value>].
-h                          : Help. Shows this text.
-b <positive integer>       : Integer whose bitfield sets the initial product state. Negative is unused (default -1)
-c <.cfg or .h5 filename>   : Full or relative path to a config file or hdf5 file from a previous simulation (which has a config file) (default = input.cfg)
-n <stl threads>            : Number of C++11 threads (Used by Eigen::Tensor)
-o <output filename base>   : Full or relative path to the output file. The seed number will be appended to this filename unless -x is passed.
-i <.cfg or .h5 filename>   : Full or relative path to a config file or hdf5 file from a previous simulation (which has a config file) (default = input.cfg)
-s <seed>                   : Positive number that seeds the random number generator (default = 1)
-t <omp threads>            : Number of OpenMP threads
-v                          : Enables trace-level verbosity
-x                          : Do not append seed to the output filename.

)");
}

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

    tools::log->info("Hostname        : {}", env::hostname);

    // print current Git status
    tools::log->info("Git branch      : {}", git::branch);
    tools::log->info("    commit hash : {}", git::commit_hash);
    tools::log->info("    revision    : {}", git::revision);
    bool noseedname = false;

    {
        using namespace settings;
        CLI::App app{"DMRG++"};
        app.description("An MPS-based algorithm to find 1D quantum-states");
        app.get_formatter()->column_width(60);
        /* clang-format off*/
        app.option_defaults()->always_capture_default();
        app.add_option("-b,--bitfield", input::bitfield, "Integer whose bitfield sets the initial product state. Negative is unused")
            ->default_val(input::bitfield);
        app.add_option("-c,--config", input::config_filename, "Path to a .cfg or .h5 file from a previous simulation")->default_val(input::config_filename);
        app.add_option("-n,--stlthreads", threading::stl_threads, "Number of C++11 threads (Used by Eigen::Tensor)")->default_val(threading::stl_threads);
        app.add_option("-o,--outfile", storage::output_filepath, "Path to the output file. The seed number gets appended by default (see -x)");
        app.add_option("-s,--seed", input::seed, "Positive number seeds the random number generator");
        app.add_option("-t,--ompthreads", threading::omp_threads, "Number of OpenMP threads");
        app.add_option("-v,--verbose", console::verbosity, "Sets verbosity level");
        app.add_flag("-x,--noseedname", noseedname, "Do not append seed to the output filename");

        /* clang-format on */
        CLI11_PARSE(app, argc, argv);
        fmt::format("exe: {}", h5pp::fs::read_symlink("/proc/self/exe"));
    }
    exit(0);

    //    cxxopts::Options options("DMRG++", "An MPS-based algorithm to find 1D quantum-states");
    //    /* clang-format off */
    //    options.add_options()
    //    ("h,help",        "Show help")
    //    ("b,bitfield",    "Integer whose bitfield sets the initial product state. Negative is unused", cxxopts::value<long>())
    //    ("c,config",      "Path to a .cfg or .h5 file from a previous simulation", cxxopts::value<std::string>()->default_value("input/input.cfg"))
    //    ("n,stlthreads",  "Number of C++11 threads (Used by Eigen::Tensor)",                           cxxopts::value<int>())
    //    ("o,outfile",     "Path to the output file. The seed number gets appended by default (see
    //    -x)",cxxopts::value<std::string>()->default_value("output/output.h5"))
    //    ("s,seed",        "Positive number seeds the random number generator",                         cxxopts::value<long>())
    //    ("t,ompthreads",  "Number of OpenMP threads",                                                  cxxopts::value<int>())
    //    ("v,verbose",     "Sets verbosity level",                                                      cxxopts::value<size_t>())
    //    ("x,noseedname",  "Do not append seed to the output filename",                                 cxxopts::value<bool>());

    //    auto in = options.parse(argc, argv);
    //    if(in["help"].count() > 0) {
    //        fmt::print(options.help());
    //        exit(0);
    //    }
    //    if(in["config"].count() > 0)    settings::input::config_filename = in["config"].as<std::string>();

    //  Try loading the given config file.
    //  Note that there is a default "input/input.config" if none was given
    Loader dmrg_config(settings::input::config_filename);
    if(dmrg_config.file_exists) {
        dmrg_config.load();
        settings::load(dmrg_config); // B2
    } else
        throw std::runtime_error(fmt::format("Could not find config file: {}", settings::input::config_filename)); // Invalid file
                                                                                                                   //
                                                                                                                   //    // Override the other settings
                                                                                                                   //    bool        noseedname = false;
    //    if(in["bitfield"].count() > 0)      settings::input::bitfield           = in["bitfield"].as<long>();
    //    if(in["stlthreads"].count() > 0)    settings::threading::stl_threads    = in["stlthreads"].as<int>();
    //    if(in["outfile"].count() > 0)       settings::storage::output_filepath  = in["outfile"].as<std::string>();
    //    if(in["seed"].count() > 0)          settings::input::seed               = in["seed"].as<long>();
    //    if(in["ompthreads"].count() > 0)    settings::threading::omp_threads    = in["ompthreads"].as<int>();
    //    if(in["verbose"].count() > 0)       settings::console::verbosity        = in["verbose"].as<size_t>();
    //    if(in["noseedname"].count() > 0)    noseedname                          = in["noseedname"].as<bool>();

    tools::log = tools::Logger::setLogger("DMRG++ main", settings::console::verbosity, settings::console::timestamp);
    /* clang-format on */

    // C: Generate the correct output filename based on given seeds
    if(not noseedname) {
        settings::storage::output_filepath = filename_append_number(settings::storage::output_filepath, settings::input::seed);
        settings::storage::output_filepath = filename_append_number(settings::storage::output_filepath, settings::input::bitfield);
    }

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
