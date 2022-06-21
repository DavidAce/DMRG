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
