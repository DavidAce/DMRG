#include "threading.h"
#include "math/tenx/threads.h"
#include "settings.h"
#include "tools/common/log.h"
#include <omp.h>

#if defined(OPENBLAS_AVAILABLE)
    #include <openblas/cblas.h>
    #include <openblas/openblas_config.h>
#endif

#if defined(MKL_AVAILABLE)
    #include <mkl.h>
    #include <mkl_service.h>
#endif

namespace settings {
    void configure_threads() {
        // Set the number of threads to be used
        tools::log->info("OpenMP | omp_threads {} | max active levels {} | dynamic {} | proc bind {} | num procs {}", omp_get_max_threads(),
                         omp_get_max_active_levels(), omp_get_dynamic(), omp_get_proc_bind(), omp_get_num_procs());

        /*
         * num_threads = omp_threads + num_std_threads = Eigen::nbThreads() + std_threads
         * If num_threads <  omp_threads --> num_threads  = omp_max_threads
         * If num_threads == omp_threads --> std_threads = 1 (can't be zero)
         * If num_threads >  omp_threads --> std_threads = num_threads - omp_max_threads
         * If num_threads == 0 or -1ul   --> num_threads = max_threads, and std_threads = num_threads - omp_threads (i.e. "the rest")
         */
        using namespace settings::threading;
        unsigned int omp_threads = static_cast<unsigned int>(omp_get_max_threads());
        unsigned int std_threads = omp_threads <= 1 or omp_threads == num_threads ? 1 : 0;
        num_threads              = std::clamp(num_threads, omp_threads, max_threads);
        std_threads              = std::clamp(num_threads - omp_threads + std_threads, num_threads - omp_threads, num_threads);

        std::string eigen_msg;
#if defined(EIGEN_USE_MKL_ALL)
        eigen_msg.append(" | EIGEN_USE_MKL_ALL");
#endif
#if defined(EIGEN_USE_BLAS)
        eigen_msg.append(" | EIGEN_USE_BLAS");
#endif
#if defined(EIGEN_USE_THREADS)
        eigen_msg.append(" | EIGEN_USE_THREADS");
#endif
#if defined(EIGEN_USE_THREADS)
        {
            tenx::threads::setNumThreads(std_threads);
            tools::log->info("Eigen3 | omp_threads {} | std_threads {} | max_threads {}{}", Eigen::nbThreads(), tenx::threads::num_threads, max_threads,
                             eigen_msg);
        }

#else
        if(settings::threading::num_threads > 1)
            tools::log->warn("EIGEN_USE_THREADS is not defined: "
                             "Failed to enable threading in Eigen::Tensor with stl_threads = {}",
                             settings::threading::num_threads);
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
    }
}