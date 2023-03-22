#include "threading.h"
#include "math/tenx/threads.h"
#include "settings.h"
#include "tools/common/log.h"
#include <cstdlib>
#include <omp.h>
#include <optional>
#if defined(OPENBLAS_AVAILABLE)
    #include <openblas/cblas.h>
    #include <openblas/openblas_config.h>
#endif

#if defined(MKL_AVAILABLE)
    #include <mkl.h>
    #include <mkl_service.h>
#endif

namespace settings {

    inline std::optional<std::string> get_env(std::string_view key) {
        if(key.empty()) throw std::invalid_argument("Value requested for the empty-name environment variable");
        const char *ev_val = std::getenv(key.data());
        if(ev_val == nullptr) return std::nullopt;
        return std::string(ev_val);
    }

    void configure_threads() {
        // Set the number of threads to be used
        std::string omp_proc_bind;
        switch(omp_get_proc_bind()) {
            /* clang-format off */
            case omp_proc_bind_false   : {omp_proc_bind = "false";break;}
            case omp_proc_bind_true    : {omp_proc_bind = "true";break;}
            case omp_proc_bind_primary : {omp_proc_bind = "primary";break;}
            case omp_proc_bind_close   : {omp_proc_bind = "close";break;}
            case omp_proc_bind_spread  : {omp_proc_bind = "spread";break;}
            default: omp_proc_bind = "unknown";
                /* clang-format on */
        }
        tools::log->info("OpenMP | omp_threads {} | max active levels {} | dynamic {} | proc bind [{}] | num procs {}", omp_get_max_threads(),
                         omp_get_max_active_levels(), omp_get_dynamic(), omp_proc_bind, omp_get_num_procs());

        using namespace settings::threading;
        auto omp_threads = static_cast<unsigned int>(omp_get_max_threads());
        num_threads      = std::clamp(num_threads, omp_threads, max_threads);
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
            tenx::threads::setNumThreads(num_threads);
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
        auto envcoretype = get_env("OPENBLAS_CORETYPE");
        if(envcoretype) tools::log->info("Detected environment variable: OPENBLAS_CORETYPE={}", envcoretype.value());
        auto        openblas_parallel_mode = openblas_get_parallel();
        std::string openblas_parallel_str;
        if(openblas_parallel_mode == 0) openblas_parallel_str = "sequential";
        if(openblas_parallel_mode == 1) openblas_parallel_str = "threads";
        if(openblas_parallel_mode == 2) openblas_parallel_str = "openmp";
        if(openblas_parallel_mode == 1) openblas_set_num_threads(num_threads); // Use this for blas-related threading
        tools::log->info("{} NUM_THREADS={} | parallel_mode={} | corename={} | gemm_multithread_threshold={}", openblas_get_config(),
                         openblas_get_num_threads(), openblas_parallel_str, openblas_get_corename(), OPENBLAS_GEMM_MULTITHREAD_THRESHOLD);
#endif
#if defined(MKL_AVAILABLE)
        tools::log->info("Intel MKL | max_threads {}", mkl_get_max_threads());
#endif
        if(settings::threading::show_threads) exit(0);
    }
}