#include "threading.h"
#include "math/tenx/threads.h"
#include "settings.h"
#include "tools/common/log.h"
#include <cstdlib>
#if defined(_OPENMP)
    #include <omp.h>
#endif
#include <optional>

#if defined(MKL_AVAILABLE)
    #include <mkl_service.h>
#elif defined(OPENBLAS_AVAILABLE)
    #include <openblas/cblas.h>
    #include <openblas/openblas_config.h>
#elif __has_include(<cblas.h>)
    #include <cblas.h>
    #if defined(OPENBLAS_VERSION)
        #define OPENBLAS_AVAILABLE
    #endif
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
        using namespace settings::threading;
        unsigned int omp_threads = 1;
#if defined(_OPENMP)
        std::string omp_proc_bind;
        switch(omp_get_proc_bind()) {
            /* clang-format off */
            case 0 : {omp_proc_bind = "false";break;}
            case 1 : {omp_proc_bind = "true";break;}
            case 2 : {omp_proc_bind = "primary";break;}
            case 3 : {omp_proc_bind = "close";break;}
            case 4 : {omp_proc_bind = "spread";break;}
            default: omp_proc_bind = "unknown";
                /* clang-format on */
        }
        tools::log->info("OpenMP | omp_max_threads {} | omp_max_active_levels {} | omp_dynamic {} | omp_proc_bind [{}] | omp_num_procs {}",
                         omp_get_max_threads(), omp_get_max_active_levels(), omp_get_dynamic(), omp_proc_bind, omp_get_num_procs());

        omp_threads = static_cast<unsigned int>(omp_get_max_threads());
#endif

        num_threads = std::clamp(num_threads, omp_threads, max_threads);
        std::string eigen_msg;
#if defined(EIGEN_USE_MKL_ALL)
        eigen_msg.append(" | EIGEN_USE_MKL_ALL");
#endif
#if defined(EIGEN_USE_BLAS)
        eigen_msg.append(" | EIGEN_USE_BLAS");
#endif
#if defined(EIGEN_USE_THREADS)
        eigen_msg.append(" | EIGEN_USE_THREADS");
        tenx::threads::setNumThreads(num_threads);
#else
        if(settings::threading::num_threads > 1)
            tools::log->warn("EIGEN_USE_THREADS is not defined: "
                             "Failed to enable threading in Eigen::Tensor with stl_threads = {}",
                             settings::threading::num_threads);
#endif
        tools::log->info("Eigen3 | omp_threads {} | std_threads {} | max_threads {}{}", Eigen::nbThreads(), tenx::threads::num_threads, max_threads, eigen_msg);
#if defined(OPENBLAS_AVAILABLE)
        auto envcoretype = get_env("OPENBLAS_CORETYPE");
        if(envcoretype) tools::log->info("Detected environment variable: OPENBLAS_CORETYPE={}", envcoretype.value());
        auto        openblas_parallel_mode = openblas_get_parallel();
        std::string openblas_parallel_str;
        if(openblas_parallel_mode == 0) openblas_parallel_str = "sequential";
        if(openblas_parallel_mode == 1) openblas_parallel_str = "threads";
        if(openblas_parallel_mode == 2) openblas_parallel_str = "openmp";
        if(openblas_parallel_mode == 1) openblas_set_num_threads(num_threads); // Use this for blas-related threading
        tools::log->info("{} NUM_THREADS={} | parallel_mode={} | corename={}", openblas_get_config(), openblas_get_num_threads(), openblas_parallel_str,
                         openblas_get_corename());
#endif
#if defined(MKL_AVAILABLE)
        tools::log->info("Intel MKL | max_threads {}", mkl_get_max_threads());
#endif
        if(settings::threading::show_threads) exit(0);
    }
}