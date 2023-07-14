#include "threads.h"
#include <thread>
#include <unsupported/Eigen/CXX11/Tensor>
//#if defined(_OPENMP)
//    #include <omp.h>
//#endif

namespace tenx::threads {
#if defined(EIGEN_USE_THREADS)
    std::unique_ptr<Eigen::ThreadPool>       tp;
    std::unique_ptr<Eigen::ThreadPoolDevice> dev;
    void                                     setNumThreads(unsigned int num) {
    #if defined(_OPENMP)
        auto pause = omp_pause_resource_all(omp_pause_resource_t::omp_pause_soft);
        if (pause != 0) throw std::runtime_error("omp_pause_resource_all(omp_pause_soft) failed");
    #endif

        num_threads = std::clamp<unsigned int>(num, 1, std::thread::hardware_concurrency());
        if(not dev or not tp or tp->NumThreads() != static_cast<int>(num_threads)) {
            tp  = std::make_unique<Eigen::ThreadPool>(num_threads);
            dev = std::make_unique<Eigen::ThreadPoolDevice>(tp.get(), num_threads);
        }
    }

    Eigen::ThreadPoolDevice &getDevice() {
        setNumThreads(num_threads);
        return *dev;
    }

#else
    void                                  setNumThreads([[maybe_unused]] int num) {}
    std::unique_ptr<Eigen::DefaultDevice> dev;
    Eigen::DefaultDevice                 &getDevice() {
        if(not dev) dev = std::make_unique<Eigen::DefaultDevice>();
        return *dev;
    }

#endif

}
