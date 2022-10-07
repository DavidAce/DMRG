#include "threads.h"
#include <memory>
#include <omp.h>
#include <thread>
#include <unsupported/Eigen/CXX11/Tensor>
namespace tenx::threads {
#if defined(EIGEN_USE_THREADS)
    std::unique_ptr<Eigen::ThreadPool>       tp;
    std::unique_ptr<Eigen::ThreadPoolDevice> dev;
    void                                     setNumThreads(unsigned int num) {
        num_threads = std::clamp<unsigned int>(num, static_cast<unsigned int>(1), std::thread::hardware_concurrency());
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
