#include "nmspc_tensor_omp.h"
namespace Textra::omp {
    int num_threads = 1;
#if defined(EIGEN_USE_THREADS)
    std::unique_ptr<Eigen::ThreadPool>       tp;
    std::unique_ptr<Eigen::ThreadPoolDevice> dev;
    void                                     setNumThreads(int num) {
        num_threads = num;
        if(not dev or not tp or tp->NumThreads() != num_threads) {
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
    Eigen::DefaultDevice &                getDevice() {
        if(not dev) dev = std::make_unique<Eigen::DefaultDevice>();
        return *dev;
    }

#endif

}
