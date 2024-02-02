#include "threads.h"
#include "math/cast.h"
#include <thread>
#include <unsupported/Eigen/CXX11/Tensor>
#if defined(_OPENMP)
    #include <omp.h>
#endif

namespace tenx::threads {
    template<typename T>
    requires std::is_integral_v<T>
    void setNumThreads([[maybe_unused]] T num) noexcept {
#if defined(EIGEN_USE_THREADS)
        internal::num_threads = static_cast<unsigned int>(num);
#endif
    }
    template void setNumThreads(int num) noexcept;
    template void setNumThreads(long num) noexcept;
    template void setNumThreads(unsigned int num) noexcept;
    template void setNumThreads(unsigned long num) noexcept;

#if defined(EIGEN_USE_THREADS)
    std::unique_ptr<internal::ThreadPoolWrapper> internal::singleThreadWrapper;
    std::unique_ptr<internal::ThreadPoolWrapper> internal::multiThreadWrapper;

    int getNumThreads() noexcept {
    #if defined(_OPENMP)
        if(omp_in_parallel()) { return 1; } // Avoid simultaneous parallelization
    #endif
        return static_cast<int>(internal::num_threads);
    }

    internal::ThreadPoolWrapper::ThreadPoolWrapper(int nt)
        : tp(std::make_unique<Eigen::ThreadPool>(nt)), dev(std::make_unique<Eigen::ThreadPoolDevice>(tp.get(), nt)) {}

//    internal::ThreadPoolWrapper &get() noexcept {
//        if(not internal::singleThreadWrapper) internal::singleThreadWrapper = std::make_unique<internal::ThreadPoolWrapper>(1);
//    #if defined(_OPENMP)
//        if(omp_in_parallel()) { // Avoid simultaneous parallelization
//            return *internal::singleThreadWrapper;
//        }
//    #endif
//        if(not internal::multiThreadWrapper or //
//           (internal::multiThreadWrapper and static_cast<int>(internal::num_threads) != internal::multiThreadWrapper->dev->numThreads()))
//            internal::multiThreadWrapper = std::make_unique<internal::ThreadPoolWrapper>(internal::num_threads);
//
//        return *internal::multiThreadWrapper;
//    }




    const std::unique_ptr<internal::ThreadPoolWrapper> &get() noexcept {
        if(not internal::singleThreadWrapper) internal::singleThreadWrapper = std::make_unique<internal::ThreadPoolWrapper>(1);
    #if defined(_OPENMP)
        if(omp_in_parallel()) { // Avoid simultaneous parallelization
//            return std::make_unique<internal::ThreadPoolWrapper>(1);
            return internal::singleThreadWrapper;
        }
    #endif
        if(not internal::multiThreadWrapper or //
           (internal::multiThreadWrapper and static_cast<int>(internal::num_threads) != internal::multiThreadWrapper->dev->numThreads()))
            internal::multiThreadWrapper = std::make_unique<internal::ThreadPoolWrapper>(internal::num_threads);
//        return std::make_unique<internal::ThreadPoolWrapper>(static_cast<int>(internal::num_threads));
          return internal::multiThreadWrapper;

    }

#else
    internal::DefaultDeviceWrapper::DefaultDeviceWrapper(int nt) : dev(std::make_unique<Eigen::DefaultDevice>()) {}

    void                                  setNumThreads([[maybe_unused]] int num) {}
    std::unique_ptr<Eigen::DefaultDevice> dev;
    Eigen::DefaultDevice                 &getDevice() {
        if(not dev) dev = std::make_unique<Eigen::DefaultDevice>();
        return *dev;
    }

#endif

}
