//
// Created by david on 2019-10-06.
//

#pragma once

#ifdef _OPENMP
    #include <omp.h>
#endif

// Make sure EIGEN_USE_THREADS is defined globally already
#include <thread>
#include <unsupported/Eigen/CXX11/Tensor>

namespace Textra::omp {

#if defined(_OPENMP) && defined(EIGEN_USE_THREADS)
    inline std::unique_ptr<Eigen::ThreadPool>       tp;
    inline std::unique_ptr<Eigen::ThreadPoolDevice> dev;
    inline void                                     setNumThreads(int num_threads) {
        if(not tp) tp = std::make_unique<Eigen::ThreadPool>(num_threads);
        if(not dev and tp) dev = std::make_unique<Eigen::ThreadPoolDevice>(tp.get(), num_threads);
    }
#else
    Eigen::DefaultDevice dev;
    inline void          setNumThreads([[maybe_unused]] int num_threads) {}
#endif
}
//
//class OMP {
//    public:
//    int num_threads = 1;
//#if defined(_OPENMP) && defined(EIGEN_USE_THREADS)
//    Eigen::ThreadPool       tp;
//    Eigen::ThreadPoolDevice dev;
//    explicit OMP(int num_threads_) : num_threads(num_threads_ > 1 ? num_threads_ : Eigen::nbThreads()), tp(num_threads), dev(&tp, num_threads) {}
//    explicit OMP() : num_threads(Eigen::nbThreads()), tp(num_threads), dev(&tp, num_threads) {}
//#else
//    Eigen::DefaultDevice dev;
//    OMP(int num_threads_) : num_threads(num_threads_) {}
//    explicit OMP() : num_threads(1) {}
//#endif
//};
