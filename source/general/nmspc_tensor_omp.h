//
// Created by david on 2019-10-06.
//

#pragma once

#if defined(_OPENMP) && defined(EIGEN_USE_THREADS)
    #include <omp.h>
    #include <thread>
#endif

#include <unsupported/Eigen/CXX11/Tensor>
#include <memory>
namespace Textra::omp {

#if defined(_OPENMP) && defined(EIGEN_USE_THREADS)
    inline std::unique_ptr<Eigen::ThreadPool>       tp;
    inline std::unique_ptr<Eigen::ThreadPoolDevice> dev;
    inline void                                     setNumThreads(int num_threads) {
        if(not tp) tp = std::make_unique<Eigen::ThreadPool>(num_threads);
        if(not dev and tp) dev = std::make_unique<Eigen::ThreadPoolDevice>(tp.get(), num_threads);
    }
#else
    inline std::unique_ptr<Eigen::DefaultDevice> dev = std::make_unique<Eigen::DefaultDevice>() ;
    inline void setNumThreads([[maybe_unused]] int num_threads) {
        if(not dev) dev = std::make_unique<Eigen::DefaultDevice>();
    }
#endif

}
