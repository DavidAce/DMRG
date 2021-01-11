//
// Created by david on 2019-10-06.
//

#pragma once
#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>

namespace Textra::omp {
    extern int num_threads;
#if defined(EIGEN_USE_THREADS)
    extern std::unique_ptr<Eigen::ThreadPool>       tp;
    extern std::unique_ptr<Eigen::ThreadPoolDevice> dev;
    void                                            setNumThreads(int num);
    Eigen::ThreadPoolDevice &                       getDevice();
#else
    extern std::unique_ptr<Eigen::DefaultDevice> dev;
    void                                         setNumThreads([[maybe_unused]] int num);
    Eigen::DefaultDevice &                       getDevice();
#endif

}
