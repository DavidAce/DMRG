//
// Created by david on 2019-10-06.
//

#pragma once


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_USE_THREADS
#endif

#include <unsupported/Eigen/CXX11/Tensor>

class OMP{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    size_t num_threads;
    #ifdef _OPENMP
    Eigen::ThreadPool       tp;
    Eigen::ThreadPoolDevice dev;
    explicit OMP(size_t num_threads_):
        num_threads(num_threads_),
        tp(num_threads),
        dev(&tp, num_threads)
        {}
    explicit OMP():
        num_threads(std::thread::hardware_concurrency()),
        tp(std::thread::hardware_concurrency()),
        dev(&tp, num_threads)
        {}
    #else
    Eigen::DefaultDevice dev;
    OMP(size_t num_threads_):
        num_threads(num_threads_)
        {}
    explicit OMP():
        num_threads(1)
        {}
    #endif

};
