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

class OMP{
public:
    unsigned int num_threads;
    #if defined(_OPENMP) && defined(EIGEN_USE_THREADS)
    Eigen::ThreadPool       tp;
    Eigen::ThreadPoolDevice dev;
    explicit OMP(unsigned int num_threads_):
        num_threads(num_threads_),
        tp((int)num_threads),
        dev(&tp, (int)num_threads)
        {}
    explicit OMP():
        num_threads(std::thread::hardware_concurrency()),
        tp((int) std::thread::hardware_concurrency()),
        dev(&tp, (int) num_threads)
        {}
    #else
    Eigen::DefaultDevice dev;
    OMP(unsigned int num_threads_):
        num_threads(num_threads_)
        {}
    explicit OMP():
        num_threads(1)
        {}
    #endif

};
