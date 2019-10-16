//
// Created by david on 2019-10-06.
//

#pragma once


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_USE_THREADS
#endif

#include <unsupported/Eigen/CXX11/Tensor>


namespace omp{
    #ifdef _OPENMP
    inline Eigen::ThreadPool       tp (Eigen::nbThreads());
    inline Eigen::ThreadPoolDevice dev(&tp,Eigen::nbThreads());
    #else
    inline Eigen::DefaultDevice dev;
    #endif
}

