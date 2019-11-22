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
    size_t num_threads;
    #ifdef _OPENMP
    Eigen::ThreadPool       tp;
    Eigen::ThreadPoolDevice dev;
    OMP(size_t num_threads_):
    num_threads(num_threads_),
    tp(num_threads),
    dev(&tp, num_threads)
    {}
    #else
    Eigen::DefaultDevice dev;
    OMP(size_t num_threads_):
    num_threads(num_threads_)
    {}
    #endif

};
//
//namespace omp{
//    #ifdef _OPENMP
//    inline size_t num_threads  = Eigen::nbThreads();
//    inline Eigen::ThreadPool       tp (num_threads);
//    inline Eigen::ThreadPoolDevice dev(&tp,num_threads);
//    inline void change_num_threads(size_t threads){
//        num_threads = threads;
//        tp  = Eigen::ThreadPool       (num_threads);
//        dev = Eigen::ThreadPoolDevice (&tp,num_threads);
//
//    }
//    #else
//    inline Eigen::DefaultDevice dev;
//    #endif
//}

