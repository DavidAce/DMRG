#pragma once

//#include <unsupported/Eigen/CXX11/Tensor>

#include <string>
// Handle NVCC/CUDA/SYCL
#if defined(__CUDACC__) || defined(__SYCL_DEVICE_ONLY__)
// Do not try asserts on CUDA and SYCL!
    #ifndef EIGEN_NO_DEBUG
        #define EIGEN_NO_DEBUG
    #endif

    #ifdef EIGEN_INTERNAL_DEBUGGING
        #undef EIGEN_INTERNAL_DEBUGGING
    #endif

    #ifdef EIGEN_EXCEPTIONS
        #undef EIGEN_EXCEPTIONS
    #endif

    // All functions callable from CUDA code must be qualified with __device__
    #ifdef __CUDACC__
        // Do not try to vectorize on CUDA and SYCL!
        #ifndef EIGEN_DONT_VECTORIZE
            #define EIGEN_DONT_VECTORIZE
        #endif

        #define EIGEN_DEVICE_FUNC __host__ __device__
        // We need cuda_runtime.h to ensure that that EIGEN_USING_STD_MATH macro
        // works properly on the device side
        #include <cuda_runtime.h>
    #else
        #define EIGEN_DEVICE_FUNC
    #endif

#else
    #define EIGEN_DEVICE_FUNC

#endif

#include <Eigen/src/Core/util/Macros.h>

#include <Eigen/src/Core/util/Constants.h>

#include <Eigen/src/Core/util/Meta.h>

#include <Eigen/src/Core/util/ForwardDeclarations.h>

#include <unsupported/Eigen/CXX11/src/Tensor/TensorMacros.h>

#include <unsupported/Eigen/CXX11/src/Tensor/TensorForwardDeclarations.h>