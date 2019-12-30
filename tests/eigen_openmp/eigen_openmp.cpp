

#ifdef _OPENMP
#include <omp.h>
#include <thread>
#endif

#ifdef OpenBLAS_AVAILABLE
#include <cblas.h>
#include <openblas_config.h>
#endif

#ifdef MKL_AVAILABLE
#include <mkl_service.h>
#include <mkl.h>
#endif

#include <iostream>
#include <complex>
#include <general/nmspc_omp.h>
#include <general/nmspc_tensor_extra.h>
#include <Eigen/Core>






int main(){

    #ifdef _OPENMP
        omp_set_num_threads(std::thread::hardware_concurrency());
        Eigen::setNbThreads(std::thread::hardware_concurrency());
        std::cout << "Using Eigen  with " << Eigen::nbThreads() << " threads" << std::endl;
        std::cout << "Using OpenMP with " << omp_get_max_threads() << " threads" << std::endl;

        #ifdef OpenBLAS_AVAILABLE
        openblas_set_num_threads(std::thread::hardware_concurrency());
                    std::cout << OPENBLAS_VERSION
                              << " compiled with parallel mode " << openblas_get_parallel()
                              << " for target " << openblas_get_corename()
                              << " with config " << openblas_get_config()
                              << " with multithread threshold " << OPENBLAS_GEMM_MULTITHREAD_THRESHOLD
                              << ". Running with " << openblas_get_num_threads() << " thread(s)" << std::endl;
        #endif

        #ifdef MKL_AVAILABLE
        mkl_set_num_threads(std::thread::hardware_concurrency());
        std::cout << "Using Intel MKL with " << mkl_get_max_threads() << " threads" << std::endl;
        #endif
    #endif
    OMP omp;
    Eigen::Tensor<double,2> tensorA(100,100);
    Eigen::Tensor<double,2> tensorB(100,100);
    tensorA.setRandom();
    tensorB.setRandom();
    Eigen::Tensor<double,2> tensorC(100,100);
    tensorC.device(omp.dev) = tensorA.contract(tensorB, Textra::idx({1},{0}));
    std::cout << "tensorC: \n" << tensorC(0,0) << std::endl;
}