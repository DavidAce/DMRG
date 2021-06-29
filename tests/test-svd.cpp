
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif
#include <general/nmspc_tensor_omp.h>

#ifdef OpenBLAS_AVAILABLE
#include <cblas.h>
    #include <openblas_config.h>
#endif

#ifdef MKL_AVAILABLE
#include <mkl.h>
#include <mkl_service.h>
#endif
#include <thread>

#include <Eigen/Core>
#include <h5pp/h5pp.h>
#include <math/svd.h>




TEST_CASE("Singular value decomposition in Eigen and Lapacke", "[svd]") {
    SECTION("Test all svd decompositions") {
        using reciter = h5pp::fs::recursive_directory_iterator;
        svd::settings svd_settings;
        svd_settings.threshold = 1e-8;
        svd_settings.loglevel  = 0;
        svd_settings.use_lapacke = false;
        svd::solver svd(svd_settings);

        [[maybe_unused]] Eigen::MatrixXcd U1,V1,U2,V2;
        Eigen::VectorXcd S1,S2;
        for(auto &item : reciter(std::string(TEST_MATRIX_DIR))) {
            if(item.path().filename().string().find("svdmatrix") == std::string::npos) continue;
            if(item.path().extension() != ".h5") continue;
            size_t     logLevel = 2;
            h5pp::File file(item.path().string(), h5pp::FilePermission::READONLY, logLevel);
            auto       matrix = file.readDataset<Eigen::MatrixXcd>("svdmatrix");
            svd.use_lapacke   = false;
            std::tie(U1, S1, V1) = svd.decompose(matrix);
            svd.use_lapacke   = true;
            std::tie(U2, S2, V2) = svd.decompose(matrix);
            double differenceS = std::log10((S2.array() - S1.array()).cwiseAbs().sum());
            fmt::print("S {:<32} diff {:.24f}\n", item.path().filename().string(), differenceS);
            REQUIRE(differenceS < -12);
        }
    }
}


int main(int argc, char **argv){

// Set the number of threads to be used
    [[maybe_unused]] int num_threads = 1;
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
    Eigen::setNbThreads(num_threads);
    Textra::omp::setNumThreads(num_threads);
#ifdef OpenBLAS_AVAILABLE
    openblas_set_num_threads(num_threads);
    std::cout << OPENBLAS_VERSION
              << " compiled with parallel mode " << openblas_get_parallel()
              << " for target " << openblas_get_corename()
              << " with config " << openblas_get_config()
              << " with multithread threshold " << OPENBLAS_GEMM_MULTITHREAD_THRESHOLD
              << ". Running with " << openblas_get_num_threads() << " thread(s)" << std::endl;
#endif

#ifdef MKL_AVAILABLE
    mkl_set_num_threads(num_threads);
    std::cout <<  "Using Intel MKL with " <<  mkl_get_max_threads() << " threads" << std::endl;
#endif

#endif

    return Catch::Session().run(argc, argv);
}