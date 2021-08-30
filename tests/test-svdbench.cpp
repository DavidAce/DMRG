
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#ifdef _OPENMP
    #include <omp.h>
#endif
#include <math/tenx.h>
#include <tensors/site/mps/MpsSite.h>
#include <tid/tid.h>
#include <tools/common/split.h>

#ifdef OPENBLAS_AVAILABLE
    #include <cblas.h>
    #include <openblas_config.h>
#endif

#ifdef MKL_AVAILABLE
    #include <mkl.h>
    #include <mkl_service.h>
#endif
#include <Eigen/Core>
#include <h5pp/h5pp.h>
#include <math/svd.h>
#include <thread>

TEST_CASE("Singular value decomposition in Eigen and Lapacke", "[svd]") {
    SECTION("Bench split functions") {
        svd::settings svd_settings;
        svd_settings.threshold  = 1e-8;
        svd_settings.loglevel   = 2;
        svd_settings.use_bdc    = true;
        svd_settings.switchsize = 32;

        h5pp::File h5file(fmt::format("{}/svd-benchmark.h5", TEST_MATRIX_DIR), h5pp::FilePermission::READONLY, 2);
        double     t_split_sum = 0;
        for(const auto &multisite_tensor_name : h5file.findDatasets("multisite_tensor", "/", -1)) {
            auto dims = h5file.getDatasetDimensions(multisite_tensor_name);
            auto size = dims[0] * dims[1] * dims[2];
            if(size < 131072) continue;

            auto multisite_tensor = h5file.readDataset<Eigen::Tensor<std::complex<double>, 3>>(multisite_tensor_name);

            auto spin_dims       = h5file.readAttribute<std::vector<long>>("spin_dims", multisite_tensor_name);
            auto positions       = h5file.readAttribute<std::vector<size_t>>("positions", multisite_tensor_name);
            auto center_position = h5file.readAttribute<long>("center_position", multisite_tensor_name);
            auto chi_limit       = h5file.readAttribute<long>("chi_limit", multisite_tensor_name);
            auto t_split         = h5file.readAttribute<double>("t_split", multisite_tensor_name);
            t_split_sum += t_split;
            {
                svd_settings.use_lapacke = false;
                auto t_eig               = tid::tic_scope("eig");
                auto mps_list_eig        = tools::common::split::split_mps(multisite_tensor, spin_dims, positions, center_position, chi_limit, svd_settings);
            }
            {
                svd_settings.use_lapacke = true;
                auto t_eig               = tid::tic_scope("lpk");
                auto mps_list_lpk        = tools::common::split::split_mps(multisite_tensor, spin_dims, positions, center_position, chi_limit, svd_settings);
            }
            tools::log->info("eig {:8.2e} | lpk {:8.2e} | original {:8.2e} | {:<32} | dims {}", tid::get_unscoped("eig.split").get_last_interval(),
                             tid::get_unscoped("lpk.split").get_last_interval(), t_split, multisite_tensor_name, dims);
        }
        tools::log->info("eig {:8.2e} | lpk {:8.2e} | original {:8.2e} | TOTAL", tid::get_unscoped("eig.split").get_time(),
                         tid::get_unscoped("lpk.split").get_time(), t_split_sum);
        for(const auto &t : tid::get_tree()) tools::log->info("{}", t.str());
    }
}
namespace threading {
    int stl_threads = 4;
    int omp_threads = 4;
}
int main(int argc, char **argv) {
    tools::Logger::setLogger(tools::log, "bench", 0, true);

// Take care of threading
// Set the number of threads to be used
#if defined(EIGEN_USE_THREADS)
    if(threading::stl_threads <= 0) { threading::stl_threads = (int) std::thread::hardware_concurrency(); }
    tenx::omp::setNumThreads(threading::stl_threads);
    tools::log->info("Eigen3 Tensor | stl threads {}", tenx::omp::num_threads);
#else
    if(threading::stl_threads > 1)
        tools::log->warn("EIGEN_USE_THREADS is not defined: "
                         "Failed to enable threading in Eigen::Tensor with stl_threads = {}",
                         threading::stl_threads);
#endif

#if defined(_OPENMP)
    if(threading::omp_threads <= 0) { threading::omp_threads = (int) std::thread::hardware_concurrency(); }
    omp_set_num_threads(threading::omp_threads); // Should only need this. Both Eigen (non-Tensor) and MKL listen to this
    tools::log->info("OpenMP | threads {} | max active levels {}", omp_get_max_threads(), omp_get_max_active_levels());
#endif

#if defined(OPENBLAS_AVAILABLE)
    auto        openblas_parallel_mode = openblas_get_parallel();
    std::string openblas_parallel_str;
    if(openblas_parallel_mode == 0) openblas_parallel_str = "seq";
    if(openblas_parallel_mode == 1) openblas_parallel_str = "threads";
    if(openblas_parallel_mode == 2) openblas_parallel_str = "openmp";
    if(openblas_parallel_mode == 1) openblas_set_num_threads(threading::omp_threads); // Use the omp_threads level for blas-related threading
    tools::log->info("{} threads {} | parallel mode [{}:{}] | core type {} | config {} | multithread threshold {}", OPENBLAS_VERSION,
                     openblas_get_num_threads(), openblas_parallel_mode, openblas_parallel_str, openblas_get_corename(), openblas_get_config(),
                     OPENBLAS_GEMM_MULTITHREAD_THRESHOLD);
#endif
#if defined(MKL_AVAILABLE)
    tools::log->info("Intel MKL | threads {}", mkl_get_max_threads());
#endif

#if defined(EIGEN_USE_MKL_ALL)
    tools::log->info("Eigen3 | threads {} | EIGEN_USE_MKL_ALL", Eigen::nbThreads());
#elif defined(EIGEN_USE_BLAS)
    tools::log->info("Eigen3 | threads {} | EIGEN_USE_BLAS", Eigen::nbThreads());
#else
    tools::log->info("Eigen3 | threads {} | No BLAS backend", Eigen::nbThreads());
#endif

    return Catch::Session().run(argc, argv);
}