#include <general/iter.h>
#include <h5pp/h5pp.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_mpo.h>
#include <math/tenx.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt_mps.h>

std::string solve(double tol, size_t nev, size_t ncv, std::optional<int> maxInner, size_t maxIter, const h5pp::File &h5file, std::string_view group,
                  eig::PrimmeMethod method) {
    auto                                     mpo   = h5file.readDataset<Eigen::Tensor<eig::real, 4>>(fmt::format("{}/mpo", group));
    auto                                     envL  = h5file.readDataset<Eigen::Tensor<eig::real, 3>>(fmt::format("{}/envL", group));
    auto                                     envR  = h5file.readDataset<Eigen::Tensor<eig::real, 3>>(fmt::format("{}/envR", group));
    auto                                     mpo2  = h5file.readDataset<Eigen::Tensor<eig::real, 4>>(fmt::format("{}/mpo2", group));
    auto                                     envL2 = h5file.readDataset<Eigen::Tensor<eig::real, 3>>(fmt::format("{}/envL2", group));
    auto                                     envR2 = h5file.readDataset<Eigen::Tensor<eig::real, 3>>(fmt::format("{}/envR2", group));
    std::vector<Eigen::Tensor<eig::real, 3>> mps_init;
    for(const auto &mps_name : h5file.findDatasets("mps_init_", group))
        mps_init.emplace_back(h5file.readDataset<Eigen::Tensor<eig::real, 3>>(fmt::format("{}/{}", group, mps_name)));

    auto ham  = MatVecMPO<double>(envL, envR, mpo);
    auto ham2 = MatVecMPO<double>(envL2, envR2, mpo2);
    auto dims = ham2.get_shape_mps();

    eig::solver s;
    s.config.ritz            = eig::Ritz::primme_smallest;
    s.config.form            = eig::Form::SYMM;
    s.config.shift_invert    = eig::Shinv::OFF;
    s.config.compute_eigvecs = eig::Vecs::ON;
    s.config.remove_phase    = eig::Dephase::OFF;

    s.config.tol                         = tol;
    s.config.maxIter                     = maxIter;
    s.config.maxTime                     = 30 * 60;
    s.config.maxNev                      = nev;
    s.config.maxNcv                      = ncv; // 128
    s.config.sigma                       = 1.0;
    s.config.compress                    = true;
    s.config.lib                         = eig::Lib::PRIMME;
    s.config.primme_max_inner_iterations = maxInner;
    s.config.primme_method               = method;
    s.config.logTime                     = 1;
    s.config.primme_extra                = &ham;
    for(auto &&[idx, mps] : iter::enumerate(mps_init)) s.config.initial_guess.push_back({mps.data(), static_cast<long>(idx)});

    s.setLogLevel(1);
    std::vector<eig::settings> configs;

    s.eigs(ham2);
    auto maxInnerStr = fmt::format("{:>4}", maxInner ? std::to_string(maxInner.value()) : "auto");
    auto msg =
        fmt::format(FMT_STRING("Result: method {:<32} | {:^10} | nev {:>3} | ncv {:>3} | maxInner {:>4} | tol {:8.2e} | iter {:>5} | op {:>5} | mv {:>7} | Found {:<5} | "
                               "∇fₘₐₓ {:8.2e} | res {:8.2e} | var {:>12.6e} | time {:>10.3f} s | t/op {:>10.3f} ms | t/mv {:>10.3f} ms | dims {}"),
                    eig::MethodToString(method), group, nev, ncv, maxInnerStr, tol, s.result.meta.iter,s.result.meta.num_op, s.result.meta.num_mv, s.result.meta.eigvecsR_found,
                    s.result.meta.last_grad_max, s.result.meta.last_res_norm, 1.0 + s.result.get_eigvals()[0], s.result.meta.time_total,
                    s.result.meta.time_total * 1000.0 / s.result.meta.num_op,s.result.meta.time_total * 1000.0 / s.result.meta.num_mv, dims);

    tools::log->info("{}", msg);

    return msg;
}

namespace threading {
    int stl_threads = 2;
    int omp_threads = 1;
}

int main() {

    tools::Logger::setLogger(tools::log, "bench", 0, true);

// Take care of threading
// Set the number of threads to be used
#if defined(EIGEN_USE_THREADS)
    if(threading::stl_threads <= 0) { threading::stl_threads = (int) std::thread::hardware_concurrency(); }
    tenx::omp::setNumThreads(threading::stl_threads);
    tools::log->info("Using Eigen Tensor with {} C++11 threads", tenx::omp::num_threads);
#else
    if(threading::stl_threads > 1)
        tools::log->warn("EIGEN_USE_THREADS is not defined: "
                         "Failed to enable threading in Eigen::Tensor with stl_threads = {}",
                         threading::stl_threads);
#endif

#if defined(_OPENMP)
    if(threading::omp_threads <= 0) { threading::omp_threads = (int) std::thread::hardware_concurrency(); }
    omp_set_num_threads(threading::omp_threads); // Should only need this. Both Eigen (non-Tensor) and MKL listen to this
    //    omp_set_max_active_levels(1);
    tools::log->info("Using OpenMP with {} threads and {} active levels", omp_get_max_threads(), omp_get_max_active_levels());
#endif
#if defined(OPENBLAS_AVAILABLE)
    auto        openblas_parallel_mode = openblas_get_parallel();
    std::string openblas_parallel_str;
    if(openblas_parallel_mode == 0) openblas_parallel_str = "seq";
    if(openblas_parallel_mode == 1) openblas_parallel_str = "threads";
    if(openblas_parallel_mode == 2) openblas_parallel_str = "openmp";
    if(openblas_parallel_mode == 1) openblas_set_num_threads(threading::omp_threads); // Use the omp_threads level for blas-related threading
    tools::log->info("{} compiled with parallel mode [{}] for target {} with config {} | multithread threshold {} | running with {} threads", OPENBLAS_VERSION,
                     openblas_parallel_str, openblas_get_corename(), openblas_get_config(), OPENBLAS_GEMM_MULTITHREAD_THRESHOLD, openblas_get_num_threads());
#endif
#if defined(MKL_AVAILABLE)
    tools::log->info("Using Intel MKL with {} threads", mkl_get_max_threads());
#endif

    auto filename = fmt::format("{}/primme_mps.h5", TEST_MATRIX_DIR);
    if(h5pp::fs::exists(filename)) {
        auto h5file = h5pp::File(filename, h5pp::FilePermission::READONLY);
        for(const auto &group : h5file.findGroups("mps-")) {
            if(group.find("mps-108") == std::string::npos) continue;
            auto dims = h5file.readAttribute<std::vector<long>>("dimensions", group);
            if(dims[0] * dims[1] * dims[2] < 10000) {
                tools::log->info("Skipped problem size {} = {} < 10000", dims, dims[0] * dims[1] * dims[2]);
                continue;
            }
            if(dims[0] * dims[1] * dims[2] > 64000) {
                tools::log->info("Skipped problem size {} = {} > 64000", dims, dims[0] * dims[1] * dims[2]);
                continue;
            }
            std::vector<std::string> msg;
            //        msg.push_back(solve(1e-14, 1, 16, std::nullopt, 30000, h5file, group, eig::PrimmeMethod::PRIMME_GD_plusK));
            //        msg.push_back(solve(1e-14, 1, 16, -1, 30000, h5file, group, eig::PrimmeMethod::PRIMME_DYNAMIC));
            //        msg.push_back(solve(1e-14, 1, 16, -1, 30000, h5file, group, eig::PrimmeMethod::PRIMME_GD_plusK));
            //        msg.push_back(solve(1e-14, 1, 16, std::nullopt, 30000, h5file, group, eig::PrimmeMethod::PRIMME_GD_Olsen_plusK));
            msg.push_back(solve(1e-14, 1, 16, -1, 30000, h5file, group, eig::PrimmeMethod::PRIMME_GD_Olsen_plusK));
            //        msg.push_back(solve(1e-14, 1, 16, std::nullopt, 30000, h5file, group, eig::PrimmeMethod::PRIMME_JD_Olsen_plusK));
            //        msg.push_back(solve(1e-14, 1, 16, -1, 30000, h5file, group, eig::PrimmeMethod::PRIMME_JD_Olsen_plusK));
            //        msg.push_back(solve(1e-14, 1, 16, std::nullopt, 30000, h5file, group, eig::PrimmeMethod::PRIMME_JDQMR));
            msg.push_back(solve(1e-14, 1, 16, std::nullopt, 30000, h5file, group, eig::PrimmeMethod::PRIMME_JDQMR_ETol));
            msg.push_back(solve(1e-14, 1, 16, -1, 30000, h5file, group, eig::PrimmeMethod::PRIMME_Arnoldi));
            //        msg.push_back(solve(1e-14, 1, 16, std::nullopt, 30000, h5file, group, eig::PrimmeMethod::PRIMME_LOBPCG_OrthoBasis));
            //        msg.push_back(solve(1e-14, 1, 16, -1, 30000, h5file, group, eig::PrimmeMethod::PRIMME_LOBPCG_OrthoBasis));
            //        msg.push_back(solve(1e-14, 1, 16, std::nullopt, 30000, h5file, group, eig::PrimmeMethod::PRIMME_LOBPCG_OrthoBasis_Window));
            //        msg.push_back(solve(1e-14, 1, 16, -1, 30000, h5file, group, eig::PrimmeMethod::PRIMME_LOBPCG_OrthoBasis_Window));
            tools::log->debug("Result summary:");
            for(auto &m : msg) tools::log->debug("{}", m);
        }
    }
}