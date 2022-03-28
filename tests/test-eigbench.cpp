#include <general/iter.h>
#include <h5pp/h5pp.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_mpo.h>
#include <math/tenx.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt_mps.h>

std::string solve(const eig::settings &config, const h5pp::File &h5file, std::string_view group) {
    std::string tgt;
    {
        auto dims  = h5file.readAttribute<std::array<long, 3>>("dims", group);
        auto time  = h5file.readAttribute<double>("time", group);
        auto iter  = h5file.readAttribute<long>("iter", group);
        auto eval  = h5file.readAttribute<double>("eval", group);
        auto rnorm = h5file.readAttribute<double>("rnorm", group);
        auto mv    = h5file.readAttribute<long>("mv", group);
        auto op    = h5file.readAttribute<long>("op", group);
        tgt = fmt::format("Target | iter {:>5} | op {:>5} | mv {:>7} | rnorm {:8.2e} | var {:>12.6e} | time {:>10.3f} s | dims {}", iter, op, mv, rnorm, eval,
                          time, dims);
    }

    Eigen::Tensor<eig::real, 3> mps   = h5file.readDataset<Eigen::Tensor<eig::cplx, 3>>(fmt::format("{}/mps", group)).real();
    auto                        mpo2  = h5file.readDataset<Eigen::Tensor<eig::real, 4>>(fmt::format("{}/mpo2", group));
    auto                        envL2 = h5file.readDataset<Eigen::Tensor<eig::real, 3>>(fmt::format("{}/envL2", group));
    auto                        envR2 = h5file.readDataset<Eigen::Tensor<eig::real, 3>>(fmt::format("{}/envR2", group));

    auto ham2 = MatVecMPO<double>(envL2, envR2, mpo2);
    auto dims = ham2.get_shape_mps();

    eig::solver s;
    s.config = config;
    s.config.initial_guess.push_back({mps.data(), 0});

    s.setLogLevel(1);
    s.eigs(ham2);
    auto rnorm       = *std::max_element(s.result.meta.residual_norms.begin(), s.result.meta.residual_norms.end());
    auto maxInnerStr = fmt::format("{:>4}", s.config.primme_max_inner_iterations ? std::to_string(s.config.primme_max_inner_iterations.value()) : "auto");
    auto msg         = fmt::format(
                FMT_STRING("Result: method {:<32} | {:^10} | nev {:>3} | ncv {:>3} | maxInner {:>4} | tol {:8.2e} | iter {:>5} | op {:>5} | mv {:>7} | Found {:<5} | "
                                   "rnorm {:8.2e} | var {:>12.6e} | time {:>10.3f} s | t/op {:>10.3f} ms | t/mv {:>10.3f} ms | dims {}"),
                eig::MethodToString(config.primme_method.value()), group, s.config.maxNev.value(), s.config.maxNcv.value(), maxInnerStr, s.config.tol.value(),
                s.result.meta.iter, s.result.meta.num_op, s.result.meta.num_mv, s.result.meta.eigvecsR_found, rnorm, s.result.get_eigvals()[0],
                s.result.meta.time_total, s.result.meta.time_total * 1000.0 / static_cast<double>(s.result.meta.num_op),
                s.result.meta.time_total * 1000.0 / static_cast<double>(s.result.meta.num_mv), dims);

    tools::log->info("{}", tgt);
    tools::log->info("{}", msg);

    return msg;
}

namespace threading {
    int stl_threads = 4;
    int omp_threads = 4;
}

struct Params {};

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

    std::vector<eig::settings> configs(1);
    // https://www.cs.wm.edu/~andreas/software/doc/appendix.html#c.primme_params.eps
    configs[0].tol                         = 1e-12; // 1e-12 is good. This Sets "eps" in primme, see link above.
    configs[0].maxIter                     = 20000;
    configs[0].maxNev                      = 1;
    configs[0].maxNcv                      = 4;
    configs[0].compress                    = false;
    configs[0].maxTime                     = 2 * 60 * 60; // Two hours
    configs[0].lib                         = eig::Lib::PRIMME;
    configs[0].ritz                        = eig::Ritz::SA;
    configs[0].compute_eigvecs             = eig::Vecs::ON;
    configs[0].loglevel                    = 2;
    configs[0].primme_max_inner_iterations = 0;
    configs[0].primme_method               = eig::PrimmeMethod::PRIMME_GD_plusK; // eig::PrimmeMethod::PRIMME_JDQMR;
                                                                                 //    configs[0].primme_projection
    auto filename = fmt::format("{}/eigs.h5", TEST_MATRIX_DIR);
    if(h5pp::fs::exists(filename)) {
        auto                     h5file = h5pp::File(filename, h5pp::FilePermission::READONLY);
        std::vector<std::string> msg;

        for(const auto &group : h5file.findGroups("eigs-")) {
            if(group.find("eigs-0") != std::string::npos) continue;
            tools::log->info("Running group: {}", group);
            for(const auto &config : configs) { msg.push_back(solve(config, h5file, group)); }
        }
        tools::log->debug("Result summary:");
        for(auto &m : msg) tools::log->debug("{}", m);
    }
}