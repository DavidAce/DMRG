#include <general/iter.h>
#include <h5pp/h5pp.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_mps.h>
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

    auto ham  = MatVecMps<double>(envL, envR, mpo);
    auto ham2 = MatVecMps<double>(envL2, envR2, mpo2);
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
    s.config.logTime                     = 10;
    s.config.primme_extra                = &ham;
    for(auto &&[idx, mps] : iter::enumerate(mps_init)) s.config.initial_guess.push_back({mps.data(), static_cast<long>(idx)});

    s.setLogLevel(2);
    std::vector<eig::settings> configs;

    s.eigs(ham2);
    auto maxInnerStr = fmt::format("{:>4}", maxInner ? std::to_string(maxInner.value()) : "auto");
    auto msg = fmt::format(FMT_STRING("Result: method {:<32} | {:^10} | nev {:>3} | ncv {:>3} | maxInner {:>4} | tol {:8.2e} | iter {:>5} | mv {:>7} | Found {:<5} | "
                                      "∇fₘₐₓ {:8.2e} | res {:8.2e} | var {:>12.6e} | time {:>10.3f} s | dims {}"),
                           eig::MethodToString(method), group, nev, ncv, maxInnerStr, tol, s.result.meta.iter, s.result.meta.matvecs, s.result.meta.eigvecsR_found,
                           s.result.meta.last_grad_max, s.result.meta.last_res_norm, 1.0+s.result.get_eigvals()[0], s.result.meta.time_total, dims);

    tools::log->info("{}", msg);

    return msg;
}

int main() {
    tools::Logger::setLogger(tools::log, "eig", 0, true);

    auto h5file = h5pp::File(fmt::format("{}/primme_mps.h5", TEST_MATRIX_DIR), h5pp::FilePermission::READONLY);

    for(const auto &group : h5file.findGroups("mps-")) {
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