#include <general/iter.h>
#include <h5pp/h5pp.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_mpo.h>
#include <math/omp.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt_mps.h>

std::string solve(size_t ncv = 32, size_t maxIter = 30000, std::string_view method = "PRIMME_DEFAULT_METHOD") {
    h5pp::File h5file("../../output/krylov-17.h5", h5pp::FilePermission::READONLY);
    auto       mps = h5file.readDataset<Eigen::Tensor<double, 3>>("initial_mps");
    auto       mpo = h5file.readDataset<Eigen::Tensor<double, 4>>("mpo2");
    auto       enL = h5file.readDataset<Eigen::Tensor<double, 3>>("env2L");
    auto       enR = h5file.readDataset<Eigen::Tensor<double, 3>>("env2R");

    auto dims_mpo = mpo.dimensions();
    auto dims_mps = mps.dimensions();
    auto size     = mps.size();

    eig::solver solver_primme_sa;
    solver_primme_sa.config.tol           = 1e-12;
    solver_primme_sa.config.maxIter       = maxIter;
    solver_primme_sa.config.maxTime       = 1 * 60;
    solver_primme_sa.config.maxNev        = 1;
    solver_primme_sa.config.maxNcv        = ncv; // 128
    solver_primme_sa.config.compress      = true;
    solver_primme_sa.config.lib           = eig::Lib::PRIMME;
    solver_primme_sa.config.primme_method = method;
    solver_primme_sa.setLogLevel(2);

    // Reset the matrix
    tools::log->info("Finding excited state using shifted operator [(H-E)²-λ], with λ = 1.0 | primme SA | {} | mps {} | mpo {} | residual on ...", method,
                     dims_mps, dims_mpo);
    auto hamiltonian_squared = MatVecMPO<double>(enL.data(), enR.data(), mpo.data(), dims_mps, dims_mpo);

    solver_primme_sa.eigs(hamiltonian_squared, -1, -1, eig::Ritz::SA, eig::Form::SYMM, eig::Side::R, 1.0, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF,
                          mps.data());

    auto msg = fmt::format("Result: {:<32} | Found {:<5} | iters {:<7} |  matvecs {:<7} | time {:<13} s | min {:.16f}", method,
                           solver_primme_sa.result.meta.eigvecsR_found, solver_primme_sa.result.meta.iter, solver_primme_sa.result.meta.matvecs,
                           solver_primme_sa.result.meta.time_total, solver_primme_sa.result.get_eigvals()[0]);

    tools::log->info("{}\n", msg);

    return msg;
}

int main() {
    tools::Logger::setLogger(tools::log, "eig", 0, true);

    //    solve(32, 30000, "PRIMME_DEFAULT_METHOD");
    //    solve(32, 30000, "PRIMME_DYNAMIC");
    std::vector<std::string> msg;
    //    msg.push_back(solve(16, 30000, "PRIMME_DEFAULT_MIN_TIME"));
    //    msg.push_back(solve(16, 30000, "PRIMME_DEFAULT_MIN_MATVECS"));
    //    msg.push_back(solve(16, 30000, "PRIMME_Arnoldi"));
    //    msg.push_back(solve(16, 30000, "PRIMME_GD"));
    msg.push_back(solve(16, 30000, "PRIMME_GD_plusK"));
    msg.push_back(solve(16, 30000, "PRIMME_GD_Olsen_plusK"));
    msg.push_back(solve(16, 30000, "PRIMME_JD_Olsen_plusK"));
    //    msg.push_back(solve(16, 30000, "PRIMME_RQI"));
    //    msg.push_back(solve(16, 30000, "PRIMME_JDQR"));
    msg.push_back(solve(16, 30000, "PRIMME_JDQMR"));
    msg.push_back(solve(16, 30000, "PRIMME_JDQMR_ETol"));
    //    msg.push_back(solve(16, 30000, "PRIMME_STEEPEST_DESCENT"));
    msg.push_back(solve(16, 30000, "PRIMME_LOBPCG_OrthoBasis"));
    msg.push_back(solve(16, 30000, "PRIMME_LOBPCG_OrthoBasis_Window"));

    tools::log->info("Result summary:");
    for(auto &m : msg) tools::log->info("{}", m);
}