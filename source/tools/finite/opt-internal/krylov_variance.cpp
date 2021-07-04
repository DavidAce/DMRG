
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_iter.h>
#include <general/nmspc_tensor_extra.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_mpo.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_mps.h>

tools::finite::opt::opt_mps tools::finite::opt::internal::krylov_variance_optimization(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                                                       const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                                       OptSpace optSpace) {
    using namespace internal;
    using namespace settings::precision;
    if(not tensors.model->is_reduced()) throw std::runtime_error("krylov_variance_optimization requires energy-reduced MPO²");
    if(tensors.model->is_compressed_mpo_squared()) throw std::runtime_error("krylov_variance_optimization requires non-compressed MPO²");

    auto t_eig    = tools::common::profile::get_default_prof()["t_eig"]->tic_token();
    auto dims_mps = initial_mps.get_tensor().dimensions();
    auto size     = initial_mps.get_tensor().size();

    // NOTE:
    // We can't use the MPO-shift trick here!
    // The lower-left corner of the MPO is NOT where the energy is subtracted on a double-layer MPO,
    // especially when using MPO compression.
    tools::log->debug("Finding excited state: smallest eigenstate of (H-E)² | dims {} = {}", dims_mps, size);

    eig::solver solver_full, solver_shift, solver_max, solver_primme_sa;
    solver_max.config.tol             = settings::precision::eig_tolerance;
    solver_max.config.maxIter         = 200;
    solver_max.config.maxNev          = 1;
    solver_max.config.maxNcv          = 16;
    solver_max.config.compute_eigvecs = eig::Vecs::OFF;
    solver_max.config.compress        = true;
    solver_max.setLogLevel(2);

    solver_shift.config.tol      = settings::precision::eig_tolerance;
    solver_shift.config.maxIter  = 400;
    solver_shift.config.maxTime  = 30 * 60;
    solver_shift.config.maxNev   = 1;
    solver_shift.config.maxNcv   = 64;
    solver_shift.config.compress = true;
    solver_shift.setLogLevel(2);
    //

    solver_primme_sa.config.tol      = settings::precision::eig_tolerance;
    solver_primme_sa.config.maxIter  = 800000;
    solver_primme_sa.config.maxTime  = 30 * 60;
    solver_primme_sa.config.maxNev   = 1;
    solver_primme_sa.config.maxNcv   = 64; // 128
    solver_primme_sa.config.compress = true;
    solver_primme_sa.config.lib      = eig::Lib::PRIMME;
    solver_primme_sa.setLogLevel(2);

    if(optType == OptType::REAL) {
        tools::log->trace("- Generating real-valued multisite components");
        const auto              &env2      = tensors.get_multisite_var_blk();
        Eigen::Tensor<double, 3> env2L     = env2.L.real();
        Eigen::Tensor<double, 3> env2R     = env2.R.real();
        Eigen::Tensor<double, 4> mpo2      = tensors.get_multisite_mpo_squared().real();
        Eigen::Tensor<double, 3> residual  = initial_mps.get_tensor().real();
        auto                     dims_mpo2 = mpo2.dimensions();

        if(size <= 1024) {
            //        if(size <= settings::precision::max_size_full_diag) {
            tools::log->trace("Full diagonalization of (H-E)²", dims_mps, size);
            auto matrix = tools::finite::opt::internal::get_multisite_hamiltonian_squared_matrix<double>(*tensors.model, *tensors.edges);
            solver_full.eig(matrix.data(), matrix.rows());
        } else {
            tools::log->trace("Defining reduced Hamiltonian-squared matrix-vector product");
            MatVecMPO<double> hamiltonian_squared(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);

            if(not tensors.model->is_compressed_mpo_squared()) {
                tools::log->trace("Finding largest-magnitude state");
                solver_max.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::OFF,
                                eig::Dephase::OFF);
                auto eigvals_ham_sq = eig::view::get_eigvals<double>(solver_max.result);
                auto largest_eigval = eigvals_ham_sq.cwiseAbs().maxCoeff() + 1.0; // Add one just to make sure we shift enough

                // Reset the matrix
                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | arpack LM | residual on ...",
                                  largest_eigval);
                hamiltonian_squared = MatVecMPO<double>(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                residual            = initial_mps.get_tensor().real();
                solver_shift.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF, eig::Vecs::ON,
                                  eig::Dephase::OFF, residual.data());
                tools::log->trace("Finding excited state using shifted operator H²-λmax, with λmax = {:.16f} | residual on ... OK: iters {} | counter {} | "
                                  "time {:.3f} | matprod {:.3f}s",
                                  largest_eigval, solver_shift.result.meta.iter, solver_shift.result.meta.counter, solver_shift.result.meta.time_total,
                                  solver_shift.result.meta.time_matprod);
                tools::log->info("arnoldi found: {}", solver_shift.result.meta.arnoldi_found);
                tools::log->info("eigvecs found: {}", solver_shift.result.meta.eigvecsR_found);
                if(not solver_shift.result.meta.arnoldi_found or solver_shift.result.meta.iter > solver_shift.config.maxIter) {
                    // Reset the matrix
                    tools::log->debug("Finding excited state using shifted operator [(H-E)²-λ], with λ = 1.0 | primme SA | residual on ...");
                    hamiltonian_squared = MatVecMPO<double>(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                    residual            = initial_mps.get_tensor().real();
                    solver_primme_sa.eigs(hamiltonian_squared, -1, -1, eig::Ritz::SA, eig::Form::SYMM, eig::Side::R, 1.0, eig::Shinv::OFF, eig::Vecs::ON,
                                          eig::Dephase::OFF, residual.data());
                }

            } else {
                tools::log->warn("Finding excited state using H² with ritz SM");
                solver_shift.eigs(hamiltonian_squared, -1, -1, eig::Ritz::SM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::ON,
                                  eig::Dephase::OFF, residual.data());
            }
        }
    } else if(optType == OptType::CPLX) {
        tools::log->trace("- Generating complex-valued multisite components");
        const auto &           env2      = tensors.get_multisite_var_blk();
        Eigen::Tensor<cplx, 3> env2L     = env2.L;
        Eigen::Tensor<cplx, 3> env2R     = env2.R;
        Eigen::Tensor<cplx, 4> mpo2      = tensors.get_multisite_mpo_squared();
        Eigen::Tensor<cplx, 3> residual  = initial_mps.get_tensor();
        auto                   dims_mpo2 = mpo2.dimensions();

        if(size <= settings::precision::max_size_full_diag) {
            tools::log->debug("Finding excited state using full spectrum of H²");
            auto matrix = tools::finite::opt::internal::get_multisite_hamiltonian_squared_matrix<cplx>(*tensors.model, *tensors.edges);
            solver_full.eig(matrix.data(), matrix.rows());
        } else {
            tools::log->trace("Defining reduced Hamiltonian-squared matrix-vector product");
            MatVecMPO<cplx> hamiltonian_squared(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);

            if(not tensors.model->is_compressed_mpo_squared()) {
                tools::log->trace("Finding largest-magnitude state");
                solver_max.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::OFF,
                                eig::Dephase::OFF);
                auto eigvals_ham_sq = eig::view::get_eigvals<double>(solver_max.result);
                auto largest_eigval = eigvals_ham_sq.cwiseAbs().maxCoeff() + 1.0; // Add one just to make sure we shift enough

                // Reset the matrix
                hamiltonian_squared = MatVecMPO<cplx>(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                tools::log->debug("Finding excited state using shifted operator H²-λmax, with λmax = {:.16f} | residual on ...", largest_eigval);
                solver_shift.eigs(hamiltonian_squared, -1, -1, eig::Ritz::SA, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF, eig::Vecs::ON,
                                  eig::Dephase::OFF, residual.data());
                tools::log->trace("Finding excited state using shifted operator H²-λmax, with λmax = {:.16f} | residual on ... OK: iters {} | counter {} | "
                                  "time {:.3f} | matprod {:.3f}s",
                                  largest_eigval, solver_shift.result.meta.iter, solver_shift.result.meta.counter, solver_shift.result.meta.time_total,
                                  solver_shift.result.meta.time_matprod);

            } else {
                tools::log->warn("Finding excited state using H² with ritz SM");
                solver_shift.eigs(hamiltonian_squared, -1, -1, eig::Ritz::SM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::ON,
                                  eig::Dephase::OFF, residual.data());
            }
        }
    }

    t_eig.toc();

    std::vector<opt_mps> eigvecs_mps;
    eigvecs_mps.emplace_back(initial_mps);
    internal::krylov_extract_solutions(initial_mps, tensors, solver_shift, eigvecs_mps, "shifted", false);
    internal::krylov_extract_solutions(initial_mps, tensors, solver_primme_sa, eigvecs_mps, "primme", false);
    internal::krylov_extract_solutions(initial_mps, tensors, solver_full, eigvecs_mps, "full");

    if(eigvecs_mps.empty())
        return internal::ceres_direct_optimization(tensors, initial_mps, status, optType, OptMode::VARIANCE, OptSpace::DIRECT); // The solver failed

    auto comp_variance = [](const opt_mps &lhs, const opt_mps &rhs) {
        return lhs.get_variance() < rhs.get_variance();
    };
    auto comp_overlap = [](const opt_mps &lhs, const opt_mps &rhs) {
        return lhs.get_overlap() > rhs.get_overlap();
    };

    if(eigvecs_mps.size() >= 2) {
        if(optMode == OptMode::VARIANCE) std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), comp_variance);
        if(optMode == OptMode::OVERLAP) std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), comp_overlap);
    }

    constexpr size_t max_print = settings::debug ? 32 : 8;
    for(const auto &[num, mps] : iter::enumerate(eigvecs_mps)) {
        if(num >= max_print) break;
        internal::reports::krylov_add_entry(mps);
    }
    reports::print_krylov_report();
    std::vector<opt_mps> optimized_mps;

    // Try optimizing the first eigenvectors
    for(const auto &mps : eigvecs_mps) {
        if(&mps != &eigvecs_mps.front()) {
            // Break secondary states that are too far away in energy from the first one
            if(std::abs(mps.get_energy() - eigvecs_mps.front().get_energy()) > 1e2 * eigvecs_mps.front().get_variance()) break;
            // Break secondary states that are too far away in overlap from the first one
            if(mps.get_overlap() < 0.7) break;
        }
        optimized_mps.emplace_back(internal::ceres_direct_optimization(tensors, mps, status, optType, optMode, optSpace));
    }

    // Append all the optimized mps's
    for(const auto &mps : optimized_mps) eigvecs_mps.emplace_back(mps);

    if(eigvecs_mps.size() >= 2) // Sort the results in order of increasing variance (again!)
        std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), comp_variance);
    return eigvecs_mps.front();
    //    if(eigvecs_mps.empty())
    //        return internal::ceres_direct_optimization(tensors, initial_mps, status, optType, OptMode::VARIANCE,
    //                                                   OptSpace::DIRECT); // The optimization strategy failed
    //    else
    //    return eigvecs_mps.front();
    //
}
