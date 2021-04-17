
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_iter.h>
#include <general/nmspc_tensor_extra.h>
#include <math/eig.h>
#include <math/eig/arpack_solver/matrix_product_hamiltonian.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_mps.h>

void tools::finite::opt::extract_solutions(const opt_mps &initial_mps, const class_tensors_finite &tensors, eig::solver &solver,
                                           std::vector<tools::finite::opt::opt_mps> &eigvecs_mps, const std::string &tag) {
    if(solver.result.meta.eigvals_found) {
        auto dims_mps = initial_mps.get_tensor().dimensions();
        auto eigvecs  = eig::view::get_eigvecs<Scalar>(solver.result, eig::Side::R, true);
        auto eigvals  = eig::view::get_eigvals<double>(solver.result, true);
        eigvecs_mps.reserve(eigvecs_mps.size() + static_cast<size_t>(eigvals.size()));
        for(long idx = 0; idx < eigvals.size(); idx++) {
            auto eigvec_i = Textra::asNormalized(Textra::TensorCast(eigvecs.col(idx), dims_mps));
            auto overlap  = std::abs(initial_mps.get_vector().dot(eigvecs.col(idx)));
            auto energy   = tools::finite::measure::energy(eigvec_i, tensors);
            auto eigval   = energy - initial_mps.get_energy_reduced();
            auto variance = tools::finite::measure::energy_variance(eigvec_i, tensors);
            eigvecs_mps.emplace_back(fmt::format("{:<8} eigenvector {}", tag,idx), eigvec_i, tensors.active_sites, eigval, initial_mps.get_energy_reduced(),
                                     variance, overlap, tensors.get_length());
            auto & mps = eigvecs_mps.back();
            mps.set_time(solver.result.meta.time_total);
            mps.set_counter(static_cast<size_t>(solver.result.meta.counter));
            mps.set_iter(static_cast<size_t>(solver.result.meta.iter));
            mps.is_basis_vector = true;
            mps.validate_candidate();
            mps.set_krylov_nev(solver.result.meta.nev_converged);
            mps.set_krylov_ncv(solver.result.meta.ncv);
            mps.set_krylov_tol(solver.result.meta.tol);
            mps.set_krylov_eigval(eigvals(idx));
            mps.set_krylov_ritz(solver.result.meta.ritz);
        }
    }
}

tools::finite::opt::opt_mps tools::finite::opt::internal::arpack_variance_optimization(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                                                       const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                                       OptSpace optSpace) {
    using namespace internal;
    using namespace settings::precision;
    if (not tensors.model->is_reduced()) throw std::runtime_error("arpack_variance_optimization requires energy-reduced MPO²");
    if (tensors.model->is_compressed_mpo_squared()) throw std::runtime_error("arpack_variance_optimization requires non-compressed MPO²");


    auto t_eig    = tools::common::profile::get_default_prof()["t_eig"]->tic_token();
    auto dims_mps = initial_mps.get_tensor().dimensions();
    auto size     = initial_mps.get_tensor().size();

    // NOTE:
    // We can't use the MPO-shift trick here!
    // The lower-left corner of the MPO is NOT where the energy is subtracted on a double-layer MPO,
    // especially when using MPO compression.
    tools::log->debug("Finding excited state: smallest eigenstate of (H-E)² | dims {} = {}", dims_mps, size);

    eig::solver solver_full, solver_shft, solver_max;
    solver_max.config.eigThreshold = 1e-12;
    solver_max.config.eigMaxIter   = 1000;
    solver_max.config.eigMaxNev    = 1;
    solver_max.config.eigMaxNcv    = 16;
    solver_max.config.compress     = true;
    solver_max.setLogLevel(2);

    solver_shft.config.eigThreshold = std::clamp(0.1 * initial_mps.get_variance(), settings::precision::eig_threshold, 1e-10,);
    solver_shft.config.eigMaxIter   = 1000;
    solver_shft.config.eigMaxNev    = 1;
    solver_shft.config.eigMaxNcv    = 128;
    solver_shft.config.compress     = true;
    solver_shft.setLogLevel(2);

    eig::solver solver_shft2;
    solver_shft2.config.eigThreshold = 1e-12;
    solver_shft2.config.eigMaxIter   = 200;
    solver_shft2.config.eigMaxNev    = 1;
    solver_shft2.config.eigMaxNcv    = 512;
    solver_shft2.config.compress     = true;
    solver_shft2.setLogLevel(2);
//
    eig::solver solver_shft3;
    solver_shft3.config.eigThreshold = 1e-12;
    solver_shft3.config.eigMaxIter   = 100;
    solver_shft3.config.eigMaxNev    = 1;
    solver_shft3.config.eigMaxNcv    = 1024;
    solver_shft3.config.compress     = true;
    solver_shft3.setLogLevel(2);
//
//    eig::solver solver_shft4;
//    solver_shft4.config.eigThreshold = 1e-10;
//    solver_shft4.config.eigMaxIter   = 100000;
//    solver_shft4.config.eigMaxNev    = 1;
//    solver_shft4.config.eigMaxNcv    = 512;
//    solver_shft4.config.compress     = true;
//    solver_shft4.setLogLevel(2);
//
//    eig::solver solver_shft5;
//    solver_shft5.config.eigThreshold = 1e-12;
//    solver_shft5.config.eigMaxIter   = 100000;
//    solver_shft5.config.eigMaxNev    = 1;
//    solver_shft5.config.eigMaxNcv    = 128;
//    solver_shft5.config.compress     = true;
//    solver_shft5.setLogLevel(2);
//
//    eig::solver solver_shft6;
//    solver_shft6.config.eigThreshold = 1e-12;
//    solver_shft6.config.eigMaxIter   = 100000;
//    solver_shft6.config.eigMaxNev    = 1;
//    solver_shft6.config.eigMaxNcv    = 256;
//    solver_shft6.config.compress     = true;
//    solver_shft6.setLogLevel(2);
//
//    eig::solver solver_shft7;
//    solver_shft7.config.eigThreshold = 1e-12;
//    solver_shft7.config.eigMaxIter   = 100000;
//    solver_shft7.config.eigMaxNev    = 1;
//    solver_shft7.config.eigMaxNcv    = 512;
//    solver_shft7.config.compress     = true;
//    solver_shft7.setLogLevel(2);

//    eig::solver solver_shft8;
//    solver_shft8.config.eigThreshold = 1e-12;
//    solver_shft8.config.eigMaxIter   = 100000;
//    solver_shft8.config.eigMaxNev    = 1;
//    solver_shft8.config.eigMaxNcv    = 16;
//    solver_shft8.config.compress     = true;
//    solver_shft8.setLogLevel(2);
//
//
//    eig::solver solver_shft9;
//    solver_shft9.config.eigThreshold = 1e-12;
//    solver_shft9.config.eigMaxIter   = 100000;
//    solver_shft9.config.eigMaxNev    = 1;
//    solver_shft9.config.eigMaxNcv    = 32;
//    solver_shft9.config.compress     = true;
//    solver_shft9.setLogLevel(2);
//
//    eig::solver solver_shft10;
//    solver_shft10.config.eigThreshold = 1e-12;
//    solver_shft10.config.eigMaxIter   = 100000;
//    solver_shft10.config.eigMaxNev    = 1;
//    solver_shft10.config.eigMaxNcv    = 64;
//    solver_shft10.config.compress     = true;
//    solver_shft10.setLogLevel(2);


    if(optType == OptType::REAL) {
        tools::log->trace("- Generating real-valued multisite components");
        const auto &             env2      = tensors.get_multisite_var_blk();
        Eigen::Tensor<double, 3> env2L     = env2.L.real();
        Eigen::Tensor<double, 3> env2R     = env2.R.real();
        Eigen::Tensor<double, 4> mpo2      = tensors.get_multisite_mpo_squared().real();
        Eigen::Tensor<double, 3> residual  = initial_mps.get_tensor().real();
        auto                     dims_mpo2 = mpo2.dimensions();

        if(size <= 128) {
//        if(size <= settings::precision::max_size_full_diag) {
            tools::log->trace("Full diagonalization of (H-E)²", dims_mps, size);
            auto matrix = tools::finite::opt::internal::get_multisite_hamiltonian_squared_matrix<double>(*tensors.model, *tensors.edges);
            solver_full.eig(matrix.data(), matrix.rows());
        } else {
            tools::log->trace("Defining reduced Hamiltonian-squared matrix-vector product");
            MatrixProductHamiltonian<double> hamiltonian_squared(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);

            if(not tensors.model->is_compressed_mpo_squared()) {
                tools::log->trace("Finding largest-magnitude state");
                solver_max.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF,
                                    eig::Vecs::OFF, eig::Dephase::OFF);
                auto eigvals_ham_sq = eig::view::get_eigvals<double>(solver_max.result);
                auto largest_eigval = eigvals_ham_sq.cwiseAbs().maxCoeff() + 1.0; // Add one just to make sure we shift enough

                // Reset the matrix
                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
                solver_shft.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                                 eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                tools::log->trace("Finding excited state using shifted operator H²-λmax, with λmax = {:.16f} | residual on ... OK: iters {} | counter {} | "
                                  "time {:.3f} | matprod {:.3f}s",
                                  largest_eigval, solver_shft.result.meta.iter, solver_shft.result.meta.counter, solver_shft.result.meta.time_total,
                                  solver_shft.result.meta.time_matprod);

                // Reset the matrix
//                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
//                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//                residual  = initial_mps.get_tensor().real();
//                solver_shft2.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
//                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
//                // Reset the matrix
//                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
//                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//                residual  = initial_mps.get_tensor().real();
//                solver_shft3.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
//                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
//                // Reset the matrix
//                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
//                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//                residual  = initial_mps.get_tensor().real();
//                solver_shft4.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
//                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
//                // Reset the matrix
//                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
//                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//                residual  = initial_mps.get_tensor().real();
//                solver_shft5.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
//                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
//                // Reset the matrix
//                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
//                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//                residual  = initial_mps.get_tensor().real();
//                solver_shft6.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
//                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
//                // Reset the matrix
//                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
//                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//                residual  = initial_mps.get_tensor().real();
//                solver_shft7.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
//                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
//                // Reset the matrix
//                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
//                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//                residual  = initial_mps.get_tensor().real();
//                solver_shft8.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
//                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
//                // Reset the matrix
//                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
//                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//                residual  = initial_mps.get_tensor().real();
//                solver_shft9.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
//                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
//                // Reset the matrix
//                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
//                hamiltonian_squared = MatrixProductHamiltonian<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
//                residual  = initial_mps.get_tensor().real();
//                solver_shft10.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
//                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
#pragma message "print also the actual H² eigval"
            } else {
                tools::log->warn("Finding excited state using H² with ritz SM");
                solver_shft.eigs(hamiltonian_squared, -1, -1, eig::Ritz::SM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF,
                                 eig::Vecs::ON, eig::Dephase::OFF, residual.data());
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
            MatrixProductHamiltonian<cplx> hamiltonian_squared(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);

            if(not tensors.model->is_compressed_mpo_squared()) {
                tools::log->trace("Finding largest-magnitude state");
                solver_max.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF,
                                    eig::Vecs::OFF, eig::Dephase::OFF);
                auto eigvals_ham_sq = eig::view::get_eigvals<double>(solver_max.result);
                auto largest_eigval = eigvals_ham_sq.cwiseAbs().maxCoeff() + 1.0; // Add one just to make sure we shift enough

                // Reset the matrix
                hamiltonian_squared = MatrixProductHamiltonian<cplx> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                tools::log->debug("Finding excited state using shifted operator H²-λmax, with λmax = {:.16f} | residual on ...", largest_eigval);
                solver_shft.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                                 eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                tools::log->trace("Finding excited state using shifted operator H²-λmax, with λmax = {:.16f} | residual on ... OK: iters {} | counter {} | "
                                  "time {:.3f} | matprod {:.3f}s",
                                  largest_eigval, solver_shft.result.meta.iter, solver_shft.result.meta.counter, solver_shft.result.meta.time_total,
                                  solver_shft.result.meta.time_matprod);

            } else {
                tools::log->warn("Finding excited state using H² with ritz SM");
                solver_shft.eigs(hamiltonian_squared, -1, -1, eig::Ritz::SM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF,
                                 eig::Vecs::ON, eig::Dephase::OFF, residual.data());
            }
        }
    }

    t_eig.toc();

    std::vector<opt_mps> eigvecs_mps;
    eigvecs_mps.emplace_back(initial_mps);
    extract_solutions(initial_mps, tensors, solver_shft, eigvecs_mps, "shifted");
    extract_solutions(initial_mps, tensors, solver_shft2, eigvecs_mps, "shifted2");
    extract_solutions(initial_mps, tensors, solver_shft3, eigvecs_mps, "shifted3");
//    extract_solutions(initial_mps, tensors, solver_shft4, eigvecs_mps, "shifted4");
//    extract_solutions(initial_mps, tensors, solver_shft5, eigvecs_mps, "shifted5");
//    extract_solutions(initial_mps, tensors, solver_shft6, eigvecs_mps, "shifted6");
//    extract_solutions(initial_mps, tensors, solver_shft7, eigvecs_mps, "shifted7");
//    extract_solutions(initial_mps, tensors, solver_shft8, eigvecs_mps, "shifted8");
//    extract_solutions(initial_mps, tensors, solver_shft9, eigvecs_mps, "shifted9");
//    extract_solutions(initial_mps, tensors, solver_shft10, eigvecs_mps, "shifted10");
    extract_solutions(initial_mps, tensors, solver_full, eigvecs_mps, "full");

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

    constexpr size_t max_print = settings::debug ? 32 : 20;
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
    for(const auto & mps: optimized_mps) eigvecs_mps.emplace_back(mps);

    if(eigvecs_mps.size() >= 2) // Sort the results in order of increasing variance (again!)
        std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), comp_variance);

    if(eigvecs_mps.empty())
        return internal::ceres_direct_optimization(tensors, initial_mps, status, optType, OptMode::VARIANCE,
                                                   OptSpace::DIRECT); // The optimization strategy failed
    else
        return eigvecs_mps.front();
}
