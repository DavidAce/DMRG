
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_iter.h>
#include <general/nmspc_tensor_extra.h>
#include <math/eig_primme.h>
#include <math/eig_primme/primme_solver/matrix_product_hamiltonian.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_mps.h>

tools::finite::opt::opt_mps tools::finite::opt::internal::primme_variance_optimization(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                                                       const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                                       OptSpace optSpace) {
    using namespace internal;
    using namespace settings::precision;
    if(not tensors.model->is_reduced()) throw std::runtime_error("primme_variance_optimization requires energy-reduced MPO²");
    if(tensors.model->is_compressed_mpo_squared()) throw std::runtime_error("primme_variance_optimization requires non-compressed MPO²");

    auto tic      = tid::tic_scope("var");
    auto dims_mps = initial_mps.get_tensor().dimensions();
    auto size     = initial_mps.get_tensor().size();

    // NOTE:
    // We can't use the MPO-shift trick here!
    // The lower-left corner of the MPO is NOT where the energy is subtracted on a double-layer MPO,
    // especially when using MPO compression.
    tools::log->debug("Finding excited state: smallest eigenstate of (H-E)² | dims {} = {}", dims_mps, size);

    eig::primme::solver solver_full, solver_shft, solver_max;
    solver_max.config.tol      = settings::precision::eig_tolerance;
    solver_max.config.maxIter  = 600;
    solver_max.config.maxNev   = 1;
    solver_max.config.maxNcv   = 16;
    solver_max.config.compress = true;
    solver_max.setLogLevel(2);

    solver_shft.config.tol          = settings::precision::eig_tolerance;
    solver_shft.config.maxIter      = 1000;
    solver_shft.config.maxNev       = 1;
    solver_shft.config.maxNcv       = 64;
    solver_shft.config.compress     = true;
    solver_shft.config.ncv_x_factor = 8;
    solver_shft.config.iter_ncv_x   = {250, 200, 150, 100, 50};
    solver_shft.config.time_tol_x10 = {4 * 3600, 2 * 3600, 30 * 60, 10 * 60};

    solver_shft.setLogLevel(2);
    //
    //    eig::solver solver_shft2;
    //    solver_shft2.config.tol = 1e-12;
    //    solver_shft2.config.maxIter   = 200;
    //    solver_shft2.config.maxNev    = 1;
    //    solver_shft2.config.maxNcv    = 512;
    //    solver_shft2.config.compress     = true;
    //    solver_shft2.setLogLevel(2);
    ////
    //    eig::solver solver_shft3;
    //    solver_shft3.config.tol = 1e-12;
    //    solver_shft3.config.maxIter   = 100;
    //    solver_shft3.config.maxNev    = 1;
    //    solver_shft3.config.maxNcv    = 1024;
    //    solver_shft3.config.compress     = true;
    //    solver_shft3.setLogLevel(2);
    //
    //    eig::solver solver_shft4;
    //    solver_shft4.config.tol = 1e-10;
    //    solver_shft4.config.maxIter   = 100000;
    //    solver_shft4.config.maxNev    = 1;
    //    solver_shft4.config.maxNcv    = 512;
    //    solver_shft4.config.compress     = true;
    //    solver_shft4.setLogLevel(2);
    //
    //    eig::solver solver_shft5;
    //    solver_shft5.config.tol = 1e-12;
    //    solver_shft5.config.maxIter   = 100000;
    //    solver_shft5.config.maxNev    = 1;
    //    solver_shft5.config.maxNcv    = 128;
    //    solver_shft5.config.compress     = true;
    //    solver_shft5.setLogLevel(2);
    //
    //    eig::solver solver_shft6;
    //    solver_shft6.config.tol = 1e-12;
    //    solver_shft6.config.maxIter   = 100000;
    //    solver_shft6.config.maxNev    = 1;
    //    solver_shft6.config.maxNcv    = 256;
    //    solver_shft6.config.compress     = true;
    //    solver_shft6.setLogLevel(2);
    //
    //    eig::solver solver_shft7;
    //    solver_shft7.config.tol = 1e-12;
    //    solver_shft7.config.maxIter   = 100000;
    //    solver_shft7.config.maxNev    = 1;
    //    solver_shft7.config.maxNcv    = 512;
    //    solver_shft7.config.compress     = true;
    //    solver_shft7.setLogLevel(2);

    //    eig::solver solver_shft8;
    //    solver_shft8.config.tol = 1e-12;
    //    solver_shft8.config.maxIter   = 100000;
    //    solver_shft8.config.maxNev    = 1;
    //    solver_shft8.config.maxNcv    = 16;
    //    solver_shft8.config.compress     = true;
    //    solver_shft8.setLogLevel(2);
    //
    //
    //    eig::solver solver_shft9;
    //    solver_shft9.config.tol = 1e-12;
    //    solver_shft9.config.maxIter   = 100000;
    //    solver_shft9.config.maxNev    = 1;
    //    solver_shft9.config.maxNcv    = 32;
    //    solver_shft9.config.compress     = true;
    //    solver_shft9.setLogLevel(2);
    //
    //    eig::solver solver_shft10;
    //    solver_shft10.config.tol = 1e-12;
    //    solver_shft10.config.maxIter   = 100000;
    //    solver_shft10.config.maxNev    = 1;
    //    solver_shft10.config.maxNcv    = 64;
    //    solver_shft10.config.compress     = true;
    //    solver_shft10.setLogLevel(2);

    if(optType == OptType::REAL) {
        tools::log->trace("- Generating real-valued multisite components");
        const auto              &env2      = tensors.get_multisite_var_blk();
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
                solver_max.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::OFF,
                                eig::Dephase::OFF);
                auto eigvals_ham_sq = eig::view::get_eigvals<double>(solver_max.result);
                auto largest_eigval = eigvals_ham_sq.cwiseAbs().maxCoeff() + 1.0; // Add one just to make sure we shift enough

                // Reset the matrix
                hamiltonian_squared = MatrixProductHamiltonian<double>(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...", largest_eigval);
                solver_shft.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF, eig::Vecs::ON,
                                 eig::Dephase::OFF, residual.data());
                tools::log->trace("Finding excited state using shifted operator H²-λmax, with λmax = {:.16f} | residual on ... OK: iters {} | counter {} | "
                                  "time {:.3f} | matprod {:.3f}s",
                                  largest_eigval, solver_shft.result.meta.iter, solver_shft.result.meta.counter, solver_shft.result.meta.time_total,
                                  solver_shft.result.meta.time_matprod);

                // Reset the matrix
                //                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...",
                //                largest_eigval); hamiltonian_squared = MatVecMPO<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                //                residual  = initial_mps.get_tensor().real();
                //                solver_shft2.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                //                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                //                // Reset the matrix
                //                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...",
                //                largest_eigval); hamiltonian_squared = MatVecMPO<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                //                residual  = initial_mps.get_tensor().real();
                //                solver_shft3.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                //                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                //                // Reset the matrix
                //                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...",
                //                largest_eigval); hamiltonian_squared = MatVecMPO<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                //                residual  = initial_mps.get_tensor().real();
                //                solver_shft4.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                //                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                //                // Reset the matrix
                //                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...",
                //                largest_eigval); hamiltonian_squared = MatVecMPO<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                //                residual  = initial_mps.get_tensor().real();
                //                solver_shft5.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                //                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                //                // Reset the matrix
                //                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...",
                //                largest_eigval); hamiltonian_squared = MatVecMPO<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                //                residual  = initial_mps.get_tensor().real();
                //                solver_shft6.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                //                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                //                // Reset the matrix
                //                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...",
                //                largest_eigval); hamiltonian_squared = MatVecMPO<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                //                residual  = initial_mps.get_tensor().real();
                //                solver_shft7.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                //                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                //                // Reset the matrix
                //                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...",
                //                largest_eigval); hamiltonian_squared = MatVecMPO<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                //                residual  = initial_mps.get_tensor().real();
                //                solver_shft8.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                //                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                //                // Reset the matrix
                //                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...",
                //                largest_eigval); hamiltonian_squared = MatVecMPO<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                //                residual  = initial_mps.get_tensor().real();
                //                solver_shft9.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                //                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
                //                // Reset the matrix
                //                tools::log->debug("Finding excited state using shifted operator [(H-E)²-λmax], with λmax = {:.16f} | residual on ...",
                //                largest_eigval); hamiltonian_squared = MatVecMPO<double> (env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                //                residual  = initial_mps.get_tensor().real();
                //                solver_shft10.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF,
                //                                  eig::Vecs::ON, eig::Dephase::OFF, residual.data());
            } else {
                tools::log->warn("Finding excited state using H² with ritz SM");
                solver_shft.eigs(hamiltonian_squared, -1, -1, eig::Ritz::SM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::ON,
                                 eig::Dephase::OFF, residual.data());
            }
        }
    } else if(optType == OptType::CPLX) {
        tools::log->trace("- Generating complex-valued multisite components");
        const auto            &env2      = tensors.get_multisite_var_blk();
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
                solver_max.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::OFF,
                                eig::Dephase::OFF);
                auto eigvals_ham_sq = eig::view::get_eigvals<double>(solver_max.result);
                auto largest_eigval = eigvals_ham_sq.cwiseAbs().maxCoeff() + 1.0; // Add one just to make sure we shift enough

                // Reset the matrix
                hamiltonian_squared = MatrixProductHamiltonian<cplx>(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo2);
                tools::log->debug("Finding excited state using shifted operator H²-λmax, with λmax = {:.16f} | residual on ...", largest_eigval);
                solver_shft.eigs(hamiltonian_squared, -1, -1, eig::Ritz::LM, eig::Form::SYMM, eig::Side::R, largest_eigval, eig::Shinv::OFF, eig::Vecs::ON,
                                 eig::Dephase::OFF, residual.data());
                tools::log->trace("Finding excited state using shifted operator H²-λmax, with λmax = {:.16f} | residual on ... OK: iters {} | counter {} | "
                                  "time {:.3f} | matprod {:.3f}s",
                                  largest_eigval, solver_shft.result.meta.iter, solver_shft.result.meta.counter, solver_shft.result.meta.time_total,
                                  solver_shft.result.meta.time_matprod);

            } else {
                tools::log->warn("Finding excited state using H² with ritz SM");
                solver_shft.eigs(hamiltonian_squared, -1, -1, eig::Ritz::SM, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::ON,
                                 eig::Dephase::OFF, residual.data());
            }
        }
    }

    t_eig.toc();

    std::vector<opt_mps> eigvecs_mps;
    eigvecs_mps.emplace_back(initial_mps);
    internal::primme_extract_solutions(initial_mps, tensors, solver_shft, eigvecs_mps, "shifted", false);
    //    init::primme_extract_solutions(initial_mps, tensors, solver_shft2, eigvecs_mps, "shifted2");
    //    init::primme_extract_solutions(initial_mps, tensors, solver_shft3, eigvecs_mps, "shifted3");
    //    extract_solutions(initial_mps, tensors, solver_shft4, eigvecs_mps, "shifted4");
    //    extract_solutions(initial_mps, tensors, solver_shft5, eigvecs_mps, "shifted5");
    //    extract_solutions(initial_mps, tensors, solver_shft6, eigvecs_mps, "shifted6");
    //    extract_solutions(initial_mps, tensors, solver_shft7, eigvecs_mps, "shifted7");
    //    extract_solutions(initial_mps, tensors, solver_shft8, eigvecs_mps, "shifted8");
    //    extract_solutions(initial_mps, tensors, solver_shft9, eigvecs_mps, "shifted9");
    //    extract_solutions(initial_mps, tensors, solver_shft10, eigvecs_mps, "shifted10");
    internal::primme_extract_solutions(initial_mps, tensors, solver_full, eigvecs_mps, "full");

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
        internal::reports::primme_add_entry(mps);
    }
    reports::print_primme_report();
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
    //        return init::ceres_direct_optimization(tensors, initial_mps, status, optType, OptMode::VARIANCE,
    //                                                   OptSpace::DIRECT); // The optimization strategy failed
    //    else
    //    return eigvecs_mps.front();
    //
}
