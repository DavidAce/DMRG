
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <math/eig.h>
#include <math/eig/arpack_solver/matrix_product_hamiltonian.h>
#include <math/linalg/tensor.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_mps.h>

tools::finite::opt::opt_mps tools::finite::opt::internal::arpack_variance_optimization(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                                                       const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                                       OptSpace optSpace) {
    using namespace internal;
    using namespace settings::precision;
    auto t_eig  = tools::common::profile::get_default_prof()["t_eig"]->tic_token();
    auto dims_mps = initial_mps.get_tensor().dimensions();
    auto size     = initial_mps.get_tensor().size();
    auto nev      = std::clamp<eig::size_type>(size / 16, std::min(size, 4l), 8);
    auto ncv      = static_cast<eig::size_type>(16);
    tools::log->info("Excited state variance optimization with arpack | dims {} = {}", dims_mps, size);

    tools::log->trace("Defining eigenvalue solver");
    // NOTE:
    // We can't use the MPO-shift trick here!
    // The lower-left corner of the MPO is NOT where the energy is subtracted on a double-layer MPO,
    // especially when using MPO compression.

    eig::solver solver_full, solver_part;
    solver_part.config.eigThreshold = 1e-10;
    solver_part.config.eigMaxIter   = 400;
    if(optType == OptType::REAL) {
        tools::log->trace("- Generating real-valued multisite components");
        const auto &             env2     = tensors.get_multisite_var_blk();
        Eigen::Tensor<double, 3> env2L    = env2.L.real();
        Eigen::Tensor<double, 3> env2R    = env2.R.real();
        Eigen::Tensor<double, 4> mpo2     = tensors.get_multisite_mpo_squared().real();
        auto                     dims_mpo = mpo2.dimensions();
        tools::log->trace("Defining reduced Hamiltonian-squared matrix-vector product");
        MatrixProductHamiltonian<double> matrix(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo);

        tools::log->trace("Finding excited state");
        Eigen::Tensor<double, 3> residual = initial_mps.get_tensor().real();
        solver_part.eigs(matrix, nev, ncv, eig::Ritz::SA, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF, residual.data());
    }
    if(optType == OptType::CPLX) {
        tools::log->trace("- Generating cplx-valued multisite components");
        const auto &env2     = tensors.get_multisite_var_blk();
        const auto &env2L    = env2.L;
        const auto &env2R    = env2.R;
        const auto &mpo2     = tensors.get_multisite_mpo_squared();
        const auto  dims_mpo = mpo2.dimensions();
        tools::log->trace("Defining reduced Hamiltonian-squared matrix-vector product");
        MatrixProductHamiltonian<Scalar> matrix(env2L.data(), env2R.data(), mpo2.data(), dims_mps, dims_mpo);

        // Since we use reduced energy mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
        // would otherwise cause trouble for Arpack. This equates to subtracting sigma * identity from the bottom corner of the mpo.
        // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
        Eigen::Tensor<Scalar, 3> residual = initial_mps.get_tensor();
        tools::log->trace("Finding excited state");
        solver_part.eigs(matrix, nev, ncv,  eig::Ritz::SR, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF,residual.data());
    }
    t_eig.toc();

    std::vector<opt_mps> eigvecs_mps;
    if(solver_part.result.meta.eigvals_found){
        auto eigvecs = eig::view::get_eigvecs<Scalar>(solver_part.result);
        auto eigvals = eig::view::get_eigvals<double>(solver_part.result);
        eigvecs_mps.reserve(eigvecs_mps.size() + static_cast<size_t>(eigvals.size()));
        for(long idx = 0; idx < eigvals.size(); idx++) {
            auto eigvec_i = Textra::TensorCast(eigvecs.col(idx), dims_mps);
            auto overlap  = std::abs(initial_mps.get_vector().dot(eigvecs.col(idx)));
            auto energy   = tools::finite::measure::energy(eigvec_i, tensors);
            auto eigval   = energy - initial_mps.get_energy_reduced();
            auto variance = tools::finite::measure::energy_variance(eigvec_i, tensors);
            eigvecs_mps.emplace_back(fmt::format("eigenvector {}", idx), eigvec_i, tensors.active_sites, eigval, initial_mps.get_energy_reduced(), variance,
                                     overlap, tensors.get_length());

            eigvecs_mps.back().set_time(tools::common::profile::get_default_prof()["t_eig"]->get_last_interval());
            eigvecs_mps.back().set_counter(static_cast<size_t>(solver_part.result.meta.counter));
            eigvecs_mps.back().set_iter(static_cast<size_t>(solver_part.result.meta.iter));
            eigvecs_mps.back().is_basis_vector = true;
            eigvecs_mps.back().validate_candidate();
            tools::log->info("eigenvector {:<2} overlap {:.16f} | E {:>20.16f} | E/L {:>20.16f} | eigval {:>20.16f} {:>20.16f} | Er {:>20.16f} | variance: "
                             "{:>20.16f} | iters {}",
                             idx, eigvecs_mps.back().get_overlap(), eigvecs_mps.back().get_energy(), eigvecs_mps.back().get_energy_per_site(),
                             eigvecs_mps.back().get_eigval(), eigvals(idx), eigvecs_mps.back().get_energy_reduced(), std::log10(eigvecs_mps.back().get_variance()),
                             eigvecs_mps.back().get_iter());
        }
    }

    if(solver_full.result.meta.eigvals_found){
        auto eigvecs = eig::view::get_eigvecs<Scalar>(solver_full.result);
        auto eigvals = eig::view::get_eigvals<double>(solver_full.result);
        eigvecs_mps.reserve(eigvecs_mps.size() + static_cast<size_t>(eigvals.size()));
        for(long idx = 0; idx < std::min(8l,eigvals.size()); idx++) {
            auto eigvec_i = Textra::TensorCast(eigvecs.col(idx), dims_mps);
            auto overlap  = std::abs(initial_mps.get_vector().dot(eigvecs.col(idx)));
            auto energy   = tools::finite::measure::energy(eigvec_i, tensors);
            auto eigval   = energy - initial_mps.get_energy_reduced();
            auto variance = tools::finite::measure::energy_variance(eigvec_i, tensors);
            eigvecs_mps.emplace_back(fmt::format("eigenvector {}", idx), eigvec_i, tensors.active_sites, eigval, initial_mps.get_energy_reduced(), variance,
                                     overlap, tensors.get_length());

            eigvecs_mps.back().set_time(tools::common::profile::get_default_prof()["t_eig"]->get_last_interval());
            eigvecs_mps.back().set_counter(static_cast<size_t>(solver_full.result.meta.counter));
            eigvecs_mps.back().set_iter(static_cast<size_t>(solver_full.result.meta.iter));
            eigvecs_mps.back().is_basis_vector = true;
            eigvecs_mps.back().validate_candidate();
            tools::log->info("eigenvector {:<2} overlap {:.16f} | E {:>20.16f} | E/L {:>20.16f} | eigval {:>20.16f} {:>20.16f} | Er {:>20.16f} | variance: "
                             "{:>20.16f} | iters {}",
                             idx, eigvecs_mps.back().get_overlap(), eigvecs_mps.back().get_energy(), eigvecs_mps.back().get_energy_per_site(),
                             eigvecs_mps.back().get_eigval(), eigvals(idx), eigvecs_mps.back().get_energy_reduced(), std::log10(eigvecs_mps.back().get_variance()),
                             eigvecs_mps.back().get_iter());
        }
    }


    if(eigvecs_mps.size() >= 2){
        std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), [](const opt_mps &lhs, const opt_mps &rhs){ return lhs.get_variance() < rhs.get_variance();});
        auto super_vector = (eigvecs_mps[0].get_vector() + eigvecs_mps[1].get_vector()).normalized();
        auto super_tensor = Textra::TensorCast(super_vector, dims_mps);
        auto overlap  = std::abs(initial_mps.get_vector().dot(super_vector));
        auto energy   = tools::finite::measure::energy(super_tensor, tensors);
        auto eigval   = energy - initial_mps.get_energy_reduced();
        auto variance = tools::finite::measure::energy_variance(super_tensor, tensors);
        auto super_mps = opt_mps(fmt::format("superposition 0+1"), super_tensor, tensors.active_sites, eigval, initial_mps.get_energy_reduced(), variance,
                           overlap, tensors.get_length());

    }

    auto current_level = tools::log->level();
    tools::log->set_level(spdlog::level::trace);


    std::vector<opt_mps> optimized_mps;
    if(not eigvecs_mps.empty()){
        // Sort the results in order of increasing variance
        std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), [](const opt_mps &lhs, const opt_mps &rhs){ return lhs.get_variance() < rhs.get_variance();});

        // Try optimizing the first few eigenvectors
        for(const auto & mps : eigvecs_mps){
            if(&mps != &eigvecs_mps.front()){
                // Only break secondary states, the first always gets optimized
                if(std::abs(mps.get_energy() - eigvecs_mps.front().get_energy()) > 1e2 * eigvecs_mps.front().get_variance()) break; // Don't bother if the energy is too far away
            }

            optimized_mps.emplace_back(internal::ceres_direct_optimization(tensors, mps, status, optType, optMode, optSpace));
        }
        if(optimized_mps.size() >= 2) // Sort the results in order of increasing variance (again!)
            std::sort(optimized_mps.begin(), optimized_mps.end(), [](const opt_mps &lhs, const opt_mps &rhs){ return lhs.get_variance() < rhs.get_variance();});

    }
    tools::finite::opt::internal::reports::print_bfgs_report();
    tools::finite::opt::internal::reports::print_time_report();
    tools::log->set_level(current_level);
    if(optimized_mps.empty()) return internal::ceres_direct_optimization(tensors, initial_mps, status, optType, optMode, optSpace);
    else return optimized_mps.front();


}
