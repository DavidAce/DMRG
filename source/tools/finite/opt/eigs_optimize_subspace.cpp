#include "math/tenx.h"
// -- (textra first)
#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/bfgs_callback.h"
#include "tools/finite/opt/bfgs_subspace_functor.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include <ceres/gradient_problem.h>

using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;

opt_mps tools::finite::opt::internal::eigs_optimize_subspace(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status,
                                                             OptMeta &meta) {
    tools::log->trace("Optimizing subspace");
    auto t_sub = tid::tic_scope("subspace");
    initial_mps.validate_initial_mps();
    if(meta.optMode != OptMode::SUBSPACE)
        throw except::runtime_error("Wrong optimization mode [{}]. Expected [{}]", enum2sv(meta.optMode), enum2sv(OptMode::SUBSPACE));

    /*
     * Subspace optimization
     *
     *
     * In subspace optimization we consider a local set of sites l of the L-site system,
     * usually this corresponds to l = 1 or up to l = 8 adjacent sites.
     *
     * A subspace in this context is a truncated basis, i.e. a small subset of eigenvectors
     * of the local "effective" Hamiltonian. The subspace is a set of k eigenvectors [x]
     * which have significant overlap with the current state |y⟩, i.e. we define k such that
     *
     *          ε > 1 - Σ_i^k |⟨x_i|y⟩|²,
     *
     * where ε is a small number that controls error of the truncation. A value ε ~1e-10 is
     * reasonable. Note that the truncation implies that k is smaller than the local
     * Hilbert space dimension.
     *
     * After having found a good subspace, the objective is to find a linear combination
     * of eigenvectors which minimizes the energy variance.
     *
     * It is worth noting some observations. Let {x} be the set of all eigenvectors
     * to the local effective Hamiltonian H_local.
     * Then, when the DMRG process is fully converged:
     *      - only one x_i has overlap ⟨x_i|y⟩ = 1
     *        Since the sum of all overlaps must add to 1, the rest have <x_j|y> = 0 when i != j.
     *      - This x is also the one that minimizes the energy variance.
     *
     * However, before the DMRG process has converged this is not true. Instead:
     *      - we have ⟨x_i|y⟩ > 0 for several i.
     *      - a linear combination of several x can have lower variance than any
     *        single x.
     *
     * Fully diagonalizing H_local yields all K eigenvectors {x}, but if H_local is too big this operation
     * becomes prohibitively expensive. Instead we resort to finding a subset with k << K eigenvectors [x],
     * whose eigenvalues are the k energies closest to the current energy. Usually the eigenvectors
     * which have some overlap ⟨x_i|y⟩ > 0 are found in the subset [x] if k is large enough.
     *
     * Subspace optimization steps
     *
     * Step 1)  Find a subspace [x], i.e. take a set of k eigenvectors of the local effective Hamiltonian.
     *          Empirically, eigenvectors whose eigenvalues (energy) are closest to the current energy,
     *          tend to have nonzero overlap with the current vector |y⟩.
     *          On an iterative solver we keep increasing "nev" (number of requested eigenvectors) until
     *          the subspace error ε is small enough.
     *          If any eigenvectors have to removed, (e.g. due to memory/performance constraints),
     *          then sort the eigenvectors in order of decreasing overlap ⟨x_i|y⟩, and start deleting
     *          from the end.
     *
     * Step 2)  Project the squared effective K*K Hamiltonian, down to the k*k subspace, H².
     *          Using BFGS, find the linear combination |w⟩ of eigenvectors that minimizes the variance.
     *
     *              min_w Var H = ⟨H²⟩ - ⟨H⟩² = ⟨w|H²|w⟩ - ⟨E⟩²
     *
     *          where E are the energy eigenvalues from step 1.
     *
     */

    // Handy references
    const auto &model = *tensors.model;
    const auto &edges = *tensors.edges;

    /*
     *  Step 1) Find the subspace.
     *  The subspace is a set of eigenstates obtained from full or partial diagonalization
     */

    std::vector<opt_mps> subspace;
    switch(meta.optType) {
        case OptType::CPLX: subspace = internal::subspace::find_subspace<cplx>(tensors, settings::precision::target_subspace_error, meta); break;
        case OptType::REAL: subspace = internal::subspace::find_subspace<real>(tensors, settings::precision::target_subspace_error, meta); break;
    }

    tools::log->trace("Subspace found with {} eigenvectors", subspace.size());

    /*
     * Filter the eigenvectors
     *
     */

    internal::subspace::filter_subspace(subspace, settings::precision::max_subspace_size);

    /*
     *
     * Step 2) Optimize variance in the subspace of k eigenvectors
     *
     */

    // We need the eigenvalues in a convenient format as well
    auto eigvals = internal::subspace::get_eigvals(subspace);

    Eigen::VectorXcd subspace_vector = internal::subspace::get_vector_in_subspace(subspace, initial_mps.get_vector());
    long             subspace_size   = subspace_vector.size();
    reports::bfgs_add_entry("Subspace", "init", initial_mps, subspace_size);
    tools::log->trace("Starting BFGS optimization of subspace");

    /*
     *
     *  Start the BFGS optimization process for the subspace
     *
     */
    opt_mps optimized_mps;
    optimized_mps.set_sites(initial_mps.get_sites());
    optimized_mps.set_length(initial_mps.get_length());
    optimized_mps.set_energy_shift(initial_mps.get_energy_shift());

    opt_mps subspace_mps;
    auto    options = internal::bfgs_default_options;
    auto    summary = ceres::GradientProblemSolver::Summary();
    auto    t_bfgs  = tid::tic_scope("bfgs");
    switch(meta.optType) {
        case OptType::CPLX: {
            auto  H2_subspace               = subspace::get_hamiltonian_squared_in_subspace<cplx>(model, edges, subspace);
            auto  subspace_vector_as_2xreal = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(subspace_vector.data()), 2 * subspace_vector.size());
            auto *functor                   = new bfgs_subspace_functor<cplx>(tensors, status, H2_subspace, eigvals);
            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running BFGS subspace cplx");
            ceres::Solve(options, problem, subspace_vector_as_2xreal.data(), &summary);
            subspace_vector =
                Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<cplx *>(subspace_vector_as_2xreal.data()), subspace_vector_as_2xreal.size() / 2).normalized();
            // Copy the results from the functor
            optimized_mps.set_name("bfgs cplx");
            optimized_mps.set_tensor(subspace::get_vector_in_fullspace(subspace, subspace_vector), initial_mps.get_tensor().dimensions());
            optimized_mps.set_mv(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_grad_max(functor->get_max_grad_norm());
            break;
        }
        case OptType::REAL: {
            Eigen::MatrixXd   H2_subspace          = subspace::get_hamiltonian_squared_in_subspace<real>(model, edges, subspace);
            Eigen::VectorXd   subspace_vector_cast = subspace_vector.real(); // Copy into a temporary
            auto             *functor              = new bfgs_subspace_functor<real>(tensors, status, H2_subspace, eigvals);
            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running BFGS subspace real");
            ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
            subspace_vector = subspace_vector_cast.normalized().cast<cplx>();
            optimized_mps.set_name("bfgs real");
            optimized_mps.set_tensor(subspace::get_vector_in_fullspace(subspace, subspace_vector), initial_mps.get_tensor().dimensions());
            optimized_mps.set_mv(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_grad_max(functor->get_max_grad_norm());
            break;
        }
    }
    reports::time_add_entry();
    t_bfgs.toc();
    // Copy and set the rest of the tensor metadata
    optimized_mps.normalize();
    optimized_mps.set_energy(tools::finite::measure::energy(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_variance(tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_overlap(std::abs(initial_mps.get_vector().dot(optimized_mps.get_vector())));
    optimized_mps.set_eigs_resid(tools::finite::measure::residual_norm(optimized_mps.get_tensor(), tensors.get_multisite_mpo_squared(),
                                                                       tensors.get_multisite_env_var_blk().L, tensors.get_multisite_env_var_blk().R));
    optimized_mps.set_iter(summary.iterations.size());
    optimized_mps.set_time(summary.total_time_in_seconds);
    optimized_mps.set_overlap(std::abs(initial_mps.get_vector().dot(optimized_mps.get_vector())));

    if constexpr(settings::debug) {
        auto t_dbg = tid::tic_scope("debug");
        // Check that Ceres results are correct
        double energy_check   = tools::finite::measure::energy(optimized_mps.get_tensor(), tensors);
        double variance_check = tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors);
        if(std::abs(1.0 - std::abs(optimized_mps.get_energy() / energy_check)) > 1e-3)
            tools::log->warn("Energy mismatch: Ceres: {:.16f} | DMRG {:.16f}", optimized_mps.get_energy(), energy_check);
        if(std::abs(1.0 - std::abs(optimized_mps.get_variance() / variance_check)) > 1e-3)
            tools::log->warn("Variance mismatch: Ceres: {:8.2e} | DMRG {:8.2e}", optimized_mps.get_variance(), variance_check);
    }

    int    hrs = static_cast<int>(summary.total_time_in_seconds / 3600);
    int    min = static_cast<int>(std::fmod(summary.total_time_in_seconds, 3600) / 60);
    double sec = std::fmod(std::fmod(summary.total_time_in_seconds, 3600), 60);
    tools::log->debug("Finished BFGS in {:0<2}:{:0<2}:{:0<.1f} seconds and {} iters. Exit status: {}. Message: {}", hrs, min, sec, summary.iterations.size(),
                      ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
    reports::bfgs_add_entry("Subspace", "opt", optimized_mps, subspace_size);

    // Return the best theta
    return optimized_mps;
}

//
// opt_mps tools::finite::opt::internal::eigs_optimize_subspace(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status,
//                                                             OptMeta &meta) {
//    tools::log->trace("Optimizing subspace");
//    auto t_sub = tid::tic_scope("subspace");
//    initial_mps.validate_basis_vector();
//    if(meta.optMode != OptMode::SUBSPACE)
//        throw except::runtime_error("Wrong optimization mode [{}]. Expected [{}]", enum2sv(meta.optMode), enum2sv(OptMode::SUBSPACE));
//
//    /*
//     * Subspace optimization
//     *
//     *
//     * In subspace optimization we consider a local set of sites l of the L-site system,
//     * usually this corresponds to l = 1 or up to l = 8 adjacent sites.
//     *
//     * A subspace in this context is a truncated basis, i.e. a small subset of eigenvectors
//     * of the local "effective" Hamiltonian. The subspace is a set of k eigenvectors [x]
//     * which have significant overlap with the current state |y⟩, i.e. we define k such that
//     *
//     *          ε > 1 - Σ_i^k |⟨x_i|y⟩|²,
//     *
//     * where ε is a small number that controls error of the truncation. A value ε ~1e-10 is
//     * reasonable. Note that the truncation implies that k is smaller than the local
//     * Hilbert space dimension.
//     *
//     * After having found a good subspace, the objective is to find a linear combination
//     * of eigenvectors which minimizes the energy variance.
//     *
//     * It is worth noting some observations. Let {x} be the set of all eigenvectors
//     * to the local effective Hamiltonian H_local.
//     * Then, when the DMRG process is fully converged:
//     *      - only one x_i has overlap ⟨x_i|y⟩ = 1
//     *        Since the sum of all overlaps must add to 1, the rest have <x_j|y> = 0 when i != j.
//     *      - This x is also the one that minimizes the energy variance.
//     *
//     * However, before the DMRG process has converged this is not true. Instead:
//     *      - we have ⟨x_i|y⟩ > 0 for several i.
//     *      - a linear combination of several x can have lower variance than any
//     *        single x.
//     *
//     * Fully diagonalizing H_local yields all K eigenvectors {x}, but if H_local is too big this operation
//     * becomes prohibitively expensive. Instead we resort to finding a subset with k << K eigenvectors [x],
//     * whose eigenvalues are the k energies closest to the current energy. Usually the eigenvectors
//     * which have some overlap ⟨x_i|y⟩ > 0 are found in the subset [x] if k is large enough.
//     *
//     * Subspace optimization steps
//     *
//     * Step 1)  Find a subspace [x], i.e. take a set of k eigenvectors of the local effective Hamiltonian.
//     *          Empirically, eigenvectors whose eigenvalues (energy) are closest to the current energy,
//     *          tend to have nonzero overlap with the current vector |y>.
//     *          On an iterative solver we keep increasing "nev" (number of requested eigenvectors) until
//     *          the subspace error ε is small enough.
//     *          If any eigenvectors have to removed, (e.g. due to memory/performance constraints),
//     *          then sort the eigenvectors in order of decreasing overlap ⟨x_i|y⟩, and start deleting
//     *          from the end.
//     *
//     * Step 2)  Project the squared effective K*K Hamiltonian, down to the k*k subspace, H².
//     *          Using BFGS, find the linear combination |w⟩ of eigenvectors that minimizes the variance.
//     *
//     *              min_w Var H = ⟨H²⟩ - ⟨H⟩² = ⟨w|H²|w⟩ - ⟨E⟩²
//     *
//     *          where E are the energy eigenvalues from step 1.
//     *
//     */
//
//
//    // Handy references
//    const auto &state = *tensors.state;
//    const auto &model = *tensors.model;
//    const auto &edges = *tensors.edges;
//
//    /*
//     *  Step 1) Find the subspace.
//     *  The subspace is a set of eigenstates obtained from full or partial diagonalization
//     */
//
//    std::vector<opt_mps> subspace;
//    switch(meta.optType) {
//        case OptType::CPLX: subspace = internal::subspace::find_subspace<cplx>(tensors, settings::precision::target_subspace_error, meta); break;
//        case OptType::REAL: subspace = internal::subspace::find_subspace<double>(tensors, settings::precision::target_subspace_error, meta); break;
//    }
//
//    tools::log->trace("Subspace found with {} eigenvectors", subspace.size());
//
//    /*
//     * Filter the eigenvectors
//     *
//     */
//
//    internal::subspace::filter_subspace(subspace, settings::precision::max_subspace_size);
//
//    /*
//     *
//     * Step 2) Project the squared hamiltonian to the subspace
//     *
//     */
//
//    // Construct H² as a matrix (expensive operation!)
//    // Also make sure you do this before prepending the current state. All eigenvectors here should be basis vectors
//    Eigen::MatrixXcd H2_subspace = subspace::get_hamiltonian_squared_in_subspace<cplx>(model, edges, subspace);
//    if(meta.optType == OptType::REAL) H2_subspace = H2_subspace.real();
//
//    // Find the best eigenvectors and compute their variance
//    auto eigvecs_top_idx = internal::subspace::get_idx_to_eigvec_with_highest_overlap(subspace, 5, status.energy_llim, status.energy_ulim);
//    for(auto &idx : eigvecs_top_idx) {
//        auto &eigvec = *std::next(subspace.begin(), static_cast<long>(idx));
//        eigvec.set_variance(tools::finite::measure::energy_variance(eigvec.get_tensor(), tensors));
//    }
//
//    // We need the eigenvalues in a convenient format as well
//    auto eigvals = internal::subspace::get_eigvals(subspace);
//
//    // All the eigenvectors with significant overlap should be considered as initial guesses, including the current theta
//    // So we append it here to the list of eigenvectors
//    subspace.emplace_back(initial_mps);
//    eigvecs_top_idx.emplace_back(subspace.size() - 1);
//
//    tools::log->debug("Optimizing with {} initial guesses", eigvecs_top_idx.size());
//    std::vector<opt_mps> optimized_results;
//    size_t               eigvec_count = 0;
//    for(auto &idx : eigvecs_top_idx) {
//        const auto &eigvec = *std::next(subspace.begin(), static_cast<long>(idx));
//        tools::log->trace("Starting BFGS with eigvec {:<2} as initial guess: overlap {:.16f} | energy {:>20.16f} | variance: {:>8.2e} | eigvec {}", idx,
//                          eigvec.get_overlap(), eigvec.get_energy_per_site(), eigvec.get_variance(), eigvec.is_basis_vector);
//        Eigen::VectorXcd        subspace_vector = internal::subspace::get_vector_in_subspace(subspace, idx);
//        long                    subspace_size   = subspace_vector.size();
//        [[maybe_unused]] double norm;
//        reports::bfgs_add_entry("Subspace", "init", eigvec, subspace_size);
//        if constexpr(settings::debug) {
//            if(tools::log->level() == spdlog::level::trace) {
//                // Check using explicit matrix
//                auto             t_dbg               = tid::tic_scope("debug");
//                Eigen::VectorXcd initial_vector_sbsp = internal::subspace::get_vector_in_subspace(subspace, initial_mps.get_vector());
//                real             overlap_sbsp        = std::abs(subspace_vector.dot(initial_vector_sbsp));
//                Eigen::VectorXcd Hv                  = eigvals.asDiagonal() * subspace_vector;
//                Eigen::VectorXcd H2v                 = H2_subspace.template selfadjointView<Eigen::Upper>() * subspace_vector;
//                cplx             vHv                 = subspace_vector.dot(Hv);
//                cplx             vH2v                = subspace_vector.dot(H2v);
//                real             vv                  = subspace_vector.squaredNorm();
//                cplx             ene                 = vHv / vv;
//                cplx             var                 = vH2v / vv - ene * ene;
//                real             ene_init_san        = std::real(ene + model.get_energy_shift()) / static_cast<double>(state.get_length());
//                real             var_init_san        = std::real(var) / static_cast<double>(state.get_length());
//                std::string      description         = fmt::format("{:<8} {:<16} {}", "Subspace", eigvec.get_name(), "matrix check");
//                reports::bfgs_add_entry(description, subspace_vector.size(), subspace_size, ene_init_san, var_init_san, overlap_sbsp,
//                                        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), subspace_vector.norm(), 1, 1,
//                                        t_dbg->get_last_interval());
//            }
//            if(tools::log->level() == spdlog::level::trace) {
//                // Check current tensor
//                auto             t_dbg          = tid::tic_scope("debug");
//                Eigen::VectorXcd theta_0        = internal::subspace::get_vector_in_fullspace(subspace, subspace_vector);
//                auto             theta_0_tensor = tenx::TensorMap(theta_0, state.active_dimensions());
//                real             energy_0       = tools::finite::measure::energy_per_site(theta_0_tensor, tensors);
//                real             variance_0     = tools::finite::measure::energy_variance(theta_0_tensor, tensors);
//                real             overlap_0      = std::abs(initial_mps.get_vector().dot(theta_0));
//                std::string      description    = fmt::format("{:<8} {:<16} {}", "Subspace", eigvec.get_name(), "fullspace check");
//                reports::bfgs_add_entry(description, theta_0.size(), subspace_size, energy_0, variance_0, overlap_0, theta_0.norm(),
//                                        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 1, 1, t_dbg->get_last_interval());
//            }
//        }
//
//        /*
//         *
//         *  Start the BFGS optimization process for the subspace
//         *
//         */
//
//        auto options = internal::bfgs_default_options;
//        auto summary = ceres::GradientProblemSolver::Summary();
//        optimized_results.emplace_back(opt_mps());
//        auto &optimized_mps = optimized_results.back();
//        optimized_mps.set_name(eigvec.get_name());
//        optimized_mps.set_sites(eigvec.get_sites());
//        optimized_mps.set_length(eigvec.get_length());
//        optimized_mps.set_energy_shift(eigvec.get_energy_shift());
//        auto t_bfgs = tid::tic_scope("bfgs");
//        switch(meta.optType) {
//            case OptType::CPLX: {
//                Eigen::VectorXd subspace_vector_cast =
//                    Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(subspace_vector.data()), 2 * subspace_vector.size());
//                auto             *functor = new bfgs_subspace_functor<std::complex<double>>(tensors, status, H2_subspace, eigvals);
//                CustomLogCallback ceres_logger(*functor);
//                options.callbacks.emplace_back(&ceres_logger);
//                ceres::GradientProblem problem(functor);
//                tools::log->trace("Running BFGS subspace cplx");
//                ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
//                subspace_vector =
//                    Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<cplx *>(subspace_vector_cast.data()), subspace_vector_cast.size() / 2).normalized();
//                // Copy the results from the functor
//                optimized_mps.set_tensor(subspace::get_vector_in_fullspace(subspace, subspace_vector), initial_mps.get_tensor().dimensions());
//                optimized_mps.set_mv(functor->get_count());
//                optimized_mps.set_delta_f(functor->get_delta_f());
//                optimized_mps.set_grad_max(functor->get_max_grad_norm());
//                tid::get("vH2") += *functor->t_H2n;
//                tid::get("vH2v") += *functor->t_nH2n;
//                tid::get("vH") += *functor->t_Hn;
//                tid::get("vHv") += *functor->t_nHn;
//                tid::get("step") += *functor->t_step;
//                break;
//            }
//            case OptType::REAL: {
//                Eigen::VectorXd   subspace_vector_cast = subspace_vector.real();
//                Eigen::MatrixXd   H2_subspace_real     = H2_subspace.real();
//                auto             *functor              = new bfgs_subspace_functor<double>(tensors, status, H2_subspace_real, eigvals);
//                CustomLogCallback ceres_logger(*functor);
//                options.callbacks.emplace_back(&ceres_logger);
//                ceres::GradientProblem problem(functor);
//                tools::log->trace("Running BFGS subspace real");
//                ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
//                subspace_vector = subspace_vector_cast.normalized().cast<cplx>();
//                optimized_mps.set_tensor(subspace::get_vector_in_fullspace(subspace, subspace_vector), initial_mps.get_tensor().dimensions());
//                optimized_mps.set_mv(functor->get_count());
//                optimized_mps.set_delta_f(functor->get_delta_f());
//                optimized_mps.set_grad_max(functor->get_max_grad_norm());
//                tid::get("vH2") += *functor->t_H2n;
//                tid::get("vH2v") += *functor->t_nH2n;
//                tid::get("vH") += *functor->t_Hn;
//                tid::get("vHv") += *functor->t_nHn;
//                tid::get("step") += *functor->t_step;
//                break;
//            }
//        }
//        reports::time_add_entry();
//        t_bfgs.toc();
//        // Copy and set the rest of the tensor metadata
//        optimized_mps.normalize();
//        optimized_mps.set_energy(tools::finite::measure::energy(optimized_mps.get_tensor(), tensors));
//        optimized_mps.set_variance(tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors));
//        optimized_mps.set_overlap(std::abs(initial_mps.get_vector().dot(optimized_mps.get_vector())));
//        optimized_mps.set_iter(summary.iterations.size());
//        optimized_mps.set_time(summary.total_time_in_seconds);
//        optimized_mps.set_overlap(std::abs(initial_mps.get_vector().dot(optimized_mps.get_vector())));
//
//        if constexpr(settings::debug) {
//            auto t_dbg = tid::tic_scope("debug");
//            // Check that Ceres results are correct
//            double energy_check   = tools::finite::measure::energy(optimized_mps.get_tensor(), tensors);
//            double variance_check = tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors);
//            if(std::abs(1.0 - std::abs(optimized_mps.get_energy() / energy_check)) > 1e-3)
//                tools::log->warn("Energy mismatch: Ceres: {:.16f} | DMRG {:.16f}", optimized_mps.get_energy(), energy_check);
//            if(std::abs(1.0 - std::abs(optimized_mps.get_variance() / variance_check)) > 1e-3)
//                tools::log->warn("Variance mismatch: Ceres: {:8.2e} | DMRG {:8.2e}", optimized_mps.get_variance(), variance_check);
//        }
//
//        int    hrs = static_cast<int>(summary.total_time_in_seconds / 3600);
//        int    min = static_cast<int>(std::fmod(summary.total_time_in_seconds, 3600) / 60);
//        double sec = std::fmod(std::fmod(summary.total_time_in_seconds, 3600), 60);
//        tools::log->debug("Finished BFGS in {:0<2}:{:0<2}:{:0<.1f} seconds and {} iters. Exit status: {}. Message: {}", hrs, min, sec,
//                          summary.iterations.size(), ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
//        //    std::cout << summary.FullReport() << "\n";
//        if(tools::log->level() <= spdlog::level::debug) { reports::bfgs_add_entry("Subspace", "opt", optimized_mps, subspace_size); }
//        eigvec_count++;
//    }
//
//    // Sort thetas in ascending order in variance
//    std::sort(optimized_results.begin(), optimized_results.end(),
//              [](const opt_mps &lhs, const opt_mps &rhs) { return lhs.get_variance() < rhs.get_variance(); });
//
//    // Return the best theta
//    return optimized_results.front();
//}
