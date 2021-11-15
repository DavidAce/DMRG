#include <math/tenx.h>
// -- (textra first)
#include "../opt_meta.h"
#include "../opt_mps.h"
#include "ceres_subspace_functor.h"
#include "opt-internal.h"
#include "report.h"
#include <algorithms/AlgorithmStatus.h>
#include <ceres/gradient_problem.h>
#include <config/settings.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;

template<typename Scalar>
std::vector<opt_mps> internal::subspace::find_subspace(const TensorsFinite &tensors, double target_subspace_error, const OptMeta &meta) {
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    tools::log->trace("Finding subspace");
    auto t_find     = tid::tic_scope("find");
    auto dbl_length = static_cast<double>(state.get_length());

    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;

    // If the mps is small enough you can afford full diag.
    if(tensors.state->active_problem_size() <= settings::precision::max_size_full_diag) {
        std::tie(eigvecs, eigvals) = find_subspace_full<Scalar>(tensors);
    } else {
        double eigval_target;
        double energy_target = tools::finite::measure::energy(tensors);
        if(model.is_reduced()) {
            eigval_target = tools::finite::measure::energy_minus_energy_reduced(tensors);
            tools::log->trace("Energy reduce = {:.16f} | per site = {:.16f}", model.get_energy_reduced(), model.get_energy_per_site_reduced());
            tools::log->trace("Energy target = {:.16f} | per site = {:.16f}", energy_target, energy_target / dbl_length);
            tools::log->trace("Eigval target = {:.16f} | per site = {:.16f}", eigval_target, eigval_target / dbl_length);
            tools::log->trace("Eigval target + Energy reduce = Energy: {:.16f} + {:.16f} = {:.16f}", eigval_target / dbl_length,
                              model.get_energy_per_site_reduced(), energy_target / dbl_length);
        } else {
            eigval_target = energy_target;
        }
        std::tie(eigvecs, eigvals) = find_subspace_part<Scalar>(tensors, eigval_target, target_subspace_error, meta);
    }
    tools::log->trace("Eigval range         : {:.16f} --> {:.16f}", eigvals.minCoeff(), eigvals.maxCoeff());
    tools::log->trace("Energy range         : {:.16f} --> {:.16f}", eigvals.minCoeff() + model.get_energy_reduced(),
                      eigvals.maxCoeff() + model.get_energy_reduced());
    tools::log->trace("Energy range per site: {:.16f} --> {:.16f}", eigvals.minCoeff() / dbl_length + model.get_energy_per_site_reduced(),
                      eigvals.maxCoeff() / dbl_length + model.get_energy_per_site_reduced());
    reports::print_eigs_report();

    if constexpr(std::is_same<Scalar, double>::value) {
        tenx::subtract_phase(eigvecs);
        double trunc = eigvecs.imag().cwiseAbs().sum();
        if(trunc > 1e-12) tools::log->warn("truncating imag of eigvecs, sum: {}", trunc);
        eigvecs = eigvecs.real();
    }

    const auto     &multisite_mps  = state.get_multisite_mps();
    const auto      multisite_vec  = Eigen::Map<const Eigen::VectorXcd>(multisite_mps.data(), multisite_mps.size());
    auto            energy_reduced = model.get_energy_reduced();
    Eigen::VectorXd overlaps       = (multisite_vec.adjoint() * eigvecs).cwiseAbs().real();

    double eigvec_time = 0;
    for(const auto &item : reports::eigs_log) { eigvec_time += item.ham_time + item.lu_time + item.eig_time; }

    std::vector<opt_mps> subspace;
    subspace.reserve(static_cast<size_t>(eigvals.size()));
    for(long idx = 0; idx < eigvals.size(); idx++) {
        // Important to normalize the eigenvectors that we get from the solver: they are not always well normalized when we get them!
        auto eigvec_i = tenx::TensorCast(eigvecs.col(idx).normalized(), state.active_dimensions());
        subspace.emplace_back(fmt::format("eigenvector {}", idx), eigvec_i, tensors.active_sites, eigvals(idx), energy_reduced, std::nullopt, overlaps(idx),
                              tensors.get_length());
        subspace.back().set_time(eigvec_time);
        subspace.back().set_mv(reports::eigs_log.size());
        subspace.back().set_iter(reports::eigs_log.size());
        subspace.back().is_basis_vector = true;

        subspace.back().set_krylov_idx(idx);
        subspace.back().set_krylov_eigval(eigvals(idx));
        subspace.back().set_krylov_ritz(enum2sv(meta.optRitz));
        subspace.back().set_optmode(meta.optMode);
        subspace.back().set_optspace(meta.optSpace);

        subspace.back().validate_basis_vector();
    }
    return subspace;
}

template std::vector<opt_mps> internal::subspace::find_subspace<cplx>(const TensorsFinite &tensors, double target_subspace_error, const OptMeta &meta);

template std::vector<opt_mps> internal::subspace::find_subspace<real>(const TensorsFinite &tensors, double target_subspace_error, const OptMeta &meta);

opt_mps tools::finite::opt::internal::ceres_optimize_subspace(const TensorsFinite &tensors, const opt_mps &initial_mps, const std::vector<opt_mps> &subspace,
                                                              const Eigen::MatrixXcd &H2_subspace, const AlgorithmStatus &status, OptMeta &meta) {
    auto t_sub = tid::tic_scope("subspace");
    initial_mps.validate_basis_vector();

    // Handy references
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;

    opt_mps optimized_mps;
    auto    eigvec_best_idx = internal::subspace::get_idx_to_eigvec_with_highest_overlap(subspace, status.energy_llim_per_site, status.energy_ulim_per_site);
    if(eigvec_best_idx)
        optimized_mps = *std::next(subspace.begin(), static_cast<long>(eigvec_best_idx.value()));
    else
        optimized_mps = *subspace.begin();
    auto             fullspace_dims  = optimized_mps.get_tensor().dimensions();
    Eigen::VectorXcd subspace_vector = internal::subspace::get_vector_in_subspace(subspace, optimized_mps.get_vector());

    tools::log->debug("Optimizing with initial guess: {}", optimized_mps.get_name());

    // We need the eigenvalues in a convenient format as well
    auto eigvals = internal::subspace::get_eigvals(subspace);

    /*
     *
     *  Sanity checks
     *
     */

    if(tools::log->level() <= spdlog::level::debug) {
        if constexpr(settings::debug) {
            if(tools::log->level() == spdlog::level::trace) {
                // Check using explicit matrix
                auto             t_dbg               = tid::tic_scope("debug");
                Eigen::VectorXcd initial_vector_sbsp = internal::subspace::get_vector_in_subspace(subspace, initial_mps.get_vector());
                real             overlap_sbsp        = std::abs(subspace_vector.dot(initial_vector_sbsp));
                Eigen::VectorXcd Hv                  = eigvals.asDiagonal() * subspace_vector;
                Eigen::VectorXcd H2v                 = H2_subspace.template selfadjointView<Eigen::Upper>() * subspace_vector;
                cplx             vHv                 = subspace_vector.dot(Hv);
                cplx             vH2v                = subspace_vector.dot(H2v);
                real             vv                  = subspace_vector.squaredNorm();
                cplx             ene                 = vHv / vv;
                cplx             var                 = vH2v / vv - ene * ene;
                real             ene_init_san        = std::real(ene + model.get_energy_reduced()) / static_cast<double>(state.get_length());
                real             var_init_san        = std::real(var) / static_cast<double>(state.get_length());
                std::string      description         = fmt::format("{:<8} {:<16} {}", "Subspace", optimized_mps.get_name(), "matrix check");
                reports::bfgs_add_entry(description, subspace_vector.size(), subspace_vector.size(), ene_init_san, var_init_san, overlap_sbsp,
                                        subspace_vector.norm(), std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 1, 1,
                                        t_dbg->get_last_interval());
            }
            if(tools::log->level() == spdlog::level::trace) {
                // Check current tensor
                auto             t_dbg          = tid::tic_scope("debug");
                Eigen::VectorXcd theta_0        = internal::subspace::get_vector_in_fullspace(subspace, subspace_vector);
                auto             theta_0_tensor = tenx::TensorMap(theta_0, state.active_dimensions());
                real             energy_0       = tools::finite::measure::energy_per_site(theta_0_tensor, tensors);
                real             variance_0     = tools::finite::measure::energy_variance(theta_0_tensor, tensors);
                real             overlap_0      = std::abs(initial_mps.get_vector().dot(theta_0));
                std::string      description    = fmt::format("{:<8} {:<16} {}", "Subspace", initial_mps.get_name(), "fullspace check");
                reports::bfgs_add_entry(description, theta_0.size(), subspace_vector.size(), energy_0, variance_0, overlap_0, theta_0.norm(),
                                        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 1, 1, t_dbg->get_last_interval());
            }
        }
    }

    /*
     *
     *  Start the LBFGS optimization process for the subspace
     *
     */

    auto options = internal::ceres_default_options;
    auto summary = ceres::GradientProblemSolver::Summary();
    auto t_lbfgs = tid::tic_scope("lbfgs");
    switch(meta.optType) {
        case OptType::CPLX: {
            Eigen::VectorXd subspace_vector_cast = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(subspace_vector.data()), 2 * subspace_vector.size());
            auto           *functor              = new ceres_subspace_functor<std::complex<double>>(tensors, status, H2_subspace, eigvals);
            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS subspace cplx");
            ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
            subspace_vector = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<cplx *>(subspace_vector_cast.data()), subspace_vector_cast.size() / 2).normalized();
            // Copy the results from the functor
            optimized_mps.set_tensor(subspace::get_tensor_in_fullspace(subspace, subspace_vector, fullspace_dims));
            optimized_mps.set_mv(functor->get_count());
            tid::get("vH2") += *functor->t_H2n;
            tid::get("vH2v") += *functor->t_nH2n;
            tid::get("vH") += *functor->t_Hn;
            tid::get("vHv") += *functor->t_nHn;
            tid::get("step") += *functor->t_step;
            break;
        }
        case OptType::REAL: {
            Eigen::VectorXd   subspace_vector_cast = subspace_vector.real();
            Eigen::MatrixXd   H2_subspace_real     = H2_subspace.real();
            auto             *functor              = new ceres_subspace_functor<double>(tensors, status, H2_subspace_real, eigvals);
            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS subspace real");
            ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
            subspace_vector = subspace_vector_cast.normalized().cast<cplx>();
            optimized_mps.set_tensor(subspace::get_tensor_in_fullspace(subspace, subspace_vector, fullspace_dims));
            optimized_mps.set_mv(functor->get_count());
            tid::get("vH2") += *functor->t_H2n;
            tid::get("vH2v") += *functor->t_nH2n;
            tid::get("vH") += *functor->t_Hn;
            tid::get("vHv") += *functor->t_nHn;
            tid::get("step") += *functor->t_step;
            break;
        }
    }
    reports::time_add_opt_entry();
    t_lbfgs.toc();
    // Copy and set the rest of the tensor metadata
    optimized_mps.normalize();
    optimized_mps.set_energy(tools::finite::measure::energy(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_variance(tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_overlap(std::abs(initial_mps.get_vector().dot(optimized_mps.get_vector())));
    optimized_mps.set_iter(summary.iterations.size());
    optimized_mps.set_time(summary.total_time_in_seconds);

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
    tools::log->debug("Finished LBFGS in {:0<2}:{:0<2}:{:0<.1f} seconds and {} iters. Exit status: {}. Message: {}", hrs, min, sec, summary.iterations.size(),
                      ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
    //    std::cout << summary.FullReport() << "\n";
    if(tools::log->level() <= spdlog::level::debug) { reports::bfgs_add_entry("Subspace", "opt", optimized_mps, subspace_vector.size()); }

    //    if(optSpace == OptSpace::SUBSPACE) {
    //        tools::log->trace("SUBSPACE optimization done. Starting fine tuning with DIRECT optimization");
    //        optimized_mps = ceres_direct_optimization(tensors, optimized_mps, status, optType, optMode, optSpace);
    //    }
    // Return the optimized result
    return optimized_mps;
}

opt_mps tools::finite::opt::internal::ceres_subspace_optimization(const TensorsFinite &tensors, const AlgorithmStatus &status, OptMeta &meta) {
    std::vector<size_t> sites(tensors.active_sites.begin(), tensors.active_sites.end());
    opt_mps             initial_mps("current state", tensors.get_multisite_mps(), sites,
                                    tools::finite::measure::energy(tensors) - tensors.model->get_energy_reduced(), // Eigval
                                    tensors.model->get_energy_reduced(),                                           // Energy reduced for full system
                                    tools::finite::measure::energy_variance(tensors),
                                    1.0, // Overlap
                                    tensors.get_length());
    return ceres_subspace_optimization(tensors, initial_mps, status, meta);
}

opt_mps tools::finite::opt::internal::ceres_subspace_optimization(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status,
                                                                  OptMeta &meta) {
    tools::log->trace("Optimizing in SUBSPACE mode");
    auto t_sub = tid::tic_scope("subspace");
    initial_mps.validate_basis_vector();

    /*
     * Subspace optimization
     *
     * Introduction
     * In subspace optimization we consider a "local" subsection of the "global" L-site system,
     * usually this corresponds to n=2 local sites but in multisite dmrg we may consider
     * any n=[2,8] adjacent sites.
     *
     * The point here is to minimize the global energy variance Var H_global, by tuning only the parameters of
     * the local wavefunction |psi>_local corresponding to the local Hamiltonian H_local for n sites.
     * In other words, we seek the local vector which minimizes the global energy variance.
     *
     * It is worth noting some properties which hold when the DMRG process is fully converged:
     *      - If {x} are all the eigenvectors of H_local, only one of them has overlap <x_i|y> = 1
     *        where y is the current local vector for n sites.
     *        Since the sum of all overlaps must add to 1, the rest have <x_j|y> = 0 when i != j.
     *      - This particular eigenvector is also the one that minimizes the Var H_global,
     *        and no further optimization should be needed.
     *
     * However, before the DMRG process has converged this is not true. Instead:
     *      - If {x} are all the eigenvectors of H_local we have <x_i|y> > 0 for several i.
     *      - None {x} minimizes Var H_global, but a linear combination of several x could.
     *
     * Fully diagonalizing H_local yields all K eigenvectors {x}, but if H_local is too big this operation
     * becomes prohibitively expensive. Instead we resort to finding a subset with k << K eigenvectors [x],
     * whose eigenvalues are the k energies closest to the current energy. Usually the eigenvectors
     * which have some overlap <x_i|y> > 0 are found in the subset [x] if k is large enough.
     *
     * In Subspace Optimization we find a linear combination eigenvectors in the subspace [x] which
     * minimizes Var H_global, hence the name.
     *
     *
     *
     * Subspace optimization steps
     *
     * Step 0)  Find a subspace, i.e. a set of k eigenvectors [x] to the local Hamiltonian H_local with
     *          energy eigenvalue closest to the current energy. Note that H_local is a K * K matrix,
     *          and k << K. The set [x] is sorted in order of descending overlap <x_i|y>,
     *          where y is the current vector.
     *
     *
     * Step 1)
     *      - (O) If  OptMode::OVERLAP:
     *            Find the index for the best eigvec inside of the energy window:
     *              - OA) idx >= 0: eigvec found in window. Return that eigvec
     *              - OB) idx = -1: No eigvec found in window. Return old tensor
     *       -(V) If OptMode::SUBSPACE
     *          Make sure k is as small as possible. I.e. filter out eigenvectors from [x] down
     *          to an even smaller set of "relevant" eigvecs for doing subspace optimization.
     *          Allowing a maximum of k == 64 eigvecs keeps ram below 2GB when the problem
     *          size == 4096 (the linear size of H_local and the mps multisite_tensor).
     *          This means that we filter out
     *              * eigvecs outside of the energy window (if one is enabled)
     *              * eigvecs with little or no overlap to the current state.
     *          The filtering process collects the eigvec state with best overlap until
     *          either of these two criteria is met:
     *              * More than N=max_accept eigvecs have been collected
     *              * The subspace error "eps" is low enough¹
     *          The filtering returns the eigvecs sorted in decreasing overlap.
     *          ¹ We define eps = 1 - Σ_i |<x_i|y>|². A value eps ~1e-10 is reasonable.
     *
     * Step 2)
     *          Find the index for the best eigvec inside of the energy window:
     *              - VA) idx = -1: No eigvec found in window. Return old tensor
     *              - VB) We have many eigvecs, including the current tensor.
     *                    They should be viewed as good starting guesses for LBFGS optimization.
     *                    Optimize them one by one, and keep the
     *
     * Step 2)  Find the best overlapping state among the relevant eigvecs.
     * Step 3)  We can now make different decisions based on the overlap.
     *          A)  If best_overlap_idx == -1
     *              No state is in energy window -> discard! Return old multisite_mps_tensor.
     *          B)  If overlap_high <= best_overlap.
     *              This can happen if the environments have been modified just slightly since the last time considered
     *              these sites, but the signal is still clear -- we are still targeting the same state.
     *              However we can't be sure that the contributions from nearby states is just noise. Instead of just
     *              keeping the state we should optimize its variance. This is important in the later stages when variance
     *              is low and we don't want to ruin those last decimals.
     *              We just need to decide which initial guess to use.
     *                  B1) If best_overlap_variance <= theta_variance: set theta_initial = best_overlap_theta.
     *                  B2) Else, set theta_initial = multisite_mps_tensor.
     *          C)  If overlap_cat <= best_overlap and best_overlap < overlap_high
     *              This can happen for one reasons:
     *                  1) There are a few eigvec states with significant overlap (superposition)
     *              It's clear that we need to optimize, but we have to think carefully about the initial guess.
     *              Right now it makes sense to always choose best overlap theta, since that forces the algorithm to
     *              choose a particular state and not get stuck in superposition. Choosing the old theta may just entrench
     *              the algorithm into a local minima.
     *          D)  If 0 <= best_overlap and best_overlap < overlap_cat
     *              This can happen for three reasons, most often early in the simulation.
     *                  1) There are several eigvec states with significant overlap (superposition)
     *                  2) The highest overlapping states were outside of the energy window, leaving just these eigvecs.
     *                  3) The energy targeting of states has failed for some reason, perhaps the spectrum is particularly dense_lu.
     *              In any case, it is clear we are lost Hilbert space.
     *              Also, the subspace_error is no longer a good measure of how useful the subspace is to us, since it's only
     *              measuring how well the old state can be described, but the old state is likely very different from what
     *              we're looking for.
     *              So to address all three cases, do DIRECT optimization with best_overlap_theta as initial guess.
     *
     * In particular, notice that we never use the eigvec that happens to have the best variance.
     */
    // Handy references
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    const auto &edges = *tensors.edges;

    /*
     *  Step 0) Find the subspace.
     *  The subspace set of eigenstates obtained from full or partial diagonalization
     */

    std::vector<opt_mps> subspace;
    switch(meta.optType) {
        case OptType::CPLX: subspace = internal::subspace::find_subspace<cplx>(tensors, settings::precision::target_subspace_error, meta); break;
        case OptType::REAL: subspace = internal::subspace::find_subspace<double>(tensors, settings::precision::target_subspace_error, meta); break;
    }

    tools::log->trace("Subspace found with {} eigenvectors", subspace.size());

    /*
     *
     * Step 1) (Overlap mode) Return best overlap
     *
     */
    if(meta.optMode == OptMode::OVERLAP) {
        auto max_overlap_idx = internal::subspace::get_idx_to_eigvec_with_highest_overlap(subspace, status.energy_llim_per_site, status.energy_ulim_per_site);
        if(max_overlap_idx) {
            // (OA)
            auto &eigvec_max_overlap = *std::next(subspace.begin(), static_cast<long>(max_overlap_idx.value()));
            eigvec_max_overlap.set_variance(tools::finite::measure::energy_variance(eigvec_max_overlap.get_tensor(), tensors));
            if(tools::log->level() == spdlog::level::trace) {
                tools::log->trace("ceres_subspace_optimization: eigvec {:<2} has highest overlap {:.16f} | energy {:>20.16f} | variance {:>8.2e}",
                                  max_overlap_idx.value(), eigvec_max_overlap.get_overlap(), eigvec_max_overlap.get_energy_per_site(),
                                  eigvec_max_overlap.get_variance());
            }
            if(eigvec_max_overlap.get_overlap() < 0.1)
                tools::log->debug("ceres_subspace_optimization: Overlap fell below < 0.1: {:20.16f}", eigvec_max_overlap.get_overlap());
            return eigvec_max_overlap;
        } else {
            // (OB)
            tools::log->warn("ceres_subspace_optimization: No overlapping states in energy range. Returning old tensor");
            return initial_mps;
        }
    }

    /*
     *  Step 1)  (Variance mode) Filter the eigenvectors
     *
     */

    internal::subspace::filter_subspace(subspace, settings::precision::max_subspace_size);

    /*
     *
     * Step 2) (Variance mode)
     *
     */

    // Construct H² as a matrix (expensive operation!)
    // Also make sure you do this before prepending the current state. All eigenvectors here should be basis vectors
    Eigen::MatrixXcd H2_subspace = internal::get_multisite_hamiltonian_squared_subspace_matrix<cplx>(model, edges, subspace);
    if(meta.optType == OptType::REAL) H2_subspace = H2_subspace.real();

    // Find the best eigenvectors and compute their variance
    auto eigvecs_top_idx = internal::subspace::get_idx_to_eigvec_with_highest_overlap(subspace, 5, status.energy_llim_per_site, status.energy_ulim_per_site);
    for(auto &idx : eigvecs_top_idx) {
        auto &eigvec = *std::next(subspace.begin(), static_cast<long>(idx));
        eigvec.set_variance(tools::finite::measure::energy_variance(eigvec.get_tensor(), tensors));
    }

    // We need the eigenvalues in a convenient format as well
    auto eigvals = internal::subspace::get_eigvals(subspace);

    // All the eigenvectors with significant overlap should be considered as initial guesses, including the current theta
    // So we append it here to the list of eigenvectors
    subspace.emplace_back(initial_mps);
    eigvecs_top_idx.emplace_back(subspace.size() - 1);

    tools::log->debug("Optimizing with {} initial guesses", eigvecs_top_idx.size());
    std::vector<opt_mps> optimized_results;
    size_t               eigvec_count = 0;
    for(auto &idx : eigvecs_top_idx) {
        const auto &eigvec = *std::next(subspace.begin(), static_cast<long>(idx));
        tools::log->trace("Starting LBFGS with eigvec {:<2} as initial guess: overlap {:.16f} | energy {:>20.16f} | variance: {:>8.2e} | eigvec {}", idx,
                          eigvec.get_overlap(), eigvec.get_energy_per_site(), eigvec.get_variance(), eigvec.is_basis_vector);
        Eigen::VectorXcd        subspace_vector = internal::subspace::get_vector_in_subspace(subspace, idx);
        long                    subspace_size   = subspace_vector.size();
        [[maybe_unused]] double norm;
        reports::bfgs_add_entry("Subspace", "init", eigvec, subspace_size);
        if constexpr(settings::debug) {
            if(tools::log->level() == spdlog::level::trace) {
                // Check using explicit matrix
                auto             t_dbg               = tid::tic_scope("debug");
                Eigen::VectorXcd initial_vector_sbsp = internal::subspace::get_vector_in_subspace(subspace, initial_mps.get_vector());
                real             overlap_sbsp        = std::abs(subspace_vector.dot(initial_vector_sbsp));
                Eigen::VectorXcd Hv                  = eigvals.asDiagonal() * subspace_vector;
                Eigen::VectorXcd H2v                 = H2_subspace.template selfadjointView<Eigen::Upper>() * subspace_vector;
                cplx             vHv                 = subspace_vector.dot(Hv);
                cplx             vH2v                = subspace_vector.dot(H2v);
                real             vv                  = subspace_vector.squaredNorm();
                cplx             ene                 = vHv / vv;
                cplx             var                 = vH2v / vv - ene * ene;
                real             ene_init_san        = std::real(ene + model.get_energy_reduced()) / static_cast<double>(state.get_length());
                real             var_init_san        = std::real(var) / static_cast<double>(state.get_length());
                std::string      description         = fmt::format("{:<8} {:<16} {}", "Subspace", eigvec.get_name(), "matrix check");
                reports::bfgs_add_entry(description, subspace_vector.size(), subspace_size, ene_init_san, var_init_san, overlap_sbsp,
                                        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), subspace_vector.norm(), 1, 1,
                                        t_dbg->get_last_interval());
            }
            if(tools::log->level() == spdlog::level::trace) {
                // Check current tensor
                auto             t_dbg          = tid::tic_scope("debug");
                Eigen::VectorXcd theta_0        = internal::subspace::get_vector_in_fullspace(subspace, subspace_vector);
                auto             theta_0_tensor = tenx::TensorMap(theta_0, state.active_dimensions());
                real             energy_0       = tools::finite::measure::energy_per_site(theta_0_tensor, tensors);
                real             variance_0     = tools::finite::measure::energy_variance(theta_0_tensor, tensors);
                real             overlap_0      = std::abs(initial_mps.get_vector().dot(theta_0));
                std::string      description    = fmt::format("{:<8} {:<16} {}", "Subspace", eigvec.get_name(), "fullspace check");
                reports::bfgs_add_entry(description, theta_0.size(), subspace_size, energy_0, variance_0, overlap_0, theta_0.norm(),
                                        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), 1, 1, t_dbg->get_last_interval());
            }
        }

        /*
         *
         *  Start the LBFGS optimization process for the subspace
         *
         */

        auto options = internal::ceres_default_options;
        auto summary = ceres::GradientProblemSolver::Summary();
        optimized_results.emplace_back(opt_mps());
        auto &optimized_mps = optimized_results.back();
        optimized_mps.set_name(eigvec.get_name());
        optimized_mps.set_sites(eigvec.get_sites());
        optimized_mps.set_length(eigvec.get_length());
        optimized_mps.set_energy_reduced(eigvec.get_energy_reduced());
        auto t_lbfgs = tid::tic_scope("lbfgs");
        switch(meta.optType) {
            case OptType::CPLX: {
                Eigen::VectorXd subspace_vector_cast =
                    Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(subspace_vector.data()), 2 * subspace_vector.size());
                auto             *functor = new ceres_subspace_functor<std::complex<double>>(tensors, status, H2_subspace, eigvals);
                CustomLogCallback ceres_logger(*functor);
                options.callbacks.emplace_back(&ceres_logger);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running LBFGS subspace cplx");
                ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
                subspace_vector =
                    Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<cplx *>(subspace_vector_cast.data()), subspace_vector_cast.size() / 2).normalized();
                // Copy the results from the functor
                optimized_mps.set_tensor(subspace::get_vector_in_fullspace(subspace, subspace_vector), initial_mps.get_tensor().dimensions());
                optimized_mps.set_mv(functor->get_count());
                optimized_mps.set_delta_f(functor->get_delta_f());
                optimized_mps.set_max_grad(functor->get_max_grad_norm());
                tid::get("vH2") += *functor->t_H2n;
                tid::get("vH2v") += *functor->t_nH2n;
                tid::get("vH") += *functor->t_Hn;
                tid::get("vHv") += *functor->t_nHn;
                tid::get("step") += *functor->t_step;
                break;
            }
            case OptType::REAL: {
                Eigen::VectorXd   subspace_vector_cast = subspace_vector.real();
                Eigen::MatrixXd   H2_subspace_real     = H2_subspace.real();
                auto             *functor              = new ceres_subspace_functor<double>(tensors, status, H2_subspace_real, eigvals);
                CustomLogCallback ceres_logger(*functor);
                options.callbacks.emplace_back(&ceres_logger);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running LBFGS subspace real");
                ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
                subspace_vector = subspace_vector_cast.normalized().cast<cplx>();
                optimized_mps.set_tensor(subspace::get_vector_in_fullspace(subspace, subspace_vector), initial_mps.get_tensor().dimensions());
                optimized_mps.set_mv(functor->get_count());
                optimized_mps.set_delta_f(functor->get_delta_f());
                optimized_mps.set_max_grad(functor->get_max_grad_norm());
                tid::get("vH2") += *functor->t_H2n;
                tid::get("vH2v") += *functor->t_nH2n;
                tid::get("vH") += *functor->t_Hn;
                tid::get("vHv") += *functor->t_nHn;
                tid::get("step") += *functor->t_step;
                break;
            }
        }
        reports::time_add_opt_entry();
        t_lbfgs.toc();
        // Copy and set the rest of the tensor metadata
        optimized_mps.normalize();
        optimized_mps.set_energy(tools::finite::measure::energy(optimized_mps.get_tensor(), tensors));
        optimized_mps.set_variance(tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors));
        optimized_mps.set_overlap(std::abs(initial_mps.get_vector().dot(optimized_mps.get_vector())));
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
        tools::log->debug("Finished LBFGS in {:0<2}:{:0<2}:{:0<.1f} seconds and {} iters. Exit status: {}. Message: {}", hrs, min, sec,
                          summary.iterations.size(), ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
        //    std::cout << summary.FullReport() << "\n";
        if(tools::log->level() <= spdlog::level::debug) { reports::bfgs_add_entry("Subspace", "opt", optimized_mps, subspace_size); }
        eigvec_count++;
    }

    // Sort thetas in ascending order in variance
    std::sort(optimized_results.begin(), optimized_results.end(),
              [](const opt_mps &lhs, const opt_mps &rhs) { return lhs.get_variance() < rhs.get_variance(); });

    // Return the best theta
    return optimized_results.front();
}
