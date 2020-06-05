//
// Created by david on 2019-07-15.
//
#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include "ceres_subspace_functor.h"
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <math/nmspc_random.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt_tensor.h>
using namespace tools::finite::opt;

template<typename Scalar>
std::vector<opt_tensor> internal::subspace::find_candidates(const class_tensors_finite &tensors, double subspace_error_threshold, OptMode optMode,
                                                            OptSpace optSpace) {
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    const auto &edges = *tensors.edges;
    tools::log->trace("Finding subspace");
    MatrixType<Scalar> H_local          = tools::finite::opt::internal::get_multisite_hamiltonian_matrix<Scalar>(model, edges);
    const auto &       multisite_tensor = state.get_multisite_tensor();

    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;

    // If multitheta is small enough you can afford full diag.
    if(multisite_tensor.size() <= settings::precision::max_size_full_diag) {
        std::tie(eigvecs, eigvals) = find_subspace_full(H_local, multisite_tensor);
    } else {
        double energy_target;
        if(model.is_reduced()) energy_target = tools::finite::measure::energy_minus_energy_reduced(multisite_tensor, tensors);
        else
            energy_target = tools::finite::measure::energy(multisite_tensor, tensors);
        tools::log->trace("Energy target + energy reduced = energy per site: {} + {} = {}", energy_target / static_cast<double>(state.get_length()),
                          model.get_energy_reduced() / static_cast<double>(state.get_length()),
                          (energy_target + model.get_energy_reduced()) / static_cast<double>(state.get_length()));

        std::tie(eigvecs, eigvals) = find_subspace_part(H_local, multisite_tensor, energy_target, subspace_error_threshold, optMode, optSpace);
    }
    tools::log->trace("Eigenvalue range: {} --> {}", (eigvals.minCoeff() + model.get_energy_reduced()) / static_cast<double>(state.get_length()),
                      (eigvals.maxCoeff() + model.get_energy_reduced()) / static_cast<double>(state.get_length()));
    reports::print_eigs_report();

    if constexpr(std::is_same<Scalar, double>::value) {
        Textra::subtract_phase(eigvecs);
        tools::log->trace("truncating imag of eigvecs, sum: {}", eigvecs.imag().cwiseAbs().sum());
        eigvecs = eigvecs.real();
    }

    const auto &    multisite_mps_tensor = state.get_multisite_tensor();
    const auto      multisite_mps_vector = Eigen::Map<const Eigen::VectorXcd>(multisite_mps_tensor.data(), multisite_mps_tensor.size());
    auto            energy_reduced       = model.get_energy_reduced() / static_cast<double>(state.get_length());
    Eigen::VectorXd overlaps             = (multisite_mps_vector.adjoint() * eigvecs).cwiseAbs().real();

    std::vector<opt_tensor> candidate_tensors;
    for(long idx = 0; idx < eigvals.size(); idx++) {
        auto tensor_map = Textra::MatrixTensorMap(eigvecs.col(idx), state.active_dimensions());
        candidate_tensors.emplace_back(fmt::format("eigenvector {}", idx), tensor_map, tensors.active_sites, eigvals(idx), energy_reduced, std::nullopt,
                                       overlaps(idx));
        candidate_tensors.back().is_basis_vector = true;
    }
    return candidate_tensors;
}

template std::vector<opt_tensor> internal::subspace::find_candidates<cplx>(const class_tensors_finite &tensors, double subspace_error_threshold,
                                                                           OptMode optMode, OptSpace optSpace);

template std::vector<opt_tensor> internal::subspace::find_candidates<real>(const class_tensors_finite &tensors, double subspace_error_threshold,
                                                                           OptMode optMode, OptSpace optSpace);

opt_tensor tools::finite::opt::internal::ceres_subspace_optimization(const class_tensors_finite &tensors, const class_algorithm_status &status, OptType optType,
                                                                     OptMode optMode, OptSpace optSpace) {
    std::vector<size_t> sites(tensors.active_sites.begin(), tensors.active_sites.end());
    opt_tensor          initial_tensor("Current state", tensors.state->get_multisite_tensor(), sites,
                              0.0, // Energy reduced
                              tensors.model->get_energy_reduced(), tools::finite::measure::energy_variance_per_site(tensors),
                              1.0 // Overlap
    );
    return ceres_subspace_optimization(tensors, initial_tensor, status, optType, optMode, optSpace);
}

opt_tensor tools::finite::opt::internal::ceres_subspace_optimization(const class_tensors_finite &tensors, const opt_tensor &initial_tensor,
                                                                     const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                     OptSpace optSpace) {
    tools::log->trace("Optimizing in SUBSPACE mode");
    tools::common::profile::t_opt_sub->tic();

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
     * Step 1)  Make sure k is as small as possible. I.e. filter out eigenvectors from [x] down
     *          to an even smaller set of "relevant" candidates for doing subspace optimization.
     *          Allowing a maximum of k == 64 candidates keeps ram below 2GB when the problem
     *          size == 4096 (the linear size of H_local and the mps multisite_tensor).
     *          This means that we filter out
     *              * candidates outside of the energy window (if one is enabled)
     *              * candidates with little or no overlap to the current state.
     *          The filtering process collects the candidate state with best overlap until
     *          either of these two criteria is met:
     *              * More than N=max_accept candidates have been collected
     *              * The subspace error "eps" is low enough¹
     *          The filtering returns the candidates sorted in decreasing overlap.
     *          ¹ We define eps = 1 - Σ_i |<x_i|y>|². A value eps ~1e-10 is reasonable.
     *
     * Step 2)
     *      - (O) If  OptMode::OVERLAP:
     *            Find the index for the best candidate inside of the energy window:
     *              - OA) idx >= 0: Candidate found in window. Return that candidate
     *              - OB) idx = -1: No candidate found in window. Return old tensor
     *
     *      -(V) If OptMode::SUBSPACE or OptMode::SUBSPACE_AND_DIRECT
     *            Find the index for the best candidate inside of the energy window:
     *              - VA) idx = -1: No candidate found in window. Return old tensor
     *              - VB) We have many candidates, including the current tensor.
     *                    They should be viewed as good starting guesses for LBFGS optimization.
     *                    Optimize them one by one, and keep the
     *
     * Step 2)  Find the best overlapping state among the relevant candidates.
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
     *                  1) There are a few candidate states with significant overlap (superposition)
     *              It's clear that we need to optimize, but we have to think carefully about the initial guess.
     *              Right now it makes sense to always choose best overlap theta, since that forces the algorithm to
     *              choose a particular state and not get stuck in superposition. Choosing the old theta may not always
     *              may just entrench the algorithm into a local minima.
     *          D)  If 0 <= best_overlap and best_overlap < overlap_cat
     *              This can happen for three reasons, most often early in the simulation.
     *                  1) There are several candidate states with significant overlap (superposition)
     *                  2) The highest overlapping states were outside of the energy window, leaving just these candidates.
     *                  3) The energy targeting of states has failed for some reason, perhaps the spectrum is particularly dense_lu.
     *              In any case, it is clear we are lost Hilbert space.
     *              Also, the subspace_error is no longer a good measure of how useful the subspace is to us, since it's only
     *              measuring how well the old state can be described, but the old state is likely very different from what
     *              we're looking for.
     *              So to address all three cases, do DIRECT optimization with best_overlap_theta as initial guess.
     *
     * In particular, notice that we never use the candidate that happens to have the best variance.
     */
    // Handy references
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    const auto &edges = *tensors.edges;

    /*
     *  Step 0) Find the subspace
     */

    std::vector<opt_tensor> candidate_list;
    switch(optType) {
        case OptType::CPLX:
            candidate_list = internal::subspace::find_candidates<Scalar>(tensors, settings::precision::min_subspace_error, optMode, optSpace);
            break;
        case OptType::REAL:
            candidate_list = internal::subspace::find_candidates<double>(tensors, settings::precision::min_subspace_error, optMode, optSpace);
            break;
    }

    tools::log->trace("Subspace found with {} eigenvectors", candidate_list.size());

    /*
     *  Step 1) Filter the candidates
     *
     */

    internal::subspace::filter_candidates(candidate_list, settings::precision::min_subspace_error, 128);

    /*
     *
     * Step 2) (Overlap mode)
     *
     */
    if(optMode == OptMode::OVERLAP) {
        auto max_overlap_idx = internal::subspace::get_idx_to_candidate_with_highest_overlap(candidate_list, status.energy_lbound, status.energy_ubound);
        if(max_overlap_idx) {
            // (OA)
            const auto &candidate_max_overlap = *std::next(candidate_list.begin(), static_cast<long>(max_overlap_idx.value()));
            if(tools::log->level() == spdlog::level::trace) {
                tools::log->trace("ceres_subspace_optimization: OA - Candidate {:<2} has highest overlap {:.16f} | energy: {:>20.16f}", max_overlap_idx.value(),
                                  candidate_max_overlap.get_overlap(), candidate_max_overlap.get_energy());
            }
            state.tag_active_sites_have_been_updated(true);
            return candidate_max_overlap;
        } else {
            // (OB)
            tools::log->warn("ceres_subspace_optimization: OB - No overlapping states in energy range. Returning old tensor");
            state.tag_active_sites_have_been_updated(false);
            return initial_tensor;
        }
    }

    /*
     *
     * Step 2) (Variance mode)
     *
     */

    // Construct H² as a matrix (expensive operation!)
    // Also make sure you do this before prepending the current state. All candidates here should be basis vectors
    Eigen::MatrixXcd H2_subspace = internal::get_multisite_hamiltonian_squared_subspace_matrix<Scalar>(model, edges, candidate_list);
    if(optType == OptType::REAL) H2_subspace = H2_subspace.real();
    //    double t_H2_subspace = tools::common::profile::t_opt->get_last_time_interval();

    // Find the best candidates and compute their variance
    auto candidate_list_top_idx = internal::subspace::get_idx_to_candidates_with_highest_overlap(candidate_list, 5, status.energy_lbound, status.energy_ubound);
    for(auto &idx : candidate_list_top_idx) {
        auto &candidate = *std::next(candidate_list.begin(), static_cast<long>(idx));
        candidate.set_variance(tools::finite::measure::energy_variance_per_site(candidate.get_tensor(), tensors));
    }

    // We need the eigenvalues in a convenient format as well
    auto eigvals = internal::subspace::get_eigvals(candidate_list);

    // All the candidates with significant overlap should be considered as initial guesses, including the current theta
    // So we append it here to the list of candidates
    candidate_list.emplace_back(initial_tensor);
    candidate_list_top_idx.emplace_back(candidate_list.size() - 1);

    tools::log->debug("Optimizing with {} initial guesses", candidate_list_top_idx.size());
    std::vector<opt_tensor> optimized_results;
    size_t                  candidate_count = 0;
    for(auto &idx : candidate_list_top_idx) {
        const auto &candidate = *std::next(candidate_list.begin(), static_cast<long>(idx));
        tools::log->trace("Starting LBFGS with candidate {:<2} as initial guess: overlap {:.16f} | energy {:>20.16f} | variance: {:>20.16f} | eigvec {}", idx,
                          candidate.get_overlap(), candidate.get_energy(), std::log10(candidate.get_variance()), candidate.is_basis_vector);
        Eigen::VectorXcd subspace_vector = internal::subspace::get_vector_in_subspace(candidate_list, idx);
        //        Eigen::VectorXcd        theta_new;
        long                    subspace_size = subspace_vector.size();
        double                  energy_new, variance_new;
        [[maybe_unused]] double norm;
        // Note that alpha_i = <theta_initial | theta_new_i> is not supposed to be squared!
        internal::reports::bfgs_add_entry("Subspace", "init", candidate, subspace_size);
        if(tools::log->level() <= spdlog::level::debug) {
            if constexpr(settings::debug) {
                if(tools::log->level() == spdlog::level::trace) {
                    // Check using explicit matrix
                    tools::common::profile::t_chk->tic();
                    Eigen::VectorXcd initial_vector_sbsp = internal::subspace::get_vector_in_subspace(candidate_list, initial_tensor.get_vector());
                    double           overlap_sbsp        = std::abs(subspace_vector.dot(initial_vector_sbsp));
                    Eigen::VectorXcd Hv                  = eigvals.asDiagonal() * subspace_vector;
                    Eigen::VectorXcd H2v                 = H2_subspace.template selfadjointView<Eigen::Upper>() * subspace_vector;
                    Scalar           vHv                 = subspace_vector.dot(Hv);
                    Scalar           vH2v                = subspace_vector.dot(H2v);
                    double           vv                  = subspace_vector.squaredNorm();
                    Scalar           ene                 = vHv / vv;
                    Scalar           var                 = vH2v / vv - ene * ene;
                    double           ene_init_san        = std::real(ene + model.get_energy_reduced()) / static_cast<double>(state.get_length());
                    double           var_init_san        = std::real(var) / static_cast<double>(state.get_length());
                    tools::common::profile::t_chk->toc();
                    std::string description = fmt::format("{:<8} {:<16} {}", "Subspace", candidate.get_name(), "matrix check");
                    internal::reports::bfgs_add_entry(description, subspace_vector.size(), subspace_size, ene_init_san, var_init_san, overlap_sbsp,
                                                      subspace_vector.norm(), 1, 1, tools::common::profile::t_chk->get_last_time_interval());
                }
                if(tools::log->level() == spdlog::level::trace) {
                    // Check current tensor
                    tools::common::profile::t_chk->tic();
                    Eigen::VectorXcd theta_0        = internal::subspace::get_vector_in_fullspace(candidate_list, subspace_vector);
                    auto             theta_0_tensor = Textra::MatrixTensorMap(theta_0, state.active_dimensions());
                    double           energy_0       = tools::finite::measure::energy_per_site(theta_0_tensor, tensors);
                    double           variance_0     = tools::finite::measure::energy_variance_per_site(theta_0_tensor, tensors);
                    double           overlap_0      = std::abs(initial_tensor.get_vector().dot(theta_0));
                    tools::common::profile::t_chk->toc();
                    std::string description = fmt::format("{:<8} {:<16} {}", "Subspace", candidate.get_name(), "fullspace check");
                    internal::reports::bfgs_add_entry(description, theta_0.size(), subspace_size, energy_0, variance_0, overlap_0, theta_0.norm(),
                                                      1, 1, tools::common::profile::t_chk->get_last_time_interval());
                }
            }
        }

        /*
         *
         *  Start the LBFGS optimization process for the subspace
         *
         */

        auto                                  options = internal::ceres_default_options;
        ceres::GradientProblemSolver::Summary summary;
        using namespace tools::finite::opt::internal;
        optimized_results.emplace_back(opt_tensor());
        auto &optimized_tensor = optimized_results.back();
        optimized_tensor.set_name(candidate.get_name());
        optimized_tensor.set_sites(candidate.get_sites());
        tools::common::profile::t_opt_sub_bfgs->tic();
        switch(optType) {
            case OptType::CPLX: {
                Eigen::VectorXd subspace_vector_cast =
                    Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(subspace_vector.data()), 2 * subspace_vector.size());
                auto *                 functor = new ceres_subspace_functor<std::complex<double>>(tensors, status, H2_subspace, eigvals);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running LBFGS subspace cplx");
                ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
                subspace_vector =
                    Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar *>(subspace_vector_cast.data()), subspace_vector_cast.size() / 2).normalized();
                // Copy the results from the functor
                optimized_tensor.set_tensor(subspace::get_vector_in_fullspace(candidate_list, subspace_vector), initial_tensor.get_tensor().dimensions());
                optimized_tensor.set_energy(functor->get_energy());
                optimized_tensor.set_variance(functor->get_variance());
                optimized_tensor.set_iter(summary.iterations.size());
                optimized_tensor.set_counter(functor->get_count());
                optimized_tensor.set_time(summary.total_time_in_seconds);
                optimized_tensor.set_overlap(std::abs(initial_tensor.get_vector().dot(optimized_tensor.get_vector())));
                *tools::common::profile::t_opt_sub_vH2 += *functor->t_vH2;
                *tools::common::profile::t_opt_sub_vH2v += *functor->t_vH2v;
                *tools::common::profile::t_opt_sub_vH += *functor->t_vH;
                *tools::common::profile::t_opt_sub_vHv += *functor->t_vHv;
                break;
            }
            case OptType::REAL: {
                Eigen::VectorXd        subspace_vector_cast = subspace_vector.real();
                Eigen::MatrixXd        H2_subspace_real     = H2_subspace.real();
                auto *                 functor              = new ceres_subspace_functor<double>(tensors, status, H2_subspace_real, eigvals);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running LBFGS subspace real");
                ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
                subspace_vector = subspace_vector_cast.normalized().cast<Scalar>();
                optimized_tensor.set_tensor(subspace::get_vector_in_fullspace(candidate_list, subspace_vector), initial_tensor.get_tensor().dimensions());
                optimized_tensor.set_energy(functor->get_energy());
                optimized_tensor.set_variance(functor->get_variance());
                optimized_tensor.set_iter(summary.iterations.size());
                optimized_tensor.set_counter(functor->get_count());
                optimized_tensor.set_time(summary.total_time_in_seconds);
                optimized_tensor.set_overlap(std::abs(initial_tensor.get_vector().dot(optimized_tensor.get_vector())));
                *tools::common::profile::t_opt_sub_vH2 += *functor->t_vH2;
                *tools::common::profile::t_opt_sub_vH2v += *functor->t_vH2v;
                *tools::common::profile::t_opt_sub_vH += *functor->t_vH;
                *tools::common::profile::t_opt_sub_vHv += *functor->t_vHv;
                break;
            }
        }
        tools::common::profile::t_opt_sub_bfgs->toc();
        reports::time_add_sub_entry();

        tools::log->trace("Finished LBFGS after {} seconds ({} iters). Exit status: {}. Message: {}", summary.total_time_in_seconds, summary.iterations.size(),
                          ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
        //    std::cout << summary.FullReport() << "\n";
        if(tools::log->level() <= spdlog::level::debug) { reports::bfgs_add_entry("Subspace", "opt", optimized_tensor, subspace_size); }

        if(optSpace == OptSpace::SUBSPACE_AND_DIRECT) {
            tools::log->trace("SUBSPACE optimization done. Starting fine tuning with DIRECT optimization");
            tools::common::profile::t_opt_sub->toc(); // Suspend subspace timer
            optimized_tensor = ceres_direct_optimization(tensors, optimized_tensor, status, optType, optMode, optSpace);
            tools::common::profile::t_opt_sub->tic(); // Resume subspace timer
        }
        candidate_count++;
    }

    // Sort thetas in ascending order in variance
    std::sort(optimized_results.begin(), optimized_results.end(),
              [](const opt_tensor &lhs, const opt_tensor &rhs) { return lhs.get_variance() < rhs.get_variance(); });

    // Return the best theta
    tools::common::profile::t_opt_sub->toc();
    return optimized_results.front();
}
