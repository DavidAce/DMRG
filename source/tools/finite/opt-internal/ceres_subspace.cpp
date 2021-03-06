//
// Created by david on 2019-07-15.
//
#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include "ceres_subspace_functor.h"
#include <algorithms/class_algorithm_status.h>
#include <ceres/gradient_problem.h>
#include <config/nmspc_settings.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_mps.h>

using namespace tools::finite::opt;
using namespace tools::finite::opt::internal;

template<typename Scalar>
std::vector<opt_mps> internal::subspace::find_candidates(const class_tensors_finite &tensors, double subspace_error_threshold, OptMode optMode,
                                                         OptSpace optSpace) {
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    const auto &edges = *tensors.edges;
    tools::log->trace("Finding subspace");
    double energy_shift = 0.0;
    if(settings::precision::use_reduced_energy and settings::precision::use_shifted_mpo)
        energy_shift = tools::finite::measure::energy_minus_energy_reduced(tensors);

    MatrixType<Scalar> H_local          = tools::finite::opt::internal::get_multisite_hamiltonian_matrix<Scalar>(model, edges, energy_shift);
    const auto &       multisite_tensor = state.get_multisite_mps();
    auto               dbl_length       = static_cast<double>(state.get_length());

    Eigen::MatrixXcd eigvecs;
    Eigen::VectorXd  eigvals;

    // If multitheta is small enough you can afford full diag.
    if(multisite_tensor.size() <= settings::precision::max_size_full_diag) {
        std::tie(eigvecs, eigvals) = find_subspace_full(H_local, multisite_tensor);
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
        std::tie(eigvecs, eigvals) = find_subspace_part(H_local, multisite_tensor, eigval_target, subspace_error_threshold, optMode, optSpace);
    }
    tools::log->trace("Eigval range         : {:.16f} --> {:.16f}", eigvals.minCoeff(), eigvals.maxCoeff());
    tools::log->trace("Energy range         : {:.16f} --> {:.16f}", eigvals.minCoeff() + model.get_energy_reduced(),
                      eigvals.maxCoeff() + model.get_energy_reduced());
    tools::log->trace("Energy range per site: {:.16f} --> {:.16f}", eigvals.minCoeff() / dbl_length + model.get_energy_per_site_reduced(),
                      eigvals.maxCoeff() / dbl_length + model.get_energy_per_site_reduced());
    reports::print_eigs_report();

    if constexpr(std::is_same<Scalar, double>::value) {
        Textra::subtract_phase(eigvecs);
        double trunc = eigvecs.imag().cwiseAbs().sum();
        if(trunc > 1e-12) tools::log->warn("truncating imag of eigvecs, sum: {}", trunc);
        eigvecs = eigvecs.real();
    }

    const auto &    multisite_mps_tensor = state.get_multisite_mps();
    const auto      multisite_mps_vector = Eigen::Map<const Eigen::VectorXcd>(multisite_mps_tensor.data(), multisite_mps_tensor.size());
    auto            energy_reduced       = model.get_energy_reduced();
    Eigen::VectorXd overlaps             = (multisite_mps_vector.adjoint() * eigvecs).cwiseAbs().real();

    double candidate_time = 0;
    for(const auto &item : reports::eigs_log) { candidate_time += item.ham_time + item.lu_time + item.eig_time; }

    std::vector<opt_mps> candidate_list;
    candidate_list.reserve(static_cast<size_t>(eigvals.size()));
    for(long idx = 0; idx < eigvals.size(); idx++) {
        // Important to normalize the eigenvectors that we get from the solver: they are not always well normalized when we get them!
        auto eigvec_i = Textra::TensorCast(eigvecs.col(idx).normalized(), state.active_dimensions());
        candidate_list.emplace_back(fmt::format("eigenvector {}", idx),
                                    eigvec_i,
                                    tensors.active_sites, eigvals(idx), energy_reduced, std::nullopt,
                                    overlaps(idx), tensors.get_length());
        candidate_list.back().set_time(candidate_time);
        candidate_list.back().set_counter(reports::eigs_log.size());
        candidate_list.back().set_iter(reports::eigs_log.size());
        candidate_list.back().is_basis_vector = true;
        candidate_list.back().validate_candidate();
    }
    return candidate_list;
}

template std::vector<opt_mps> internal::subspace::find_candidates<cplx>(const class_tensors_finite &tensors, double subspace_error_threshold, OptMode optMode,
                                                                        OptSpace optSpace);

template std::vector<opt_mps> internal::subspace::find_candidates<real>(const class_tensors_finite &tensors, double subspace_error_threshold, OptMode optMode,
                                                                        OptSpace optSpace);

opt_mps tools::finite::opt::internal::ceres_optimize_subspace(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                              const std::vector<opt_mps> &candidate_list, const Eigen::MatrixXcd &H2_subspace,
                                                              const class_algorithm_status &status, OptType optType, OptMode optMode, OptSpace optSpace) {
    auto t_opt_sub = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub"]->tic_token();

    // Handy references
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;

    opt_mps optimized_mps;
    auto    candidate_best_idx =
        internal::subspace::get_idx_to_candidate_with_highest_overlap(candidate_list, status.energy_llim_per_site, status.energy_ulim_per_site);
    if(candidate_best_idx)
        optimized_mps = *std::next(candidate_list.begin(), static_cast<long>(candidate_best_idx.value()));
    else
        optimized_mps = *candidate_list.begin();
    auto             fullspace_dims  = optimized_mps.get_tensor().dimensions();
    Eigen::VectorXcd subspace_vector = internal::subspace::get_vector_in_subspace(candidate_list, optimized_mps.get_vector());

    tools::log->debug("Optimizing with initial guess: {}", optimized_mps.get_name());

    // We need the eigenvalues in a convenient format as well
    auto eigvals = internal::subspace::get_eigvals(candidate_list);

    /*
     *
     *  Sanity checks
     *
     */
    if(tools::log->level() <= spdlog::level::debug) {
        if constexpr(settings::debug) {
            if(tools::log->level() == spdlog::level::trace) {
                auto t_dbg = tools::common::profile::prof[AlgorithmType::xDMRG]["t_dbg"]->tic_token();
                // Check using explicit matrix
                Eigen::VectorXcd initial_vector_sbsp = internal::subspace::get_vector_in_subspace(candidate_list, initial_mps.get_vector());
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
                std::string      description         = fmt::format("{:<8} {:<16} {}", "Subspace", optimized_mps.get_name(), "matrix check");
                t_dbg.toc();
                internal::reports::bfgs_add_entry(description, subspace_vector.size(), subspace_vector.size(), ene_init_san, var_init_san, overlap_sbsp,
                                                  subspace_vector.norm(), std::numeric_limits<double>::quiet_NaN(),
                                                  std::numeric_limits<double>::quiet_NaN(), 1, 1,
                                                  tools::common::profile::prof[AlgorithmType::xDMRG]["t_dbg"]->get_last_interval());
            }
            if(tools::log->level() == spdlog::level::trace) {
                auto t_dbg = tools::common::profile::prof[AlgorithmType::xDMRG]["t_dbg"]->tic_token();
                // Check current tensor
                Eigen::VectorXcd theta_0        = internal::subspace::get_vector_in_fullspace(candidate_list, subspace_vector);
                auto             theta_0_tensor = Textra::TensorMap(theta_0, state.active_dimensions());
                double           energy_0       = tools::finite::measure::energy_per_site(theta_0_tensor, tensors);
                double           variance_0     = tools::finite::measure::energy_variance(theta_0_tensor, tensors);
                double           overlap_0      = std::abs(initial_mps.get_vector().dot(theta_0));
                std::string      description    = fmt::format("{:<8} {:<16} {}", "Subspace", initial_mps.get_name(), "fullspace check");
                t_dbg.toc();
                internal::reports::bfgs_add_entry(description, theta_0.size(), subspace_vector.size(), energy_0, variance_0, overlap_0, theta_0.norm(),
                                                  std::numeric_limits<double>::quiet_NaN(),
                                                  std::numeric_limits<double>::quiet_NaN(),
                                                  1, 1,
                                                  tools::common::profile::prof[AlgorithmType::xDMRG]["t_dbg"]->get_last_interval());
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
    switch(optType) {
        case OptType::CPLX: {
            auto t_opt_sub_bfgs = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_bfgs"]->tic_token();
            Eigen::VectorXd subspace_vector_cast = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(subspace_vector.data()), 2 * subspace_vector.size());
            auto *          functor              = new ceres_subspace_functor<std::complex<double>>(tensors, status, H2_subspace, eigvals);
            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS subspace cplx");
            ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
            subspace_vector =
                Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar *>(subspace_vector_cast.data()), subspace_vector_cast.size() / 2).normalized();
            // Copy the results from the functor
            optimized_mps.set_tensor(subspace::get_tensor_in_fullspace(candidate_list, subspace_vector, fullspace_dims));
            optimized_mps.set_counter(functor->get_count());
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2"] += *functor->t_H2n;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2v"] += *functor->t_nH2n;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH"] += *functor->t_Hn;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vHv"] += *functor->t_nHn;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_step"] += *functor->t_step;
            break;
        }
        case OptType::REAL: {
            auto t_opt_sub_bfgs = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_bfgs"]->tic_token();
            Eigen::VectorXd   subspace_vector_cast = subspace_vector.real();
            Eigen::MatrixXd   H2_subspace_real     = H2_subspace.real();
            auto *            functor              = new ceres_subspace_functor<double>(tensors, status, H2_subspace_real, eigvals);
            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS subspace real");
            ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
            subspace_vector = subspace_vector_cast.normalized().cast<Scalar>();
            optimized_mps.set_tensor(subspace::get_tensor_in_fullspace(candidate_list, subspace_vector, fullspace_dims));
            optimized_mps.set_counter(functor->get_count());
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2"] += *functor->t_H2n;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2v"] += *functor->t_nH2n;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH"] += *functor->t_Hn;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vHv"] += *functor->t_nHn;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_step"] += *functor->t_step;
            break;
        }
    }
    // Copy and set the rest of the tensor metadata
    optimized_mps.normalize();
    optimized_mps.set_energy(tools::finite::measure::energy(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_variance(tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_overlap(std::abs(initial_mps.get_vector().dot(optimized_mps.get_vector())));
    optimized_mps.set_iter(summary.iterations.size());
    optimized_mps.set_time(summary.total_time_in_seconds);

    if constexpr(settings::debug) {
        auto t_dbg = tools::common::profile::prof[AlgorithmType::xDMRG]["t_dbg"]->tic_token();
        // Check that Ceres results are correct
        double energy_check   = tools::finite::measure::energy(optimized_mps.get_tensor(), tensors);
        double variance_check = tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors);
        if(std::abs(1.0 - std::abs(optimized_mps.get_energy() / energy_check)) > 1e-3)
            tools::log->warn("Energy mismatch: Ceres: {:.16f} | DMRG {:.16f}", optimized_mps.get_energy(), energy_check);
        if(std::abs(1.0 - std::abs(optimized_mps.get_variance() / variance_check)) > 1e-3)
            tools::log->warn("Variance mismatch: Ceres: {:.16f} | DMRG {:.16f}", std::log10(optimized_mps.get_variance()), std::log10(variance_check));
    }

    reports::time_add_sub_entry();
    int    hrs = static_cast<int>(summary.total_time_in_seconds / 3600);
    int    min = static_cast<int>(std::fmod(summary.total_time_in_seconds, 3600) / 60);
    double sec = std::fmod(std::fmod(summary.total_time_in_seconds, 3600), 60);
    tools::log->debug("Finished LBFGS in {:0<2}:{:0<2}:{:0<.1f} seconds and {} iters. Exit status: {}. Message: {}", hrs, min, sec, summary.iterations.size(),
                      ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
    //    std::cout << summary.FullReport() << "\n";
    if(tools::log->level() <= spdlog::level::debug) { reports::bfgs_add_entry("Subspace", "opt", optimized_mps, subspace_vector.size()); }

    if(optSpace == OptSpace::SUBSPACE_AND_DIRECT) {
        tools::log->trace("SUBSPACE optimization done. Starting fine tuning with DIRECT optimization");
        t_opt_sub.toc();// Suspend subspace timer
        optimized_mps = ceres_direct_optimization(tensors, optimized_mps, status, optType, optMode, optSpace);
        t_opt_sub.tic(); // Resume subspace timer
    }
    // Return the optimized result
    return optimized_mps;

}

opt_mps tools::finite::opt::internal::ceres_subspace_optimization(const class_tensors_finite &tensors, const class_algorithm_status &status, OptType optType,
                                                                  OptMode optMode, OptSpace optSpace) {
    std::vector<size_t> sites(tensors.active_sites.begin(), tensors.active_sites.end());
    opt_mps             initial_mps("current state", tensors.get_multisite_mps(), sites,
                        tools::finite::measure::energy(tensors) - tensors.model->get_energy_reduced(), // Eigval
                        tensors.model->get_energy_reduced(),                                           // Energy reduced for full system
                        tools::finite::measure::energy_variance(tensors),
                        1.0, // Overlap
                        tensors.get_length());
    initial_mps.validate_candidate();
    return ceres_subspace_optimization(tensors, initial_mps, status, optType, optMode, optSpace);
}

opt_mps tools::finite::opt::internal::ceres_subspace_optimization(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                                  const class_algorithm_status &status, OptType optType, OptMode optMode, OptSpace optSpace) {
    tools::log->trace("Optimizing in SUBSPACE mode");
    auto t_opt_sub = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub"]->tic_token();


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
     *            Find the index for the best candidate inside of the energy window:
     *              - OA) idx >= 0: Candidate found in window. Return that candidate
     *              - OB) idx = -1: No candidate found in window. Return old tensor
     *       -(V) If OptMode::SUBSPACE or OptMode::SUBSPACE_AND_DIRECT
     *          Make sure k is as small as possible. I.e. filter out eigenvectors from [x] down
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
     *          Find the index for the best candidate inside of the energy window:
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
     *              choose a particular state and not get stuck in superposition. Choosing the old theta may just entrench
     *              the algorithm into a local minima.
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
     *  Step 0) Find the subspace.
     *  The subspace set of candidate eigenstates obtained from full or partial diagonalization
     */

    std::vector<opt_mps> candidate_list;
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
     *
     * Step 1) (Overlap mode) Return best overlap
     *
     */
    if(optMode == OptMode::OVERLAP) {
        auto max_overlap_idx =
            internal::subspace::get_idx_to_candidate_with_highest_overlap(candidate_list, status.energy_llim_per_site, status.energy_ulim_per_site);
        if(max_overlap_idx) {
            // (OA)
            //            if constexpr(settings::debug) {
            //                for(const auto & [idx, candidate] : iter::enumerate(candidate_list)) {
            //                    candidate.set_variance(tools::finite::measure::energy_variance(candidate.get_tensor(), tensors));
            //                    std::string msg = fmt::format("Candidate {:10} | overlap {:<14.12f} | energy {:<+20.16f} | variance {:<+20.16f}",
            //                    candidate.get_name(),
            //                                                  candidate.get_overlap(), candidate.get_energy_per_site(), std::log10(candidate.get_variance()));
            //                    if(idx == max_overlap_idx) msg.append("   <--- max overlap");
            //                    tools::log->trace(msg);
            //                }
            //            }
            auto &candidate_max_overlap = *std::next(candidate_list.begin(), static_cast<long>(max_overlap_idx.value()));
            candidate_max_overlap.set_variance(tools::finite::measure::energy_variance(candidate_max_overlap.get_tensor(), tensors));
            if(tools::log->level() == spdlog::level::trace) {
                tools::log->trace("ceres_subspace_optimization: Candidate {:<2} has highest overlap {:.16f} | energy {:>20.16f} | variance {:>20.16f}",
                                  max_overlap_idx.value(), candidate_max_overlap.get_overlap(), candidate_max_overlap.get_energy_per_site(),
                                  std::log10(candidate_max_overlap.get_variance()));
            }
            if(candidate_max_overlap.get_overlap() < 0.1)
                tools::log->debug("ceres_subspace_optimization: Overlap fell below < 0.1: {:20.16f}", candidate_max_overlap.get_overlap());
            return candidate_max_overlap;
        } else {
            // (OB)
            tools::log->warn("ceres_subspace_optimization: No overlapping states in energy range. Returning old tensor");
            return initial_mps;
        }
    }

    /*
     *  Step 1)  (Variance mode) Filter the candidates
     *
     */

    internal::subspace::filter_candidates(candidate_list, settings::precision::min_subspace_error, settings::precision::max_subspace_size);

    /*
     *
     * Step 2) (Variance mode)
     *
     */

    // Construct H² as a matrix (expensive operation!)
    // Also make sure you do this before prepending the current state. All candidates here should be basis vectors
    double energy_shift = 0.0;
    if(settings::precision::use_reduced_energy and settings::precision::use_shifted_mpo and not model.is_compressed_mpo_squared())
        energy_shift = tools::finite::measure::energy_minus_energy_reduced(tensors);
    Eigen::MatrixXcd H2_subspace = internal::get_multisite_hamiltonian_squared_subspace_matrix<Scalar>(model, edges, candidate_list, std::pow(energy_shift,2));
    if(optType == OptType::REAL) H2_subspace = H2_subspace.real();

    // Find the best candidates and compute their variance
    auto candidate_list_top_idx =
        internal::subspace::get_idx_to_candidates_with_highest_overlap(candidate_list, 5, status.energy_llim_per_site, status.energy_ulim_per_site);
    for(auto &idx : candidate_list_top_idx) {
        auto &candidate = *std::next(candidate_list.begin(), static_cast<long>(idx));
        candidate.set_variance(tools::finite::measure::energy_variance(candidate.get_tensor(), tensors));
    }

    // We need the eigenvalues in a convenient format as well
    auto eigvals = internal::subspace::get_eigvals(candidate_list);

    // All the candidates with significant overlap should be considered as initial guesses, including the current theta
    // So we append it here to the list of candidates
    candidate_list.emplace_back(initial_mps);
    candidate_list_top_idx.emplace_back(candidate_list.size() - 1);

    tools::log->debug("Optimizing with {} initial guesses", candidate_list_top_idx.size());
    std::vector<opt_mps> optimized_results;
    size_t               candidate_count = 0;
    for(auto &idx : candidate_list_top_idx) {
        const auto &candidate = *std::next(candidate_list.begin(), static_cast<long>(idx));
        tools::log->trace("Starting LBFGS with candidate {:<2} as initial guess: overlap {:.16f} | energy {:>20.16f} | variance: {:>20.16f} | eigvec {}", idx,
                          candidate.get_overlap(), candidate.get_energy_per_site(), std::log10(candidate.get_variance()), candidate.is_basis_vector);
        Eigen::VectorXcd        subspace_vector = internal::subspace::get_vector_in_subspace(candidate_list, idx);
        long                    subspace_size   = subspace_vector.size();
        [[maybe_unused]] double norm;
        internal::reports::bfgs_add_entry("Subspace", "init", candidate, subspace_size);
        if(tools::log->level() <= spdlog::level::debug) {
            if constexpr(settings::debug) {
                auto t_dbg = tools::common::profile::prof[AlgorithmType::xDMRG]["t_dbg"]->tic_token();
                if(tools::log->level() == spdlog::level::trace) {
                    // Check using explicit matrix
                    Eigen::VectorXcd initial_vector_sbsp = internal::subspace::get_vector_in_subspace(candidate_list, initial_mps.get_vector());
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
                    std::string description = fmt::format("{:<8} {:<16} {}", "Subspace", candidate.get_name(), "matrix check");
                    internal::reports::bfgs_add_entry(description, subspace_vector.size(), subspace_size, ene_init_san, var_init_san, overlap_sbsp,
                                                      std::numeric_limits<double>::quiet_NaN(),
                                                      std::numeric_limits<double>::quiet_NaN(),
                                                      subspace_vector.norm(), 1, 1,
                                                      tools::common::profile::prof[AlgorithmType::xDMRG]["t_dbg"]->get_last_interval());
                }
                if(tools::log->level() == spdlog::level::trace) {
                    // Check current tensor
                    Eigen::VectorXcd theta_0        = internal::subspace::get_vector_in_fullspace(candidate_list, subspace_vector);
                    auto             theta_0_tensor = Textra::TensorMap(theta_0, state.active_dimensions());
                    double           energy_0       = tools::finite::measure::energy_per_site(theta_0_tensor, tensors);
                    double           variance_0     = tools::finite::measure::energy_variance(theta_0_tensor, tensors);
                    double           overlap_0      = std::abs(initial_mps.get_vector().dot(theta_0));
                    std::string description = fmt::format("{:<8} {:<16} {}", "Subspace", candidate.get_name(), "fullspace check");
                    internal::reports::bfgs_add_entry(description, theta_0.size(), subspace_size, energy_0, variance_0, overlap_0, theta_0.norm(),
                                                      std::numeric_limits<double>::quiet_NaN(),
                                                      std::numeric_limits<double>::quiet_NaN(),
                                                      1, 1,
                                                      tools::common::profile::prof[AlgorithmType::xDMRG]["t_dbg"]->get_last_interval());
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
        optimized_results.emplace_back(opt_mps());
        auto &optimized_mps = optimized_results.back();
        optimized_mps.set_name(candidate.get_name());
        optimized_mps.set_sites(candidate.get_sites());
        optimized_mps.set_length(candidate.get_length());
        optimized_mps.set_energy_reduced(candidate.get_energy_reduced());
        switch(optType) {
            case OptType::CPLX: {
                auto t_opt_sub_bfgs = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_bfgs"]->tic_token();
                Eigen::VectorXd subspace_vector_cast =
                    Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(subspace_vector.data()), 2 * subspace_vector.size());
                auto *            functor = new ceres_subspace_functor<std::complex<double>>(tensors, status, H2_subspace, eigvals);
                CustomLogCallback ceres_logger(*functor);
                options.callbacks.emplace_back(&ceres_logger);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running LBFGS subspace cplx");
                ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
                subspace_vector =
                    Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar *>(subspace_vector_cast.data()), subspace_vector_cast.size() / 2).normalized();
                // Copy the results from the functor
                optimized_mps.set_tensor(subspace::get_vector_in_fullspace(candidate_list, subspace_vector), initial_mps.get_tensor().dimensions());
                optimized_mps.set_counter(functor->get_count());
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2"] += *functor->t_H2n;
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2v"] += *functor->t_nH2n;
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH"] += *functor->t_Hn;
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vHv"] += *functor->t_nHn;
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_step"] += *functor->t_step;
                break;
            }
            case OptType::REAL: {
                auto t_opt_sub_bfgs = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_bfgs"]->tic_token();
                Eigen::VectorXd   subspace_vector_cast = subspace_vector.real();
                Eigen::MatrixXd   H2_subspace_real     = H2_subspace.real();
                auto *            functor              = new ceres_subspace_functor<double>(tensors, status, H2_subspace_real, eigvals);
                CustomLogCallback ceres_logger(*functor);
                options.callbacks.emplace_back(&ceres_logger);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running LBFGS subspace real");
                ceres::Solve(options, problem, subspace_vector_cast.data(), &summary);
                subspace_vector = subspace_vector_cast.normalized().cast<Scalar>();
                optimized_mps.set_tensor(subspace::get_vector_in_fullspace(candidate_list, subspace_vector), initial_mps.get_tensor().dimensions());
                optimized_mps.set_counter(functor->get_count());
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2"] += *functor->t_H2n;
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2v"] += *functor->t_nH2n;
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH"] += *functor->t_Hn;
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vHv"] += *functor->t_nHn;
                *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_step"] += *functor->t_step;
                break;
            }
        }
        // Copy and set the rest of the tensor metadata
        optimized_mps.normalize();
        optimized_mps.set_energy(tools::finite::measure::energy(optimized_mps.get_tensor(), tensors));
        optimized_mps.set_variance(tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors));
        optimized_mps.set_overlap(std::abs(initial_mps.get_vector().dot(optimized_mps.get_vector())));
        optimized_mps.set_iter(summary.iterations.size());
        optimized_mps.set_time(summary.total_time_in_seconds);
        optimized_mps.set_overlap(std::abs(initial_mps.get_vector().dot(optimized_mps.get_vector())));

        if constexpr(settings::debug) {
            auto t_dbg = tools::common::profile::prof[AlgorithmType::xDMRG]["t_dbg"]->tic_token();
            // Check that Ceres results are correct
            double energy_check   = tools::finite::measure::energy(optimized_mps.get_tensor(), tensors);
            double variance_check = tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors);
            if(std::abs(1.0 - std::abs(optimized_mps.get_energy() / energy_check)) > 1e-3)
                tools::log->warn("Energy mismatch: Ceres: {:.16f} | DMRG {:.16f}", optimized_mps.get_energy(), energy_check);
            if(std::abs(1.0 - std::abs(optimized_mps.get_variance() / variance_check)) > 1e-3)
                tools::log->warn("Variance mismatch: Ceres: {:.16f} | DMRG {:.16f}", std::log10(optimized_mps.get_variance()), std::log10(variance_check));
        }

        reports::time_add_sub_entry();
        int    hrs = static_cast<int>(summary.total_time_in_seconds / 3600);
        int    min = static_cast<int>(std::fmod(summary.total_time_in_seconds, 3600) / 60);
        double sec = std::fmod(std::fmod(summary.total_time_in_seconds, 3600), 60);
        tools::log->debug("Finished LBFGS in {:0<2}:{:0<2}:{:0<.1f} seconds and {} iters. Exit status: {}. Message: {}", hrs, min, sec,
                          summary.iterations.size(), ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
        //    std::cout << summary.FullReport() << "\n";
        if(tools::log->level() <= spdlog::level::debug) { reports::bfgs_add_entry("Subspace", "opt", optimized_mps, subspace_size); }

        if(optSpace == OptSpace::SUBSPACE_AND_DIRECT) {
            tools::log->trace("SUBSPACE optimization done. Starting fine tuning with DIRECT optimization");
            t_opt_sub.toc(); // Suspend subspace timer
            optimized_mps = ceres_direct_optimization(tensors, optimized_mps, status, optType, optMode, optSpace);
            t_opt_sub.tic(); // Resume subspace timer
        }
        candidate_count++;
    }

    // Sort thetas in ascending order in variance
    std::sort(optimized_results.begin(), optimized_results.end(),
              [](const opt_mps &lhs, const opt_mps &rhs) { return lhs.get_variance() < rhs.get_variance(); });

    // Return the best theta
    return optimized_results.front();
}
