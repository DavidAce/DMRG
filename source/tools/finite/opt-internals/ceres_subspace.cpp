//
// Created by david on 2019-07-15.
//
#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include "ceres_subspace_functor.h"
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <iostream>
#include <math/arpack_extra/matrix_product_stl.h>
#include <math/class_eigsolver.h>
#include <math/nmspc_random.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>

using namespace tools::finite::opt;
using namespace tools::finite::opt::internal::subspace;

// template<typename T>
// using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar>
std::tuple<Eigen::MatrixXcd, Eigen::VectorXd> tools::finite::opt::internal::subspace::find_subspace(const class_tensors_finite &tensors,
                                                                                                    double subspace_error_threshold, OptMode optMode,
                                                                                                    OptSpace optSpace) {
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    const auto &edges = *tensors.edges;
    tools::log->trace("Finding subspace");
    using namespace eigutils::eigSetting;
    tools::common::profile::t_ham->tic();
    MatrixType<Scalar> H_local = tools::finite::opt::internal::get_multi_hamiltonian_matrix<Scalar>(model, edges);
    tools::common::profile::t_ham->toc();
    const auto &multisite_tensor = state.get_multisite_tensor();

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
        //        tools::log->debug("Energy target, per site: {}",energy_target/state.get_length());
        tools::log->trace("Energy target + energy reduced = energy per site: {} + {} = {}", energy_target / static_cast<double>(state.get_length()),
                          model.get_energy_reduced() / static_cast<double>(state.get_length()),
                          (energy_target + model.get_energy_reduced()) / static_cast<double>(state.get_length()));

        std::tie(eigvecs, eigvals) = find_subspace_part(H_local, multisite_tensor, energy_target, subspace_error_threshold, optMode, optSpace);
    }
    tools::log->trace("Eigenvalue range: {} --> {}", (eigvals.minCoeff() + model.get_energy_reduced()) / static_cast<double>(state.get_length()),
                      (eigvals.maxCoeff() + model.get_energy_reduced()) / static_cast<double>(state.get_length()));
    //    eigvals = eigvals.array() + model.get_energy_reduced();
    //    tools::log->debug("Eigenvalue range: {} --> {}", eigvals.minCoeff()/state.get_length(),eigvals.maxCoeff()/state.get_length());
    reports::print_eigs_report();

    if constexpr(std::is_same<Scalar, double>::value) {
        Textra::subtract_phase(eigvecs);
        tools::log->trace("truncating imag of eigvecs, sum: {:.16f}", eigvecs.imag().cwiseAbs().sum());
        eigvecs = eigvecs.real();
    }

    return std::make_tuple(eigvecs, eigvals);
}

template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd> tools::finite::opt::internal::subspace::find_subspace<cplx>(const class_tensors_finite &tensors,
                                                                                                                   double  subspace_error_threshold,
                                                                                                                   OptMode optMode, OptSpace optSpace);

template std::tuple<Eigen::MatrixXcd, Eigen::VectorXd> tools::finite::opt::internal::subspace::find_subspace<real>(const class_tensors_finite &tensors,
                                                                                                                   double  subspace_error_threshold,
                                                                                                                   OptMode optMode, OptSpace optSpace);

template<typename Scalar>
std::vector<tools::finite::opt::internal::candidate_tensor> tools::finite::opt::internal::subspace::find_candidates(const class_tensors_finite &tensors,
                                                                                                                    double  subspace_error_threshold,
                                                                                                                    OptMode optMode, OptSpace optSpace) {
    const auto &state = *tensors.state;
    const auto &model = *tensors.model;
    const auto &edges = *tensors.edges;
    tools::log->trace("Finding subspace");
    using namespace eigutils::eigSetting;
    tools::common::profile::t_ham->tic();
    MatrixType<Scalar> H_local = tools::finite::opt::internal::get_multi_hamiltonian_matrix<Scalar>(model, edges);
    tools::common::profile::t_ham->toc();
    const auto &multisite_tensor = state.get_multisite_tensor();

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
        //        tools::log->debug("Energy target, per site: {}",energy_target/state.get_length());
        tools::log->trace("Energy target + energy reduced = energy per site: {} + {} = {}", energy_target / static_cast<double>(state.get_length()),
                          model.get_energy_reduced() / static_cast<double>(state.get_length()),
                          (energy_target + model.get_energy_reduced()) / static_cast<double>(state.get_length()));

        std::tie(eigvecs, eigvals) = find_subspace_part(H_local, multisite_tensor, energy_target, subspace_error_threshold, optMode, optSpace);
    }
    tools::log->trace("Eigenvalue range: {} --> {}", (eigvals.minCoeff() + model.get_energy_reduced()) / static_cast<double>(state.get_length()),
                      (eigvals.maxCoeff() + model.get_energy_reduced()) / static_cast<double>(state.get_length()));
    //    eigvals = eigvals.array() + model.get_energy_reduced();
    //    tools::log->debug("Eigenvalue range: {} --> {}", eigvals.minCoeff()/state.get_length(),eigvals.maxCoeff()/state.get_length());
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

    std::vector<candidate_tensor> candidate_tensors;
    for(long idx = 0; idx < eigvals.size(); idx++) {
        auto tensor_map = Textra::MatrixTensorMap(eigvecs.col(idx), state.active_dimensions());
        // Because we work with reduced energies, the eigenvalues are all close to zero
        // Here we shift the energies back to their unreduced value.
//        auto variance = std::numeric_limits<double>::quiet_NaN();
//        if(optMode == OptMode::VARIANCE)
            // This is an expensive operation so we only do it if we have to
//            variance = tools::finite::measure::energy_variance_per_site(tensor_map, tensors);
        auto overlap = overlaps(idx);
        candidate_tensors.emplace_back(tensor_map, eigvals(idx),energy_reduced, std::nullopt, overlap);
        candidate_tensors.back().is_basis_vector = true;
    }
    return candidate_tensors;
}

template std::vector<tools::finite::opt::internal::candidate_tensor>
    tools::finite::opt::internal::subspace::find_candidates<cplx>(const class_tensors_finite &tensors, double subspace_error_threshold, OptMode optMode,
                                                                  OptSpace optSpace);

template std::vector<tools::finite::opt::internal::candidate_tensor>
    tools::finite::opt::internal::subspace::find_candidates<real>(const class_tensors_finite &tensors, double subspace_error_threshold, OptMode optMode,
                                                                  OptSpace optSpace);

Eigen::Tensor<class_state_finite::Scalar, 3> tools::finite::opt::internal::ceres_subspace_optimization(const class_tensors_finite &  tensors,
                                                                                                       const class_algorithm_status &status, OptType optType,
                                                                                                       OptMode optMode, OptSpace optSpace) {
    tools::log->trace("Optimizing in SUBSPACE mode");
    using namespace eigutils::eigSetting;
    tools::common::profile::t_opt->tic();

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

    double theta_old_energy         = tools::finite::measure::energy_per_site(tensors);
    double theta_old_variance       = tools::finite::measure::energy_variance_per_site(tensors);
    double subspace_error_threshold = settings::precision::min_subspace_error;

    const auto &multisite_mps_tensor = state.get_multisite_tensor();
    const auto  multisite_mps_vector = Eigen::Map<const Eigen::VectorXcd>(multisite_mps_tensor.data(), multisite_mps_tensor.size());
    tools::log->trace("Current energy  : {:.16f}", tools::finite::measure::energy_per_site(tensors));
    tools::log->trace("Current variance: {:.16f}", std::log10(theta_old_variance));

    /*
     *  Step 0) Find the subspace
     */

    std::vector<candidate_tensor> candidate_list;
    switch(optType) {
        case OptType::CPLX: candidate_list = subspace::find_candidates<Scalar>(tensors, subspace_error_threshold, optMode, optSpace); break;
        case OptType::REAL: candidate_list = subspace::find_candidates<double>(tensors, subspace_error_threshold, optMode, optSpace); break;
    }

    tools::log->trace("Subspace found with {} eigenvectors", candidate_list.size());

    /*
     *  Step 1) Filter the candidates
     *
     */

    subspace::filter_candidates(candidate_list, subspace_error_threshold, 128);

    /*
     *
     * Step 2) (Overlap mode)
     *
     */
    if(optMode == OptMode::OVERLAP) {
        auto max_overlap_idx = subspace::get_idx_to_candidate_with_highest_overlap(candidate_list, status.energy_lbound, status.energy_ubound);
        if(max_overlap_idx) {
            // (OA)
            const auto &candidate_max_overlap = *std::next(candidate_list.begin(), static_cast<long>(max_overlap_idx.value()));
            if(tools::log->level() == spdlog::level::trace) {
                tools::log->trace("ceres_subspace_optimization: OA - Candidate {:<2} has highest overlap {:.16f} | energy: {:>20.16f}", max_overlap_idx.value(),
                                  candidate_max_overlap.get_overlap(), candidate_max_overlap.get_energy());
            }
            state.tag_active_sites_have_been_updated(true);
            return candidate_max_overlap.get_tensor();
        } else {
            // (OB)
            tools::log->warn("ceres_subspace_optimization: OB - No overlapping states in energy range. Returning old tensor");
            state.tag_active_sites_have_been_updated(false);
            return multisite_mps_tensor;
        }
    }

    /*
     *
     * Step 2) (Variance mode)
     *
     */



    // Construct H² as a matrix (expensive operation!)
    // Also make sure you do this before prepending the current state. All candidates here should be basis vectors
    tools::common::profile::t_opt->tic();
    Eigen::MatrixXcd H2_subspace = tools::finite::opt::internal::get_multi_hamiltonian_squared_subspace_matrix_new<Scalar>(model, edges, candidate_list);
    if(optType == OptType::REAL) H2_subspace = H2_subspace.real();
    tools::common::profile::t_opt->toc();
    double t_H2_subspace = tools::common::profile::t_opt->get_last_time_interval();

    auto candidate_list_top_idx = subspace::get_idx_to_candidates_with_highest_overlap(candidate_list, 5, status.energy_lbound, status.energy_ubound);

    // We need the eigenvalues in a convenient format as well
    auto eigvals = subspace::get_eigvals(candidate_list);

    // All the candidates with significant overlap should be considered as initial guesses, including the current theta
    // So we append it here to the list of candidates
    reports::bfgs_add_entry("Current state", multisite_mps_tensor.size(), 0, theta_old_energy, std::log10(theta_old_variance), 1.0, multisite_mps_vector.norm(),
                            0, 0, 0.0);
    candidate_list.emplace_back(multisite_mps_tensor, 0.0,theta_old_energy, theta_old_variance, 1.0);
    candidate_list_top_idx.emplace_back(candidate_list.size()-1);


    tools::log->debug("Optimizing mode {} space {} using subspace size {} with {} initial guesses", optMode, optSpace, candidate_list.size(),
                      candidate_list_top_idx.size());

    if(tools::log->level() == spdlog::level::trace){
        for(const auto idx: candidate_list_top_idx){
            auto &candidate = *std::next(candidate_list.begin(), static_cast<long>(idx));
            double var = tools::finite::measure::energy_variance_per_site(candidate.get_tensor(),tensors);
            candidate.set_variance(var);
            tools::log->trace(" - Selected candidate {:<2} as initial guess: overlap {:.16f} | energy {:>20.16f} | variance: {:>20.16f} | eigvec {}", idx,
                              candidate.get_overlap(), candidate.get_energy(), std::log10(candidate.get_variance()), candidate.is_basis_vector);
        }
    }




    std::vector<std::pair<double, Eigen::Tensor<Scalar, 3>>> optimized_results;
    size_t                                                   candidate_count = 0;
    for(auto &idx : candidate_list_top_idx) {
        const auto &candidate = *std::next(candidate_list.begin(), static_cast<long>(idx));
        tools::log->trace("Starting LBFGS with candidate {:<2} as initial guess: overlap {:.16f} | energy {:>20.16f} | variance: {:>20.16f} | eigvec {}", idx,
                          candidate.get_overlap(), candidate.get_energy(), std::log10(candidate.get_variance()), candidate.is_basis_vector);
        Eigen::VectorXcd        theta_start = subspace::get_vector_in_subspace(candidate_list, idx);
        Eigen::VectorXcd        theta_new;
        double                  energy_new, variance_new, overlap_new = 0;
        [[maybe_unused]] double norm;
        // Note that alpha_i = <theta_initial | theta_new_i> is not supposed to be squared!

        if constexpr(settings::debug) {
            if(tools::log->level() == spdlog::level::trace) {
                // Check current tensor
                tools::common::profile::t_opt->tic();
                Eigen::VectorXcd theta_0        = subspace::get_vector_in_fullspace(candidate_list, theta_start);
                auto             theta_0_tensor = Textra::MatrixTensorMap(theta_0, state.active_dimensions());
                double           energy_0       = tools::finite::measure::energy_per_site(theta_0_tensor, tensors);
                double           variance_0     = tools::finite::measure::energy_variance_per_site(theta_0_tensor, tensors);
                double           overlap_0      = std::abs(multisite_mps_vector.dot(theta_0));
                tools::common::profile::t_opt->toc();
                reports::bfgs_add_entry("Candidate " + std::to_string(candidate_count), multisite_mps_tensor.size(), 0, energy_0, std::log10(variance_0),
                                        overlap_0, theta_0.norm(), 1, 1, tools::common::profile::t_opt->get_last_time_interval());
            }
            if(tools::log->level() == spdlog::level::trace) {
                // Check using explicit matrix
                tools::common::profile::t_opt->tic();
                Eigen::VectorXcd theta_0      = subspace::get_vector_in_fullspace(candidate_list, theta_start);
                double           overlap_0    = std::abs(multisite_mps_vector.dot(theta_0));
                Eigen::VectorXcd Hv           = eigvals.asDiagonal() * theta_start;
                Eigen::VectorXcd H2v          = H2_subspace.template selfadjointView<Eigen::Upper>() * theta_start;
                Scalar           vHv          = theta_start.dot(Hv);
                Scalar           vH2v         = theta_start.dot(H2v);
                double           vv           = theta_start.squaredNorm();
                Scalar           ene          = vHv / vv;
                Scalar           var          = vH2v / vv - ene * ene;
                double           ene_init_san = std::real(ene + model.get_energy_reduced()) / static_cast<double>(state.get_length());
                double           var_init_san = std::real(var) / static_cast<double>(state.get_length());
                tools::common::profile::t_opt->toc();
                reports::bfgs_add_entry("Candidate " + std::to_string(candidate_count) + " (matrix)", theta_start.size(), 0, ene_init_san,
                                        std::log10(var_init_san), overlap_0, theta_start.norm(), 1, 1,
                                        tools::common::profile::t_opt->get_last_time_interval());
            }
        }

        /*
         *
         *  Start the LBFGS optimization process for the subspace
         *
         */

        auto                                  options = ceres_default_options;
        ceres::GradientProblemSolver::Summary summary;
        tools::common::profile::t_opt->tic();
        using namespace tools::finite::opt::internal;
        size_t counter = 0, iter = 0;
        switch(optType) {
            case OptType::CPLX: {
                Eigen::VectorXd        theta_start_cast = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(theta_start.data()), 2 * theta_start.size());
                auto *                 functor          = new ceres_subspace_functor<std::complex<double>>(tensors, status, H2_subspace, eigvals);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running L-BFGS");
                ceres::Solve(options, problem, theta_start_cast.data(), &summary);
                iter         = summary.iterations.size();
                counter      = functor->get_count();
                norm         = functor->get_norm();
                energy_new   = functor->get_energy();
                variance_new = functor->get_variance();
                theta_start  = Eigen::Map<Eigen::VectorXcd>(reinterpret_cast<Scalar *>(theta_start_cast.data()), theta_start_cast.size() / 2).normalized();
                theta_new    = subspace::get_vector_in_fullspace(candidate_list, theta_start);
                break;
            }
            case OptType::REAL: {
                Eigen::VectorXd        theta_start_cast = theta_start.real();
                Eigen::MatrixXd        H2_subspace_real = H2_subspace.real();
                auto *                 functor          = new ceres_subspace_functor<double>(tensors, status, H2_subspace_real, eigvals);
                ceres::GradientProblem problem(functor);
                tools::log->trace("Running LBFGS");
                ceres::Solve(options, problem, theta_start_cast.data(), &summary);
                iter         = summary.iterations.size();
                counter      = functor->get_count();
                norm         = functor->get_norm();
                energy_new   = functor->get_energy();
                variance_new = functor->get_variance();
                theta_start  = theta_start_cast.normalized().cast<Scalar>();
                theta_new    = subspace::get_vector_in_fullspace(candidate_list, theta_start).real();
                break;
            }
        }
        tools::common::profile::t_opt->toc();
        reports::time_add_entry(tools::common::profile::t_vH2v->get_measured_time(), tools::common::profile::t_vHv->get_measured_time(),
                                tools::common::profile::t_vH2->get_measured_time(), tools::common::profile::t_vH->get_measured_time(),
                                tools::common::profile::t_op->get_measured_time());

        tools::common::profile::t_vH2v->reset();
        tools::common::profile::t_vHv->reset();
        tools::common::profile::t_vH2->reset();
        tools::common::profile::t_vH->reset();
        tools::common::profile::t_op->reset();

        if(tools::log->level() <= spdlog::level::debug) {
            // Print results of Ceres LBFGS
            overlap_new = (multisite_mps_vector.adjoint() * theta_new).cwiseAbs().sum();
            reports::bfgs_log.emplace_back("LBFGS subspace", multisite_mps_tensor.size(), options.max_lbfgs_rank, energy_new, std::log10(variance_new),
                                           overlap_new, theta_new.norm(), iter, counter, tools::common::profile::t_opt->get_last_time_interval());
        }
        if constexpr (settings::debug) {
            if(tools::log->level() == spdlog::level::trace){
                // Sanity check
                tools::common::profile::t_opt->tic();
                auto   theta_san    = Textra::MatrixTensorMap(theta_new, state.active_dimensions());
                double energy_san   = tools::finite::measure::energy_per_site(theta_san, tensors);
                double variance_san = tools::finite::measure::energy_variance_per_site(theta_san, tensors);
                tools::common::profile::t_opt->toc();
                if(std::abs((variance_san - variance_new) / variance_san) > 0.01)
                    throw std::runtime_error(fmt::format("Variance mismatch in sanity check: {:.16f} != {:.16f}", variance_san, variance_new));
                reports::bfgs_log.emplace_back("Sanity check", theta_san.size(), 0, energy_san, std::log10(variance_san), overlap_new, theta_new.norm(), 1, 1,
                                               tools::common::profile::t_opt->get_last_time_interval());
            }

        }

        tools::log->trace("Finished LBFGS after {} seconds ({} iters). Exit status: {}. Message: {}", summary.total_time_in_seconds, summary.iterations.size(),
                          ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
        //    std::cout << summary.FullReport() << "\n";

        if(optSpace == OptSpace::SUBSPACE_ONLY) {
            auto optimized_theta    = Textra::MatrixTensorMap(theta_new, state.active_dimensions());
            auto optimized_energy   = tools::finite::measure::energy_per_site(optimized_theta, tensors);
            auto optimized_variance = tools::finite::measure::energy_variance_per_site(optimized_theta, tensors);
            auto optimized_vec      = Eigen::Map<const Eigen::VectorXcd>(optimized_theta.data(), optimized_theta.size());
            auto optimized_overlap  = std::abs(multisite_mps_vector.dot(optimized_vec));
            optimized_results.emplace_back(std::make_pair(optimized_variance, optimized_theta));
            reports::bfgs_log.emplace_back("LBFGS subspace", optimized_theta.size(), 0, optimized_energy, std::log10(optimized_variance), optimized_overlap,
                                           optimized_vec.norm(), 1, 1, tools::common::profile::t_opt->get_last_time_interval());
        }
        if(optSpace == OptSpace::SUBSPACE_AND_DIRECT) {
            tools::log->trace("Fine tuning new theta after SUBSPACE optimization");
            auto optimized_theta =
                ceres_direct_optimization(tensors, Textra::MatrixTensorMap(theta_new, state.active_dimensions()), status, optType, optMode, optSpace);
            auto optimized_variance = tools::finite::measure::energy_variance_per_site(optimized_theta, tensors);
            optimized_results.emplace_back(std::make_pair(optimized_variance, optimized_theta));
        }
        candidate_count++;
    }

    // Finish up and print reports
    reports::print_bfgs_report();
    reports::print_time_report();

    // Sort thetas in ascending order in variance
    std::sort(optimized_results.begin(), optimized_results.end(), [](auto &left, auto &right) { return left.first < right.first; });
    // Return the best theta
    return optimized_results.front().second;
}
