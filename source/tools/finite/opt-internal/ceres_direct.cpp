//
// Created by david on 2019-07-09.
//

#include "ceres_direct_functor.h"
#include <algorithms/class_algorithm_status.h>
#include <ceres/gradient_problem.h>
#include <config/nmspc_settings.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/fmt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_mps.h>
tools::finite::opt::opt_mps tools::finite::opt::internal::ceres_direct_optimization(const class_tensors_finite &tensors, const class_algorithm_status &status,
                                                                                      OptType optType, OptMode optMode, OptSpace optSpace) {
    std::vector<size_t> sites(tensors.active_sites.begin(), tensors.active_sites.end());
    opt_mps             initial_state("current state", tensors.state->get_multisite_mps(), sites,
                            tools::finite::measure::energy(tensors) - tensors.model->get_energy_reduced(), // Eigval
                            tensors.model->get_energy_reduced(),                                           // Energy reduced for full system
                            tools::finite::measure::energy_variance(tensors),
                            1.0, // Overlap
                            tensors.get_length());

    return ceres_direct_optimization(tensors, initial_state, status, optType, optMode, optSpace);
}

tools::finite::opt::opt_mps tools::finite::opt::internal::ceres_direct_optimization(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                                                      const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                                      OptSpace optSpace) {
    tools::log->trace("Optimizing in DIRECT mode");
    auto t_opt_dir = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir"]->tic_token();
    if constexpr(settings::debug)
        if(initial_mps.has_nan()) throw std::runtime_error("initial_mps has nan's");
    reports::bfgs_add_entry("Direct", "init", initial_mps);
    const auto &current_mps = tensors.state->get_multisite_mps();
    const auto  current_map = Eigen::Map<const Eigen::VectorXcd>(current_mps.data(), current_mps.size());

    auto      options = internal::ceres_default_options;
    auto      summary = ceres::GradientProblemSolver::Summary();
    opt_mps   optimized_mps;
    optimized_mps.set_name(initial_mps.get_name());
    optimized_mps.set_sites(initial_mps.get_sites());
    optimized_mps.set_length(initial_mps.get_length());
    optimized_mps.set_energy_reduced(initial_mps.get_energy_reduced());

    // When we use "reduced-energy mpo's", the current energy per site is subtracted from all the MPO's at the beginning of every iteration.
    // The subtracted number is called "energy_reduced".
    // Directly after the energy subtraction, the target energy is exactly zero, and Var H = <(H-Er)²>
    // However, after some steps in the same iteration, the optimizer may have found a state with slightly different energy,
    // so the target energy is actually 0 + dE.
    // We can obtain the shift amount dE = <H-Er>, which we call "eigval", i.e. the eigenvalue of the operator <H-Er>
    // Then, technically Var H = <(H-E+|dE|)²>, and if unaccounted for, we may get Var H < 0, which is a real pain since
    // we are optimizing the logarithm of Var H.
    // Here we use functor->set_shift in order to account for the shifted energy and reach better precision.
    // Note 1: shifting only works if the mpo is not already compressed
    // Note 2: shifting only makes sense if we are using reduced-energy mpos

    switch(optType) {
        case OptType::CPLX: {
            auto t_opt_dir_bfgs = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_bfgs"]->tic_token();
            // Copy the initial guess and operate directly on it
            optimized_mps.set_tensor(initial_mps.get_tensor());
            auto *            functor = new ceres_direct_functor<std::complex<double>>(tensors, status);

            if(settings::precision::use_reduced_energy and not tensors.model->is_compressed_mpo_squared())
                functor->set_shift(-std::abs(initial_mps.get_eigval())); // Account for the shange in energy since the last energy reduction
            functor->compress(); // Compress the virtual bond between MPO² and the environments

            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS direct cplx");
            ceres::Solve(options, problem, optimized_mps.get_vector_cplx_as_2xreal().data(), &summary);
            // Copy the results from the functor
            optimized_mps.set_counter(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_grad_norm(functor->get_grad_max_norm());
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vH2"] += *functor->t_vH2;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vH2v"] += *functor->t_vH2v;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vH"] += *functor->t_vH;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vHv"] += *functor->t_vHv;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_step"] += *functor->t_step;
            break;
        }
        case OptType::REAL: {
            auto t_opt_dir_bfgs = tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_bfgs"]->tic_token();
            // Here we make a temporary
            auto              initial_state_real = initial_mps.get_vector_cplx_as_1xreal();
            auto *            functor            = new ceres_direct_functor<double>(tensors, status);

            if(settings::precision::use_reduced_energy and not tensors.model->is_compressed_mpo_squared())
                functor->set_shift(-std::abs(initial_mps.get_eigval())); // Account for the shange in energy since the last energy reduction
            functor->compress(); // Compress the virtual bond between MPO² and the environments

            CustomLogCallback ceres_logger(*functor);
            options.callbacks.emplace_back(&ceres_logger);
            ceres::GradientProblem problem(functor);
            tools::log->trace("Running LBFGS direct real");
            if constexpr(settings::debug) if(initial_state_real.hasNaN()) throw std::runtime_error("initial_state_real has nan's");
            ceres::Solve(options, problem, initial_state_real.data(), &summary);
            // Copy the results from the functor
            optimized_mps.set_counter(functor->get_count());
            optimized_mps.set_delta_f(functor->get_delta_f());
            optimized_mps.set_grad_norm(functor->get_grad_max_norm());
            optimized_mps.set_tensor_real(initial_state_real.data(), initial_mps.get_tensor().dimensions());
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vH2"] += *functor->t_vH2;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vH2v"] += *functor->t_vH2v;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vH"] += *functor->t_vH;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vHv"] += *functor->t_vHv;
            *tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_step"] += *functor->t_step;
            break;
        }
    }
    // Copy and set the rest of the tensor metadata
    optimized_mps.normalize();
    optimized_mps.set_energy(tools::finite::measure::energy(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_variance(tools::finite::measure::energy_variance(optimized_mps.get_tensor(), tensors));
    optimized_mps.set_overlap(std::abs(current_map.dot(optimized_mps.get_vector())));
    optimized_mps.set_iter(summary.iterations.size());
    optimized_mps.set_time(summary.total_time_in_seconds);

    reports::time_add_dir_entry();
    int    hrs = static_cast<int>(summary.total_time_in_seconds / 3600);
    int    min = static_cast<int>(std::fmod(summary.total_time_in_seconds, 3600) / 60);
    double sec = std::fmod(std::fmod(summary.total_time_in_seconds, 3600), 60);
    tools::log->debug("Finished LBFGS in {:0<2}:{:0<2}:{:0<.1f} seconds and {} iters. Exit status: {}. Message: {}", hrs, min, sec, summary.iterations.size(),
                      ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
    reports::bfgs_add_entry("Direct", "opt", optimized_mps);

    tools::log->trace("Returning theta from optimization mode {} space {}", optMode, optSpace);
    return optimized_mps;
}
