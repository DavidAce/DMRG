//
// Created by david on 2018-02-09.
//

#include "class_xdmrg.h"
#include "class_fdmrg.h"
#include <config/nmspc_settings.h>
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/debug.h>
#include <tools/finite/io.h>
#include <tools/finite/measure.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/opt.h>
#include <tools/finite/svd.h>

class_xdmrg::class_xdmrg(std::shared_ptr<h5pp::File> h5ppFile_) : class_algorithm_finite(std::move(h5ppFile_), AlgorithmType::xDMRG) {
    tools::log->trace("Constructing class_xdmrg");
}

void class_xdmrg::run_preprocessing() {
    tools::log->info("Running {} preprocessing", algo_name);
    class_algorithm_finite::run_preprocessing();
    tools::common::profile::t_pre->tic();
    status.energy_dens_target = settings::xdmrg::energy_density_target;
    status.energy_dens_window = settings::xdmrg::energy_density_window;
    find_energy_range();
    reset_bond_dimension_limits();
    if(settings::input::bitfield >= 0)
        reset_to_random_product_state(ResetReason::INIT);
    else
        reset_to_random_state_in_energy_window(ResetReason::INIT);
    auto spin_components = tools::finite::measure::spin_components(*tensors.state);
    tools::log->info("Initial spin components: {}", spin_components);
    tools::common::profile::t_pre->toc();
    tools::log->info("Finished {} preprocessing", algo_name);
}

void class_xdmrg::run_algorithm() {
    if(state_name.empty())
        state_name = "state_" + std::to_string(state_number);
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", algo_name, enum2str(settings::model::model_type), state_name);
    while(true) {
        single_xDMRG_step();
        print_status_update();
        write_to_file();
        copy_from_tmp();
        check_convergence();
        update_truncation_limit();     // Will update SVD threshold iff the state precision is being limited by truncation error
        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension

        try_projection();
        try_bond_dimension_quench();
        try_disorder_damping();
        try_hamiltonian_perturbation();

        // It's important not to perform the last move.
        // That last state would not get optimized
        if(tensors.state->position_is_any_edge() and not tensors.model->is_perturbed() and not tensors.model->is_damped()) {
            if(status.iter >= settings::xdmrg::max_iters) {
                stop_reason = StopReason::MAX_ITERS;
                break;
            }
            if(status.algorithm_has_succeeded) {
                stop_reason = StopReason::SUCCEEDED;
                break;
            }
            if(status.algorithm_has_to_stop) {
                stop_reason = StopReason::SATURATED;
                break;
            }
            if(status.num_resets > settings::strategy::max_resets) {
                stop_reason = StopReason::MAX_RESET;
                break;
            }
            if(settings::strategy::randomize_early and state_number == 0 and tensors.state->find_largest_chi() >= 32 and
               tools::finite::measure::energy_variance(tensors) < 1e-4) {
                stop_reason = StopReason::RANDOMIZE;
                break;
            }
        }
        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", status.step, status.iter, status.position, status.direction);
        move_center_point();
    }
    tools::log->info("Finished {} simulation of state [{}] -- stop reason: {}", algo_name, state_name, enum2str(stop_reason));
    status.algorithm_has_finished = true;
}

void class_xdmrg::single_xDMRG_step() {
    tools::log->debug("Starting xDMRG step {} | iter {} | pos {} | dir {}", status.step, status.iter, status.position, status.direction);

    using namespace tools::finite;
    using namespace tools::finite::opt;
    tools::common::profile::t_sim->tic();


    //    IDEA:
    //        Try converging to some state from a product state.
    //        We know this is a biased state because it was generated  from a product state of low entanglement.
    //        After converging, apply randomly either identity or pauli_x MPO operators on some sites, to generate
    //        a new initial guess for some other state. Make sure the parity sector is honored. Converging again now
    //        should erase the any biasing that may have occurred due to the selection of initial product state.


    // Setup optpmization mode and space.
    //      - Start from the most general condition
    //      - Override later with narrower conditions


    // Set the fastest mode by default
    opt::OptMode  optMode  = opt::OptMode::VARIANCE;
    opt::OptSpace optSpace = opt::OptSpace::DIRECT;
    opt::OptType  optType  = opt::OptType::CPLX;


    // The first decision is easy. Real or complex optimization
    if(tensors.is_real()) optType = OptType::REAL;

    if(status.chi_lim <= 12) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    if(status.chi_lim <= 8 or status.iter < 2) {
        // TODO: We may not want to do this on states post randomization
        optMode  = OptMode::OVERLAP;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    // Setup strong overrides to normal conditions, e.g.,
    //      - for experiments like perturbation or chi quench
    //      - when the algorithm has already converged

    if(tensors.state->size_2site() < settings::precision::max_size_part_diag) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    if(status.variance_mpo_has_converged) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }

    // Setup the maximum problem size
    long max_problem_size = 0;
    switch(optSpace) {
        case OptSpace::DIRECT: max_problem_size = settings::precision::max_size_direct; break;
        case OptSpace::SUBSPACE_ONLY: max_problem_size = settings::precision::max_size_part_diag; break;
        case OptSpace::SUBSPACE_AND_DIRECT: max_problem_size = settings::precision::max_size_part_diag; break;
    }

    // Make sure not to use SUBSPACE if the 2-site problem size is huge
    // When this happens, we can't use SUBSPACE even with two sites,
    // so we might as well give it to DIRECT, which handles larger problems.
    if(tensors.state->size_2site() > max_problem_size) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }

    // We can optionally make many optimization trials with different
    // number of sites. I.e., if a 2 site optimization did not give
    // any significant improvement we may try with 4 or 8 sites.
    // Here we define a list of trials, where each entry determines
    // how many sites to try.
    std::list<size_t> max_num_sites_list;
    if(status.algorithm_has_stuck_for > 0)
        max_num_sites_list = {4};
    else if(status.algorithm_has_stuck_for > 1)
        max_num_sites_list = {settings::strategy::multisite_max_sites};
    else if(chi_quench_steps > 0)
        max_num_sites_list = {settings::strategy::multisite_max_sites};
    else
        max_num_sites_list = {2};


    // Make sure the site list is sane by sorting and filtering out repeated/invalid entries.
    max_num_sites_list.sort();
    max_num_sites_list.unique();
    max_num_sites_list.remove_if([](auto &elem) { return elem > settings::strategy::multisite_max_sites; });
    if(max_num_sites_list.empty()) throw std::runtime_error("No sites selected for multisite xDMRG");

    tools::log->debug("Possible multisite step sizes: {}", max_num_sites_list);
    size_t theta_count  = 0;
    double variance_old = 1;
    std::map<double, std::tuple<Eigen::Tensor<Scalar, 3>, std::list<size_t>, size_t, double>> results;
    for(auto &max_num_sites : max_num_sites_list) {
        if(optMode == opt::OptMode::OVERLAP and optSpace == opt::OptSpace::DIRECT) throw std::logic_error("OVERLAP mode and DIRECT space are incompatible");

        auto old_num_sites = tensors.active_sites.size();
        auto old_prob_size = tensors.active_problem_size();
        tensors.activate_sites(max_problem_size, max_num_sites);

        // Check that we are not about to solve the same problem again
        if(not results.empty() and tensors.active_sites.size() == old_num_sites and tensors.active_problem_size() == old_prob_size) {
            // If we reached this point we have exhausted the number of sites available
            tools::log->debug("Can't activate more sites");
            break;
        }

        variance_old = measure::energy_variance_per_site(tensors); // Should just take value from cache

        auto   candidate_tensor = opt::find_excited_state(tensors, status, optMode, optSpace, optType);
        double variance_new = measure::energy_variance_per_site(candidate_tensor, tensors);
        results.insert({variance_new, {candidate_tensor, tensors.state->active_sites, theta_count++, tools::common::profile::t_opt->get_last_time_interval()}});

        // We can now decide if we are happy with the result or not.
        double decrease = std::log10(variance_old) / std::log10(variance_new);
        if(decrease < 0.99) {
            tools::log->debug("State improved during {} optimization. Variance decrease: {:.4f}", optSpace, decrease);
            break;
        } else {
            tools::log->debug("State did not improve during {} optimization. Variance decrease {:.4f}", optSpace, decrease);
            continue;
        }
    }
    int result_count = 0;
    for(auto &result : results)
        tools::log->debug("Result {:3} candidate {:3} | variance {:.16f} | time {:.4f} ms", result_count++, std::get<2>(result.second),
                          std::log10(result.first), 1000 * std::get<3>(result.second));


    auto        variance_new = results.begin()->first;
    const auto &multisite_tensor     = std::get<0>(results.begin()->second);
    tensors.state->active_sites      = std::get<1>(results.begin()->second);

    if(std::log10(variance_new) < std::log10(variance_old) - 1e-2) tensors.state->tag_active_sites_have_been_updated(true);

    // Truncate even more if doing chi quench
    //   if(chi_quench_steps > 0) chi_lim = chi_lim_quench_trail;

    // Do the truncation with SVD
    auto variance_before_svd = tools::finite::measure::energy_variance_per_site(multisite_tensor,tensors);
    tools::log->trace("Variance check before SVD: {:.16f}", std::log10(variance_before_svd));
    tensors.merge_multisite_tensor(multisite_tensor, status.chi_lim);
    auto variance_after_svd = tools::finite::measure::energy_variance_per_site(tensors);
    tensors.state->set_truncated_variance((variance_after_svd - variance_before_svd) / variance_after_svd);
    tools::log->trace("Variance check after  SVD: {:.16f}", std::log10(variance_after_svd));
    tools::log->debug("Variance loss due to  SVD: {:.16f}", tensors.state->get_truncated_variance());

    if(settings::precision::use_reduced_energy and tensors.state->position_is_any_edge()) {
        double site_energy = tools::finite::measure::energy_per_site(tensors);
        tools::finite::mpo::reduce_mpo_energy(*tensors.model,site_energy);
    }

    debug::check_integrity(*tensors.state);
    status.energy_dens = (tools::finite::measure::energy_per_site(tensors) - status.energy_min) / (status.energy_max - status.energy_min);


    // Update record holder
    auto var = tools::finite::measure::energy_variance_per_site(tensors);
    if(var < status.lowest_recorded_variance) status.lowest_recorded_variance = var;

    tools::common::profile::t_sim->toc();
    status.wall_time = tools::common::profile::t_tot->get_age();
    status.simu_time = tools::common::profile::t_sim->get_measured_time();
}


void class_xdmrg::single_xDMRG_step_old() {
    tools::log->debug("Starting xDMRG step {} | iter {} | pos {} | dir {}", status.step, status.iter, status.position, status.direction);

    using namespace tools::finite;
    using namespace tools::finite::opt;
    tools::common::profile::t_sim->tic();

    // Set the fastest mode by default
    opt::OptMode  optMode  = opt::OptMode::VARIANCE;
    opt::OptSpace optSpace = opt::OptSpace::DIRECT;
    opt::OptType  optType  = opt::OptType::CPLX;

    //    IDEA:
    //        Try converging to some state from a product state.
    //        We know this is a biased state because it was generated  from a product state of low entanglement.
    //        After converging, apply randomly either identity or pauli_x MPO operators on some sites, to generate
    //        a new initial guess for some other state. Make sure the parity sector is honored. Converging again now
    //        should erase the any biasing that may have occurred due to the selection of initial product state.

    // Setup normal conditions
    if(state_number == 0 and status.chi_lim <= 12) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    if(state_number == 0 and (status.chi_lim <= 8 or status.iter < 2)) {
        optMode  = OptMode::OVERLAP;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    // Setup strong overrides to normal conditions, e.g., for experiments like chi quench

    if(tensors.state->size_2site() < settings::precision::max_size_part_diag and chi_quench_steps > 0) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    if(tensors.state->size_2site() < settings::precision::max_size_part_diag and tools::finite::measure::energy_variance(tensors) > 1e-4) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    //    if(tensors.state->size_2site() < settings::precision::max_size_part_diag and tools::finite::measure::energy_variance(*tensors.state) > 1e-2){
    //        optMode  = OptMode::OVERLAP;
    //        optSpace = OptSpace::SUBSPACE_ONLY;
    //    }

    if(status.variance_mpo_has_converged) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }

    if(tensors.state->is_real()) optType = OptType::REAL;

    long threshold = 0;
    switch(optSpace) {
        case OptSpace::DIRECT: threshold = settings::precision::max_size_direct; break;
        case OptSpace::SUBSPACE_ONLY: threshold = settings::precision::max_size_part_diag; break;
        case OptSpace::SUBSPACE_AND_DIRECT: threshold = settings::precision::max_size_part_diag; break;
    }

    std::list<size_t> max_num_sites_list;
    // Generate a list of maximum number of active sites to try
    if(status.algorithm_has_stuck_for > 0)
        max_num_sites_list = {4};
    else if(status.algorithm_has_stuck_for > 1)
        max_num_sites_list = {settings::strategy::multisite_max_sites};
    else if(chi_quench_steps > 0)
        max_num_sites_list = {settings::strategy::multisite_max_sites};
    else if(optSpace == OptSpace::SUBSPACE_AND_DIRECT)
        max_num_sites_list = {2};
    else if(optSpace == OptSpace::SUBSPACE_ONLY)
        max_num_sites_list = {2};
    else
        max_num_sites_list = {2};

    // Make sure not to use SUBSPACE if the 2-site problem size is huge
    // When this happens, we can't use SUBSPACE even with two sites,
    // so we might as well give it to DIRECT, which handles larger problems.
    if(tensors.state->size_2site() > threshold) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }

    max_num_sites_list.sort();
    max_num_sites_list.unique();
    max_num_sites_list.remove_if([](auto &elem) { return elem > settings::strategy::multisite_max_sites; });
    if(max_num_sites_list.empty()) throw std::runtime_error("No sites selected for multisite xDMRG");

    tools::log->debug("Possible multisite step sizes: {}", max_num_sites_list);
    size_t theta_count  = 0;
    double variance_old = measure::energy_variance_per_site(tensors);
    std::map<double, std::tuple<Eigen::Tensor<Scalar, 3>, std::list<size_t>, size_t, double>> results;
    for(auto &max_num_sites : max_num_sites_list) {
        if(optMode == opt::OptMode::OVERLAP and optSpace == opt::OptSpace::DIRECT) throw std::logic_error("OVERLAP mode and DIRECT space are incompatible");

        auto old_num_sites = tensors.active_sites.size();
        auto old_prob_size = tensors.active_problem_size();
        tensors.activate_sites(threshold, max_num_sites);
        //        if(chi_quench_steps == 0) tensors.state->activate_sites(threshold, max_num_sites);
        //        else tensors.state->activate_truncated_sites(threshold,chi_lim_quench_ahead, max_num_sites);
        //         Reduce bond dimensions for some sites ahead
        //        if(chi_quench_steps > 0) tools::finite::mps::truncate_active_sites(*tensors.state, chi_lim_quench_ahead);

        // Check that we are not about to solve the same problem again
        if(not results.empty() and tensors.active_sites.size() == old_num_sites and tensors.active_problem_size() == old_prob_size) {
            // If we reached this point we have exhausted the number of sites available
            tools::log->debug("Can't activate more sites");
            break;
        }

        auto   theta        = opt::find_excited_state(tensors, status, optMode, optSpace, optType);
        double variance_new = measure::energy_variance_per_site(theta, tensors);
        results.insert({variance_new, {theta, tensors.state->active_sites, theta_count++, tools::common::profile::t_opt->get_last_time_interval()}});
        //
        //        if(state_is_within_energy_window(multisite_tensor)) {
        //            variance_new = measure::energy_variance_per_site(*tensors.state, multisite_tensor);
        //            results.insert({variance_new, {multisite_tensor, tensors.state->active_sites, theta_count++}});
        //        }else{
        //            tools::log->info("Rejecting state found out of energy window");
        //            multisite_tensor = tensors.state->get_multisite_tensor();
        //            variance_new = measure::energy_variance_per_site(*tensors.state);
        //            results.insert({variance_new, {multisite_tensor, tensors.state->active_sites, theta_count++}});
        //        }
        // We can now decide if we are happy with the result or not.
        double decrease = std::log10(variance_old) / std::log10(variance_new);
        if(decrease < 0.99) {
            tools::log->debug("State improved during {} optimization. Variance decrease: {:.4f}", optSpace, decrease);
            break;
        } else {
            tools::log->debug("State did not improve during {} optimization. Variance decrease {:.4f}", optSpace, decrease);
            continue;
        }
    }
    int result_count = 0;
    for(auto &result : results)
        tools::log->debug("Result {:3} candidate {:3} | variance {:.16f} | time {:.4f} ms", result_count++, std::get<2>(result.second),
                          std::log10(result.first), 1000 * std::get<3>(result.second));


    auto        variance_new = results.begin()->first;
    const auto &multisite_tensor     = std::get<0>(results.begin()->second);
    tensors.state->active_sites      = std::get<1>(results.begin()->second);

    if(std::log10(variance_new) < std::log10(variance_old) - 1e-2) tensors.state->tag_active_sites_have_been_updated(true);

    // Truncate multisite_tensor down to chi_lim
    auto chi_lim = status.chi_lim;

    // Truncate even more if doing chi quench
    //    if(chi_quench_steps > 0) chi_lim = chi_lim_quench_trail;

    // Do the truncation with SVD
    auto variance_before_svd = tools::finite::measure::energy_variance_per_site(multisite_tensor,tensors);
    tools::log->trace("Variance check before SVD: {:.16f}", std::log10(variance_before_svd));
    tensors.merge_multisite_tensor(multisite_tensor, status.chi_lim);
    auto variance_after_svd = tools::finite::measure::energy_variance_per_site(tensors);
    tensors.state->set_truncated_variance((variance_after_svd - variance_before_svd) / variance_after_svd);
    tools::log->trace("Variance check after  SVD: {:.16f}", std::log10(variance_after_svd));
    tools::log->debug("Variance loss due to  SVD: {:.16f}", tensors.state->get_truncated_variance());

    if(settings::precision::use_reduced_energy and tensors.state->position_is_any_edge()) {
        double site_energy = tools::finite::measure::energy_per_site(tensors);
        tools::finite::mpo::reduce_mpo_energy(*tensors.model,site_energy);
    }

    debug::check_integrity(*tensors.state);
    status.energy_dens = (tools::finite::measure::energy_per_site(tensors) - status.energy_min) / (status.energy_max - status.energy_min);


    // Update record holder
    auto var = tools::finite::measure::energy_variance_per_site(tensors);
    if(var < status.lowest_recorded_variance) status.lowest_recorded_variance = var;

    tools::common::profile::t_sim->toc();
    status.wall_time = tools::common::profile::t_tot->get_age();
    status.simu_time = tools::common::profile::t_sim->get_measured_time();
}

void class_xdmrg::check_convergence() {
    tools::common::profile::t_con->tic();
    if(tensors.state->position_is_any_edge()) {
        check_convergence_variance();
        check_convergence_entg_entropy();
    }

    //TODO: Move this reset block away from here

    //    status.energy_dens = (tools::finite::measure::energy_per_site(*tensors.state) - status.energy_min ) / (status.energy_max - status.energy_min);
    bool outside_of_window = std::abs(status.energy_dens - status.energy_dens_target) > status.energy_dens_window;
    if(status.iter > 2 and tensors.state->position_is_any_edge()) {
        if(outside_of_window and
           (status.variance_mpo_has_saturated or status.variance_mpo_has_converged or tools::finite::measure::energy_variance_per_site(tensors) < 1e-4)) {
            double      old_energy_dens_window = status.energy_dens_window;
            double      new_energy_dens_window = std::min(energy_window_growth_factor * status.energy_dens_window, 0.5);
            std::string reason = fmt::format("energy {:.16f} saturated outside of energy window {} ± {}", tools::finite::measure::energy_per_site(tensors),
                                             status.energy_dens_target, status.energy_dens_window);
            tools::log->info("Increasing energy window: {} --> {}", old_energy_dens_window, new_energy_dens_window);
            status.energy_dens_window = new_energy_dens_window;
            reset_to_random_state_in_energy_window(ResetReason::SATURATED,settings::strategy::target_sector);
        }
        //        else
        //        if( not     tensors.state->all_sites_updated()
        //            and     status.algorithm_has_got_stuck
        //            and     tools::finite::measure::energy_variance_per_site(*tensors.state) > 1e-4)
        //        {
        //            status.energy_dens_window = std::min(energy_window_growth_factor*status.energy_dens_window, 0.5);
        //            std::string reason = fmt::format("could not update all sites. Energy density: {}, Energy window: {} --> {}",
        //                     status.energy_dens, status.energy_dens_window, std::min(energy_window_growth_factor*status.energy_dens_window, 0.5)
        //                     );
//                     reset_to_random_state_in_energy_window(ResetReason::SATURATED,settings::strategy::target_sector);
        //        }
    }

    status.algorithm_has_converged = status.variance_mpo_has_converged and status.entanglement_has_converged;

    status.algorithm_has_saturated =
        ((status.variance_mpo_saturated_for >= min_saturation_iters and status.entanglement_saturated_for >= min_saturation_iters) or
         (tensors.state->get_iteration() > settings::xdmrg::min_iters and not tensors.state->any_sites_updated()));

    status.algorithm_has_succeeded = status.algorithm_has_converged and status.algorithm_has_saturated;

    status.algorithm_has_got_stuck = status.algorithm_has_saturated and not status.algorithm_has_succeeded;

    if(tensors.state->position_is_any_edge()) {
        status.algorithm_has_stuck_for = status.algorithm_has_got_stuck ? status.algorithm_has_stuck_for + 1 : 0;
    }

    status.algorithm_has_to_stop = status.algorithm_has_stuck_for >= max_stuck_iters;

    if(tensors.state->position_is_any_edge()) {
        tools::log->debug("Simulation has converged: {}", status.algorithm_has_converged);
        tools::log->debug("Simulation has saturated: {}", status.algorithm_has_saturated);
        tools::log->debug("Simulation has succeeded: {}", status.algorithm_has_succeeded);
        tools::log->debug("Simulation has got stuck: {}", status.algorithm_has_got_stuck);
        tools::log->debug("Simulation has stuck for: {}", status.algorithm_has_stuck_for);
        tools::log->debug("Simulation has to stop  : {}", status.algorithm_has_to_stop);
    }

    tools::common::profile::t_con->toc();
}


void class_xdmrg::reset_to_random_state_in_energy_window(ResetReason reason, std::optional<std::string> sector) {
    tools::log->info("Resetting to product state in energy window -- reason: {}", enum2str(reason));
    tools::log->info("Searching for product state in normalized energy range: {} +- {}", status.energy_dens_target, status.energy_dens_window);

    status.num_resets++;
    if(reason == ResetReason::SATURATED and status.num_resets > settings::strategy::max_resets) {
        tools::log->info("Not allowed more resets due to saturation: num resets {} > max resets {}", status.num_resets, settings::strategy::max_resets);
        return;
    }

    int  counter           = 0;
    bool outside_of_window = true;
    tensors.activate_sites(settings::precision::max_size_full_diag,2);
    while(true) {
        reset_to_random_product_state(reason, sector);
        status.energy_dens = tools::finite::measure::energy_normalized(tensors,status.energy_min,status.energy_max);
        outside_of_window  = std::abs(status.energy_dens - status.energy_dens_target) >= status.energy_dens_window;
        tools::log->info("New energy density: {:.16f} | window {} | outside of window: {}", status.energy_dens, status.energy_dens_window, outside_of_window);
        if(not outside_of_window) break;
        counter++;
        if(counter >= 200) throw std::runtime_error(fmt::format("Failed to find initial state in energy window after {}. retries: ", counter));
        if(counter % 10 == 0 and energy_window_growth_factor != 1.0) {
            double old_energy_dens_window = status.energy_dens_window;
            double new_energy_dens_window = std::min(energy_window_growth_factor * status.energy_dens_window, 0.5);

            tools::log->info("Can't find state in energy window.  Increasing energy window: {} --> {}", old_energy_dens_window, new_energy_dens_window);
            status.energy_dens_window = new_energy_dens_window;
        }
    }
    tools::log->info("Energy initial (per site) = {:.16f} | density = {:.8f} | retries = {}", tools::finite::measure::energy_per_site(tensors),
                     status.energy_dens, counter);
    clear_convergence_status();
    status.energy_ubound = status.energy_target + status.energy_dens_window * (status.energy_max - status.energy_min);
    status.energy_lbound = status.energy_target - status.energy_dens_window * (status.energy_max - status.energy_min);
    tools::log->info("Number of product state resets: {}", status.num_resets);
}

void class_xdmrg::find_energy_range() {
    tools::log->trace("Finding energy range");
    // Here we define a set of tasks for fdmrg in order to produce the lowest and highest energy eigenstates,
    // We don't want it to randomize its own model, so we implant our current model before running the tasks.

    std::list<fdmrg_task> gs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_DIM_LIMITS, fdmrg_task::INIT_RANDOM_PRODUCT_STATE,
                                      fdmrg_task::FIND_GROUND_STATE, fdmrg_task::POST_WRITE_RESULT};
    std::list<fdmrg_task> hs_tasks = {fdmrg_task::INIT_CLEAR_STATUS, fdmrg_task::INIT_BOND_DIM_LIMITS, fdmrg_task::INIT_RANDOM_PRODUCT_STATE,
                                      fdmrg_task::FIND_HIGHEST_STATE, fdmrg_task::POST_WRITE_RESULT};
    class_fdmrg fdmrg(h5pp_file);
    *fdmrg.tensors.model = *tensors.model; // Copy the model
    // Find loewst energy state
    tools::log = Logger::setLogger(std::string(enum2str(algo_type)) + "-gs", settings::console::verbosity, settings::console::timestamp);
    fdmrg.run_task_list(gs_tasks);
    status.energy_min = tools::finite::measure::energy_per_site(fdmrg.tensors);

    //Find highest energy state
    tools::log = Logger::setLogger(std::string(enum2str(algo_type)) + "-hs", settings::console::verbosity, settings::console::timestamp);
    fdmrg.run_task_list(hs_tasks);
    status.energy_max = tools::finite::measure::energy_per_site(fdmrg.tensors);

    tools::log = Logger::getLogger(std::string(enum2str(algo_type)));
    status.energy_target = status.energy_min + status.energy_dens_target * (status.energy_max - status.energy_min);
    status.energy_ubound = status.energy_target + status.energy_dens_window * (status.energy_max - status.energy_min);
    status.energy_lbound = status.energy_target - status.energy_dens_window * (status.energy_max - status.energy_min);
    tools::log->info("Energy minimum (per site) = {:.8f}", status.energy_min);
    tools::log->info("Energy maximum (per site) = {:.8f}", status.energy_max);
    tools::log->info("Energy target  (per site) = {:.8f}", status.energy_target);
    tools::log->info("Energy lbound  (per site) = {:.8f}", status.energy_lbound);
    tools::log->info("Energy ubound  (per site) = {:.8f}", status.energy_ubound);

//    exit(0);

//    if(tensors.state->get_length() != settings::model::model_size) throw std::runtime_error("find_energy_range: state length mismatch");
//
//    size_t max_sweeps_during_f_range = 4;
//    // Backup the state and current status
//    class_algorithm_status sim_status_backup = status;
//    class_state_finite     state_backup      = *tensors.state;
//
//    update_bond_dimension_limit(16);
//
//    // Find energy minimum
//    tools::finite::mps::random_product_state(*tensors.state, "random", -1, true);
//    status.iter = tensors.state->reset_iter();
//    status.step = tensors.state->reset_step();
//    while(true) {
//        //        class_algorithm_finite::single_fdmrg_step();
//        print_status_update();
//        // It's important not to perform the last moves. That last state would not get optimized
//        if(tensors.state->position_is_any_edge())
//            if(status.iter >= max_sweeps_during_f_range or tools::finite::measure::energy_variance_per_site(tensors) < 1e-8) break;
//        move_center_point();
//    }
//    double energy_min = tools::finite::measure::energy_per_site(tensors);
//    tensors.normalize_state();
//    write_to_file(StorageReason::EMIN_STATE);
//
//    // Find energy maximum
//    tools::finite::mps::random_product_state(*tensors.state, "random", -1, true);
//    status.iter = tensors.state->reset_iter();
//    status.step = tensors.state->reset_step();
//    while(true) {
//        //        class_algorithm_finite::single_fdmrg_step();
//        print_status_update();
//        // It's important not to perform the last moves. That last state would not get optimized
//        if(tensors.state->position_is_any_edge())
//            if(status.iter >= max_sweeps_during_f_range or tools::finite::measure::energy_variance_per_site(tensors) < 1e-8) break;
//        move_center_point();
//    }
//    double energy_max = tools::finite::measure::energy_per_site(tensors);
//    tensors.normalize_state();
//    write_to_file(StorageReason::EMAX_STATE);
//
//    // Recover the backup
//    status = sim_status_backup;
//    *tensors.state = state_backup;

    // Define energy targets and bounds
//    status.energy_min    = energy_min;
//    status.energy_max    = energy_max;

}

bool class_xdmrg::algo_on() { return settings::xdmrg::on; }
long class_xdmrg::chi_lim_max() { return settings::xdmrg::chi_lim_max; }
// size_t class_xDMRG::write_freq() { return settings::xdmrg::write_freq; }
size_t class_xdmrg::print_freq() { return settings::xdmrg::print_freq; }
bool   class_xdmrg::chi_lim_grow() { return settings::xdmrg::chi_lim_grow; }
long   class_xdmrg::chi_lim_init() { return settings::xdmrg::chi_lim_init; }
bool   class_xdmrg::store_wave_function() { return settings::xdmrg::store_wavefn; }
