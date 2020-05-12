//
// Created by david on 2018-02-09.
//

#include "class_xDMRG.h"
#include "class_fDMRG.h"
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <model/class_model_finite.h>
#include <edges/class_edges_finite.h>
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

class_xDMRG::class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_) : class_algorithm_finite(std::move(h5ppFile_), SimulationType::xDMRG) {
    tools::log->trace("Constructing class {}", sim_name);
}

void class_xDMRG::run_preprocessing() {
    tools::log->info("Running {} preprocessing", sim_name);
    class_algorithm_finite::run_preprocessing();
    tools::common::profile::t_pre->tic();
    sim_status.energy_dens_target = settings::xdmrg::energy_density_target;
    sim_status.energy_dens_window = settings::xdmrg::energy_density_window;
    find_energy_range();
    state->set_chi_max(chi_max());
    sim_status.chi_max = chi_max();
    update_bond_dimension_limit(chi_init());
    if(settings::input::bitfield >= 0)
        reset_to_initial_state();
    else
        reset_to_random_state_in_energy_window(settings::strategy::initial_parity_sector, false, "Initializing");
    auto spin_components = tools::finite::measure::spin_components(*state);
    tools::log->info("Initial spin components: {}", spin_components);

    tools::common::profile::t_pre->toc();
    tools::log->info("Finished {} preprocessing", sim_name);
}

void class_xDMRG::run_simulation() {
    sim_status.clear();
    state_name = "state_" + std::to_string(state_number);
    tools::log->info("Starting {} simulation of model [{}] for state [{}]", sim_name, enum2str(settings::model::model_type), state_name);
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
        if(state->position_is_any_edge() and not model->is_perturbed() and not model->is_damped()) {
            if(sim_status.iter >= settings::xdmrg::max_iters) {
                stop_reason = StopReason::MAX_ITERS;
                break;
            }
            if(sim_status.simulation_has_succeeded) {
                stop_reason = StopReason::SUCCEEDED;
                break;
            }
            if(sim_status.simulation_has_to_stop) {
                stop_reason = StopReason::SATURATED;
                break;
            }
            if(sim_status.num_resets > settings::precision::max_resets) {
                stop_reason = StopReason::MAX_RESET;
                break;
            }
            if(settings::strategy::randomize_early and state_number == 0 and state->find_largest_chi() >= 32 and
               tools::finite::measure::energy_variance(*state) < 1e-4) {
                stop_reason = StopReason::RANDOMIZE;
                break;
            }
        }
        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", sim_status.step, sim_status.iter, state->get_position(), state->get_direction());
        move_center_point();
    }
    tools::log->info("Finished {} simulation of state [{}] -- status: {}", sim_name, state_name, enum2str(stop_reason));
//    if(++sim_status.state_number >= settings::xdmrg::max_states)
//        break;
//    else {
//        run_postprocessing(); // Saves and prints full status update
//        reset_to_random_current_state(32);
//        state_name = "state_" + std::to_string(sim_status.state_number);
//    }
}

void class_xDMRG::single_xDMRG_step() {
    tools::log->debug("Starting xDMRG step {} | iter {} | pos {} | dir {}", sim_status.step, sim_status.iter, sim_status.position, sim_status.direction);

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
    if(state_number == 0 and state->get_chi_lim() <= 12) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    if(state_number == 0 and (state->get_chi_lim() <= 8 or sim_status.iter < 2)) {
        optMode  = OptMode::OVERLAP;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    // Setup strong overrides to normal conditions, e.g., for experiments like chi quench

    if(state->size_2site() < settings::precision::max_size_part_diag and chi_quench_steps > 0) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    if(state->size_2site() < settings::precision::max_size_part_diag and tools::finite::measure::energy_variance(*state) > 1e-4) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::SUBSPACE_AND_DIRECT;
    }

    //    if(state->size_2site() < settings::precision::max_size_part_diag and tools::finite::measure::energy_variance(*state) > 1e-2){
    //        optMode  = OptMode::OVERLAP;
    //        optSpace = OptSpace::SUBSPACE_ONLY;
    //    }

    if(sim_status.variance_mpo_has_converged) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }

    if(state->isReal()) optType = OptType::REAL;

    size_t threshold = 0;
    switch(optSpace) {
        case OptSpace::DIRECT: threshold = settings::precision::max_size_direct; break;
        case OptSpace::SUBSPACE_ONLY: threshold = settings::precision::max_size_part_diag; break;
        case OptSpace::SUBSPACE_AND_DIRECT: threshold = settings::precision::max_size_part_diag; break;
    }

    std::list<size_t> max_num_sites_list;
    // Generate a list of maximum number of active sites to try
    if(sim_status.simulation_has_stuck_for > 0)
        max_num_sites_list = {4};
    else if(sim_status.simulation_has_stuck_for > 1)
        max_num_sites_list = {settings::precision::max_sites_multidmrg};
    else if(chi_quench_steps > 0)
        max_num_sites_list = {settings::precision::max_sites_multidmrg};
    else if(optSpace == OptSpace::SUBSPACE_AND_DIRECT)
        max_num_sites_list = {2};
    else if(optSpace == OptSpace::SUBSPACE_ONLY)
        max_num_sites_list = {2};
    else
        max_num_sites_list = {2};

    // Make sure not to use SUBSPACE if the 2-site problem size is huge
    // When this happens, we can't use SUBSPACE even with two sites,
    // so we might as well give it to DIRECT, which handles larger problems.
    if(state->size_2site() > threshold) {
        optMode  = OptMode::VARIANCE;
        optSpace = OptSpace::DIRECT;
    }

    max_num_sites_list.sort();
    max_num_sites_list.unique();
    max_num_sites_list.remove_if([](auto &elem) { return elem > settings::precision::max_sites_multidmrg; });
    if(max_num_sites_list.empty()) throw std::runtime_error("No sites selected for multisite xDMRG");

    tools::log->debug("Possible multisite step sizes: {}", max_num_sites_list);
    size_t                                                                                    theta_count  = 0;
    double                                                                                    variance_old = measure::energy_variance_per_site(*state);
    std::map<double, std::tuple<Eigen::Tensor<Scalar, 3>, std::list<size_t>, size_t, double>> results;
    for(auto &max_num_sites : max_num_sites_list) {
        if(optMode == opt::OptMode::OVERLAP and optSpace == opt::OptSpace::DIRECT) throw std::logic_error("OVERLAP mode and DIRECT space are incompatible");

        auto old_num_sites = state->active_sites.size();
        auto old_prob_size = state->active_problem_size();
        state->activate_sites(threshold, max_num_sites);
        //        if(chi_quench_steps == 0) state->activate_sites(threshold, max_num_sites);
        //        else state->activate_truncated_sites(threshold,chi_lim_quench_ahead, max_num_sites);
        //         Reduce bond dimensions for some sites ahead
        //        if(chi_quench_steps > 0) tools::finite::mps::truncate_active_sites(*state, chi_lim_quench_ahead);

        // Check that we are not about to solve the same problem again
        if(not results.empty() and state->active_sites.size() == old_num_sites and state->active_problem_size() == old_prob_size) {
            // If we reached this point we have exhausted the number of sites available
            tools::log->debug("Can't activate more sites");
            break;
        }

        auto   theta        = opt::find_excited_state(*state, sim_status, optMode, optSpace, optType);
        double variance_new = measure::energy_variance_per_site(*state, theta);
        results.insert({variance_new, {theta, state->active_sites, theta_count++, tools::common::profile::t_opt->get_last_time_interval()}});
        //
        //        if(state_is_within_energy_window(theta)) {
        //            variance_new = measure::energy_variance_per_site(*state, theta);
        //            results.insert({variance_new, {theta, state->active_sites, theta_count++}});
        //        }else{
        //            tools::log->info("Rejecting state found out of energy window");
        //            theta = state->get_multisite_mps();
        //            variance_new = measure::energy_variance_per_site(*state);
        //            results.insert({variance_new, {theta, state->active_sites, theta_count++}});
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

    // Check the contents of results.
    //    auto[variance_new,theta] = std::make_pair(results.begin()->first,results.begin()->second);
    state->clear_cache();
    state->clear_measurements();
    auto        variance_new = results.begin()->first;
    const auto &theta        = std::get<0>(results.begin()->second);
    state->active_sites      = std::get<1>(results.begin()->second);

    if(std::log10(variance_new) < std::log10(variance_old) - 1e-2) state->tag_active_sites_have_been_updated(true);

    // Truncate theta down to chi_lim
    auto chi_lim = state->get_chi_lim();

    // Truncate even more if doing chi quench
    //    if(chi_quench_steps > 0) chi_lim = chi_lim_quench_trail;

    // Do the truncation with SVD
    auto variance_before_svd = tools::finite::measure::energy_variance_per_site(*state, theta);
    tools::log->trace("Variance check before SVD: {:.16f}", std::log10(variance_before_svd));
    opt::truncate_theta(theta, *state, chi_lim);
    auto variance_after_svd = tools::finite::measure::energy_variance_per_site(*state);
    state->set_truncated_variance((variance_after_svd - variance_before_svd) / variance_after_svd);
    tools::log->trace("Variance check after  SVD: {:.16f}", std::log10(variance_after_svd));
    tools::log->debug("Variance loss due to  SVD: {:.16f}", state->get_truncated_variance());

    // Normalize if unity was lost for some reason (numerical error buildup)
    if(std::abs(tools::finite::measure::norm(*state) - 1.0) > settings::precision::max_norm_error) {
        tools::log->warn("Norm too large: {:.18f}", tools::finite::measure::norm(*state));
        tools::finite::mps::normalize(*state);
        tools::finite::mps::rebuild_environments(*state);
    }
    if(settings::precision::use_reduced_energy and state->position_is_any_edge()) {
        tools::finite::mpo::reduce_mpo_energy(*model,*state);
    }

    debug::check_integrity(*state);
    sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min) / (sim_status.energy_max - sim_status.energy_min);

    tools::common::profile::t_sim->toc();
    sim_status.wall_time = tools::common::profile::t_tot->get_age();
    sim_status.simu_time = tools::common::profile::t_sim->get_measured_time();
}

void class_xDMRG::check_convergence() {
    tools::common::profile::t_con->tic();
    if(state->position_is_any_edge()) {
        check_convergence_variance();
        check_convergence_entg_entropy();
    }

    //    sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min ) / (sim_status.energy_max - sim_status.energy_min);
    bool outside_of_window = std::abs(sim_status.energy_dens - sim_status.energy_dens_target) > sim_status.energy_dens_window;
    if(sim_status.iter > 2 and state->position_is_any_edge()) {
        if(outside_of_window and (sim_status.variance_mpo_has_saturated or sim_status.variance_mpo_has_converged or
                                  tools::finite::measure::energy_variance_per_site(*state) < 1e-4)) {
            double      old_energy_dens_window = sim_status.energy_dens_window;
            double      new_energy_dens_window = std::min(energy_window_growth_factor * sim_status.energy_dens_window, 0.5);
            std::string reason = fmt::format("energy {:.16f} saturated outside of energy window {} ± {}", tools::finite::measure::energy_per_site(*state),
                                             sim_status.energy_dens_target, sim_status.energy_dens_window);
            tools::log->info("Increasing energy window: {} --> {}", old_energy_dens_window, new_energy_dens_window);
            sim_status.energy_dens_window = new_energy_dens_window;
            reset_to_random_state_in_energy_window(settings::strategy::initial_parity_sector, false, reason);
        }
        //        else
        //        if( not     state->all_sites_updated()
        //            and     sim_status.simulation_has_got_stuck
        //            and     tools::finite::measure::energy_variance_per_site(*state) > 1e-4)
        //        {
        //            sim_status.energy_dens_window = std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5);
        //            std::string reason = fmt::format("could not update all sites. Energy density: {}, Energy window: {} --> {}",
        //                     sim_status.energy_dens, sim_status.energy_dens_window, std::min(energy_window_growth_factor*sim_status.energy_dens_window, 0.5)
        //                     );
        //            reset_to_random_state_in_energy_window(settings::strategy::initial_parity_sector, false, reason);
        //        }
    }

    sim_status.simulation_has_converged = sim_status.variance_mpo_has_converged and sim_status.entanglement_has_converged;

    sim_status.simulation_has_saturated =
        ((sim_status.variance_mpo_saturated_for >= min_saturation_iters and sim_status.entanglement_saturated_for >= min_saturation_iters) or
         (state->get_iteration() > settings::xdmrg::min_iters and not state->any_sites_updated()));

    sim_status.simulation_has_succeeded = sim_status.simulation_has_converged and sim_status.simulation_has_saturated;

    sim_status.simulation_has_got_stuck = sim_status.simulation_has_saturated and not sim_status.simulation_has_succeeded;

    if(state->position_is_any_edge()) {
        sim_status.simulation_has_stuck_for = sim_status.simulation_has_got_stuck ? sim_status.simulation_has_stuck_for + 1 : 0;
    }

    sim_status.simulation_has_to_stop = sim_status.simulation_has_stuck_for >= max_stuck_iters;

    if(state->position_is_any_edge()) {
        tools::log->debug("Simulation has converged: {}", sim_status.simulation_has_converged);
        tools::log->debug("Simulation has saturated: {}", sim_status.simulation_has_saturated);
        tools::log->debug("Simulation has succeeded: {}", sim_status.simulation_has_succeeded);
        tools::log->debug("Simulation has got stuck: {}", sim_status.simulation_has_got_stuck);
        tools::log->debug("Simulation has stuck for: {}", sim_status.simulation_has_stuck_for);
        tools::log->debug("Simulation has to stop  : {}", sim_status.simulation_has_to_stop);
    }

    //    if (    sim_status.num_resets < settings::precision::max_resets
    //            and tools::finite::measure::energy_variance_per_site(*state) > 1e-10)
    //    {
    //        std::string reason = fmt::format("simulation has saturated with bad precision",
    //                                         sim_status.energy_dens, sim_status.energy_dens_window, sim_status.energy_dens_window);
    //        reset_to_random_state_in_energy_window(settings::strategy::initial_parity_sector, false, reason);
    //    }

    tools::common::profile::t_con->toc();
}

void class_xDMRG::inflate_initial_state() {
    tools::log->trace("Inflating bond dimension");
    // Inflate by projecting randomly. Each projection doubles the bond dimension
    for(int i = 0; i < 4; i++) {
        *state = tools::finite::ops::get_projection_to_closest_parity_sector(*state, "random");
        tools::log->debug("χ = {}", tools::finite::measure::bond_dimensions(*state));
    }
    *state = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::strategy::initial_parity_sector);
}

void class_xDMRG::reset_to_random_state_in_energy_window(const std::string &parity_sector, bool inflate, std::string reason) {
    tools::log->info("Resetting to product state -- Reason: {}", reason);
    tools::log->info("Searching for product state in normalized energy range: {} +- {} and sector {}", sim_status.energy_dens_target,
                     sim_status.energy_dens_window, parity_sector);

    sim_status.num_resets++;
    if(sim_status.num_resets > settings::precision::max_resets) {
        tools::log->info("Not allowed more resets: num resets {} > max resets {}", sim_status.num_resets, settings::precision::max_resets);
        return;
    }

    int  counter           = 0;
    bool outside_of_window = true;

    while(true) {
        reset_to_random_product_state(parity_sector);
        if(inflate) inflate_initial_state();
        sim_status.energy_dens = (tools::finite::measure::energy_per_site(*state) - sim_status.energy_min) / (sim_status.energy_max - sim_status.energy_min);
        outside_of_window      = std::abs(sim_status.energy_dens - sim_status.energy_dens_target) >= sim_status.energy_dens_window;
        tools::log->info("New energy density: {:.16f} | window {} | outside of window: {}", sim_status.energy_dens, sim_status.energy_dens_window,
                         outside_of_window);
        if(not outside_of_window) break;
        counter++;
        if(counter >= 200) throw std::runtime_error(fmt::format("Failed to find initial state in energy window after {}. retries: ", counter));
        if(counter % 10 == 0 and energy_window_growth_factor != 1.0) {
            double old_energy_dens_window = sim_status.energy_dens_window;
            double new_energy_dens_window = std::min(energy_window_growth_factor * sim_status.energy_dens_window, 0.5);

            tools::log->info("Can't find state in energy window.  Increasing energy window: {} --> {}", old_energy_dens_window, new_energy_dens_window);
            sim_status.energy_dens_window = new_energy_dens_window;
        }
    }
    tools::log->info("Energy initial (per site) = {:.16f} | density = {:.8f} | retries = {}", tools::finite::measure::energy_per_site(*state),
                     sim_status.energy_dens, counter);
    clear_saturation_status();
    sim_status.energy_ubound = sim_status.energy_target + sim_status.energy_dens_window * (sim_status.energy_max - sim_status.energy_min);
    sim_status.energy_lbound = sim_status.energy_target - sim_status.energy_dens_window * (sim_status.energy_max - sim_status.energy_min);

    tools::log->info("Number of product state resets: {}", sim_status.num_resets);
}

bool class_xDMRG::state_is_within_energy_window(Eigen::Tensor<Scalar, 3> &theta) {
    double energy_std = std::sqrt(tools::finite::measure::energy_variance(*state));
    double energy_old = tools::finite::measure::energy(*state);
    double energy_new = tools::finite::measure::energy(*state, theta);
    double energy_dif = std::abs(energy_old - energy_new);
    // We can accept a new candidate state if it is less than +-2 standard deviation away from the current energy
    if(energy_dif < 1 * energy_std) {
        sim_status.energy_target = energy_new / static_cast<double>(state->get_length());
        sim_status.energy_ubound = std::min(sim_status.energy_ubound, sim_status.energy_target + 1 * energy_std / static_cast<double>(state->get_length()));
        sim_status.energy_lbound = std::max(sim_status.energy_lbound, sim_status.energy_target - 1 * energy_std / static_cast<double>(state->get_length()));
        sim_status.energy_ubound = std::max(sim_status.energy_ubound, sim_status.energy_lbound);
        sim_status.energy_lbound = std::min(sim_status.energy_ubound, sim_status.energy_lbound);
        sim_status.energy_dens   = (sim_status.energy_target - sim_status.energy_min) / (sim_status.energy_max - sim_status.energy_min);
        sim_status.energy_dens_target = (sim_status.energy_target - sim_status.energy_min) / (sim_status.energy_max - sim_status.energy_min);
        sim_status.energy_dens_window = (sim_status.energy_ubound - sim_status.energy_lbound) / (sim_status.energy_max - sim_status.energy_min);
        tools::log->info("Energy maximum (per site) = {:.16f}", sim_status.energy_max);
        tools::log->info("Energy minimum (per site) = {:.16f}", sim_status.energy_min);
        tools::log->info("Energy target  (per site) = {:.16f}", sim_status.energy_target);
        tools::log->info("Energy ubound  (per site) = {:.16f}", sim_status.energy_ubound);
        tools::log->info("Energy lbound  (per site) = {:.16f}", sim_status.energy_lbound);
        return true;
    } else {
        return false;
    }
}

void class_xDMRG::find_energy_range() {
    tools::log->trace("Finding energy range");
    class_fDMRG fDMRG(h5pp_file);
    fDMRG.ritz = StateRitz::SR;
    fDMRG.run_simulation();
    fDMRG.state->do_all_measurements();

    exit(0);

    if(state->get_length() != settings::model::model_size) throw std::runtime_error("find_energy_range: state length mismatch");



    size_t max_sweeps_during_f_range = 4;
    // Backup the state and current sim_status
    class_simulation_status sim_status_backup = sim_status;
    class_state_finite      state_backup      = *state;

    update_bond_dimension_limit(16);

    // Find energy minimum
    tools::finite::mps::random_product_state(*state, "random", -1, true);
    sim_status.iter = state->reset_iter();
    sim_status.step = state->reset_step();
    while(true) {
//        class_algorithm_finite::single_fDMRG_step();
        print_status_update();
        // It's important not to perform the last moves. That last state would not get optimized
        if(state->position_is_any_edge())
            if(sim_status.iter >= max_sweeps_during_f_range or tools::finite::measure::energy_variance_per_site(*state) < 1e-8) break;
        move_center_point();
    }
    double energy_min = tools::finite::measure::energy_per_site(*state);
    tools::finite::mps::normalize(*state);
    write_to_file(StorageReason::EMIN_STATE);

    // Find energy maximum
    tools::finite::mps::random_product_state(*state, "random", -1, true);
    sim_status.iter = state->reset_iter();
    sim_status.step = state->reset_step();
    while(true) {
//        class_algorithm_finite::single_fDMRG_step();
        print_status_update();
        // It's important not to perform the last moves. That last state would not get optimized
        if(state->position_is_any_edge())
            if(sim_status.iter >= max_sweeps_during_f_range or tools::finite::measure::energy_variance_per_site(*state) < 1e-8) break;
        move_center_point();
    }
    double energy_max = tools::finite::measure::energy_per_site(*state);
    tools::finite::mps::normalize(*state);
    write_to_file(StorageReason::EMAX_STATE);

    // Recover the backup
    sim_status = sim_status_backup;
    *state     = state_backup;

    // Define energy targets and bounds
    sim_status.energy_min    = energy_min;
    sim_status.energy_max    = energy_max;
    sim_status.energy_target = sim_status.energy_min + sim_status.energy_dens_target * (sim_status.energy_max - sim_status.energy_min);
    sim_status.energy_ubound = sim_status.energy_target + sim_status.energy_dens_window * (sim_status.energy_max - sim_status.energy_min);
    sim_status.energy_lbound = sim_status.energy_target - sim_status.energy_dens_window * (sim_status.energy_max - sim_status.energy_min);
    tools::log->info("Energy minimum (per site) = {:.8f}", sim_status.energy_min);
    tools::log->info("Energy maximum (per site) = {:.8f}", sim_status.energy_max);
    tools::log->info("Energy target  (per site) = {:.8f}", sim_status.energy_target);
    tools::log->info("Energy lbound  (per site) = {:.8f}", sim_status.energy_lbound);
    tools::log->info("Energy ubound  (per site) = {:.8f}", sim_status.energy_ubound);
}

bool   class_xDMRG::sim_on() { return settings::xdmrg::on; }
long   class_xDMRG::chi_max() { return settings::xdmrg::chi_max; }
//size_t class_xDMRG::write_freq() { return settings::xdmrg::write_freq; }
size_t class_xDMRG::print_freq() { return settings::xdmrg::print_freq; }
bool   class_xDMRG::chi_grow() { return settings::xdmrg::chi_grow; }
long   class_xDMRG::chi_init() { return settings::xdmrg::chi_init; }
bool   class_xDMRG::store_wave_function() { return settings::xdmrg::store_wavefn; }
