//
// Created by david on 2018-01-31.
//

#include "class_fDMRG.h"
#include <simulation/nmspc_settings.h>
#include <state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/ops.h>

using namespace std;
using namespace Textra;

class_fDMRG::class_fDMRG(std::shared_ptr<h5pp::File> h5ppFile_) : class_algorithm_finite(std::move(h5ppFile_), "fDMRG", SimulationType::fDMRG) {
    tools::log->trace("Constructing class_fDMRG");
    settings::fdmrg::min_sweeps = std::max(settings::fdmrg::min_sweeps, (size_t)(std::log2(chi_max())));
}

void class_fDMRG::run_simulation() {
    tools::log->info("Starting {} simulation", sim_name);
    while(true) {
        single_DMRG_step(ritz);
        print_status_update();
        write_to_file();
        copy_from_tmp();
        check_convergence();
        update_truncation_limit();     // Will update SVD threshold iff the state precision is being limited by truncation error
        update_bond_dimension_limit(); // Will update bond dimension if the state precision is being limited by bond dimension

        try_projection();

        // It's important not to perform the last move.
        // That last state would not get optimized
        if(state->position_is_any_edge()) {
            if(sim_status.iter >= settings::xdmrg::max_sweeps) {
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
            if(settings::strategy::randomize_early and sim_status.state_number == 0 and state->find_largest_chi() >= 32 and
               tools::finite::measure::energy_variance(*state) < 1e-4) {
                stop_reason = StopReason::RANDOMIZE;
                break;
            }
        }
        tools::log->trace("Finished step {}, iter {}, pos {}, dir {}", sim_status.step, sim_status.iter, state->get_position(), state->get_direction());
        move_center_point();
    }
    tools::log->info("Finished {} simulation -- reason: {}", sim_name, enum2str(stop_reason));
}

void class_fDMRG::check_convergence() {
    tools::common::profile::t_con->tic();

    if(state->position_is_any_edge()) {
        check_convergence_variance();
        check_convergence_entg_entropy();
    }

    if(sim_status.iter <= settings::fdmrg::min_sweeps) {
        clear_saturation_status();
    }

    if(sim_status.variance_mpo_has_converged and sim_status.entanglement_has_converged
       //           and sim_status.variance_mpo_saturated_for >= min_saturation_iters
       //           and sim_status.entanglement_saturated_for >= min_saturation_iters
    ) {
        tools::log->debug("Simulation has converged");
        sim_status.simulation_has_converged = true;
    }

    if(sim_status.chi_lim_has_reached_chi_max and
       (sim_status.variance_mpo_saturated_for >= max_saturation_iters or sim_status.entanglement_saturated_for >= max_saturation_iters)) {
        tools::log->debug("Simulation has to stop");
        sim_status.simulation_has_to_stop = true;
    }

    if(state->position_is_any_edge() and sim_status.variance_mpo_has_saturated and not sim_status.variance_mpo_has_converged and
       not sim_status.simulation_has_converged and not has_projected) {
        tools::log->info("Projecting to {} due to saturation", settings::strategy::target_parity_sector);
        *state        = tools::finite::ops::get_projection_to_closest_parity_sector(*state, settings::strategy::target_parity_sector);
        has_projected = true;
    }

    tools::common::profile::t_con->toc();
}

bool   class_fDMRG::sim_on() { return settings::fdmrg::on; }
long   class_fDMRG::chi_max() { return settings::fdmrg::chi_max; }
size_t class_fDMRG::write_freq() { return settings::fdmrg::write_freq; }
size_t class_fDMRG::print_freq() { return settings::fdmrg::print_freq; }
bool   class_fDMRG::chi_grow() { return settings::fdmrg::chi_grow; }
long   class_fDMRG::chi_init() { return settings::fdmrg::chi_init; }
bool   class_fDMRG::store_wave_function() { return settings::fdmrg::store_wavefn; }
