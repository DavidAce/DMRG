#include "AlgorithmFinite.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "debug/info.h"
#include "general/iter.h"
#include "math/num.h"
#include "math/tenx/span.h"
#include "qm/spin.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/h5/storage_info.h"
#include "tools/common/log.h"
#include "tools/finite/env.h"
#include "tools/finite/h5.h"
#include "tools/finite/measure.h"
#include "tools/finite/ops.h"
#include "tools/finite/print.h"
#include <h5pp/h5pp.h>

AlgorithmFinite::AlgorithmFinite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type) : AlgorithmBase(std::move(h5ppFile_), algo_type) {
    tools::log->trace("Constructing class_algorithm_finite");
    tensors.initialize(algo_type, settings::model::model_type, settings::model::model_size, 0);
}

// We need to make a destructor manually for the enclosing class "ModelFinite"
// that encloses "class_model_base". Otherwise, unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_model_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
// class_algorithm_finite::~class_algorithm_finite() = default;

void AlgorithmFinite::run()
/*!
 * \brief Dispatches finite DMRG stages.
 * This function manages the stages of simulation differently depending on whether
 * the data already existed in the hdf5 output file or not.
 *
 * We start by searching the HDF5 group tree looking for
 *
 * These are the scenarios:
 * 1) The hdf5 file existed already and contains
 *      a) nothing recognizable (previous crash?)        -- delete the file -> run full simulation from scratch.
 *      b) a converged simulation but no MPS             -- run full simulation from scratch.
 *      c) a not-yet-converged MPS                       -- resume simulation, reset the number of sweeps first.
 *      d) a converged MPS                               -- not much to do... run postprocessing
 * 2) The hdf5 file did not exist                        -- run full simulation from scratch.

        D1: output file exists: check option in settings::storage::create_mode:
            - case OPEN:        Load config and simulation state -> continue simulation
            - case RENAME:      Rename output file -> go to D2
            - case TRUNCATE:    Delete output file -> go to D2

        D2: output file does not exist: Check if coming from A2
            - Yes: load simulation state from given hdf5 file -> continue simulation
            - No: start new simulation

 */
{
    tools::log->info("Starting {}", status.algo_type_sv());
    auto t_tot  = tid::get_unscoped("t_tot").tic_token();
    auto t_algo = tid::tic_scope(status.algo_type_sv());
    // We may want to resume this simulation.
    auto finished_exists = h5file->linkExists("common/finished_all");
    auto algo_exists     = h5file->linkExists(status.algo_type_sv());
    auto policy_is_resume =
        settings::storage::file_collision_policy == FileCollisionPolicy::RESUME or settings::storage::file_collision_policy == FileCollisionPolicy::REVIVE;
    tools::log->debug("common/finished_all exists: {} | should resume: {}", finished_exists, policy_is_resume);
    if(finished_exists and algo_exists and policy_is_resume) {
        try {
            tools::log->info("Attempting resume");
            resume();
        } catch(const except::state_error &ex) {
            throw except::resume_error("Failed to resume state from file [{}]: {}", h5file->getFilePath(), ex.what());
        } catch(const except::load_error &ex) {
            throw except::resume_error("Failed to load simulation from file [{}]: {}", h5file->getFilePath(), ex.what());
        } catch(const std::exception &ex) { throw except::runtime_error("Failed to resume from file [{}]: {}", h5file->getFilePath(), ex.what()); }
    } else {
        run_default_task_list();
    }
}

void AlgorithmFinite::run_postprocessing() {
    tools::log->info("Running default postprocessing for {}", status.algo_type_sv());
    auto tic = tid::tic_scope("post");
    if(settings::strategy::project_final_state) {
        tensors.project_to_nearest_axis(settings::strategy::target_axis, svd::config(status.bond_lim, status.trnc_lim));
        tensors.rebuild_edges();
    }
    write_to_file(StorageEvent::BOND_INCREASE, CopyPolicy::OFF); // To get checkpoint/chi_# with the current result (which would otherwise be missing
    write_to_file(StorageEvent::PROJ_STATE, CopyPolicy::OFF);    // To compare the finished state to a projected one
    write_to_file(StorageEvent::LAST_STATE, CopyPolicy::FORCE);  // For final mps
    print_status_full();
    run_fes_analysis();
    tools::log->info("Finished default postprocessing for {}", status.algo_type_sv());
}

void AlgorithmFinite::move_center_point(std::optional<long> num_moves) {
    if(not num_moves.has_value()) {
        if(tensors.active_sites.empty())
            num_moves = 1;
        else {
            long posL_edge     = 0;
            long posR_edge     = tensors.get_length<long>() - 1;
            long posL_active   = static_cast<long>(tensors.active_sites.front());
            long posR_active   = static_cast<long>(tensors.active_sites.back());
            long num_active    = static_cast<long>(tensors.active_sites.size());
            bool reached_edgeR = tensors.state->get_direction() == 1 and posR_edge == posR_active;
            bool reached_edgeL = tensors.state->get_direction() == -1 and posL_edge == posL_active;

            if(reached_edgeL or reached_edgeR) {
                // In this case we have just updated from here to the edge. No point in updating
                // closer and closer to the edge. Just move until reaching the edge without flip,
                // then once more to flip, then once more to move the center position back from negative.
                //
                // Reminder: moving to the outward edge means reaching pos == -1 and dir == -1, so all sites are "B".
                // Moving again only flips the sign of dir. We move once more to get pos == 0, dir == 1.
                // On the other edge this is not necessary!
                num_moves = std::max<long>(1, num_active - 1) + 1;   // to the edge without flipping, +1 to flip
                if(reached_edgeL) num_moves = num_moves.value() + 1; // , +1 to get pos == 0
            } else if(settings::strategy::multisite_mps_move == MultisiteMove::ONE)
                num_moves = 1ul;
            else if(settings::strategy::multisite_mps_move == MultisiteMove::MID) {
                num_moves = std::max<long>(1, num_active / 2);
            } else if(settings::strategy::multisite_mps_move == MultisiteMove::MAX) {
                num_moves = std::max<long>(1, num_active - 1); // Move so that the center point moves out of the active region
            } else
                throw except::logic_error("Could not determine how many sites to move");
        }
    }

    if(num_moves <= 0) throw except::runtime_error("Cannot move center point {} sites", num_moves.value());
    tools::log->trace("Moving center point {} steps in direction {}", num_moves.value(), tensors.state->get_direction());
    tensors.clear_cache();
    tensors.clear_measurements();
    try {
        long moves = 0;
        while(num_moves > moves++) {
            if(tensors.position_is_outward_edge()) status.iter++;
            tools::log->trace("Moving center position | step {} | pos {} | dir {} ", status.step, tensors.get_position<long>(), tensors.state->get_direction());
            status.step += tensors.move_center_point();
            // Do not go past the edge if you aren't there already!
            // It's important to stay at the inward edge, so we can do convergence checks and so on
            if(tensors.position_is_inward_edge()) break;
        }
        tensors.clear_active_sites();
    } catch(std::exception &e) {
        tools::finite::print::dimensions(tensors);
        throw except::runtime_error("Failed to move center point: {}", e.what());
    }
    status.position  = tensors.state->get_position<long>();
    status.direction = tensors.state->get_direction();
}

void AlgorithmFinite::shift_mpo_energy() {
    if(not settings::precision::use_mpo_energy_shift) return;
    if(not tensors.position_is_inward_edge()) return;
    tensors.shift_mpo_energy();    // Avoid catastrophic cancellation by shifting energy on each mpo by E/L
    tensors.rebuild_mpo_squared(); // The shift clears our squared mpo's. So we have to rebuild them. Compression is retained.
    tensors.rebuild_edges();       // The shift modified all our mpo's. So we have to rebuild all the edges.
    if constexpr(settings::debug) tensors.assert_validity();
}

void AlgorithmFinite::try_moving_sites() {
    if(not settings::strategy::move_sites_when_stuck) return;
    if(not tensors.position_is_inward_edge()) return;
    if(status.algorithm_has_stuck_for == 0) return;
    // Definitions:
    //  pos : the index of a slot on the lattice, which can not be moved. (long)
    //  site: the index of a particle on the lattice, which can be moved. (size_t)
    //  dir : Direction on which to move the position: +1l = right, -1l = left.
    auto prefer_eigs_backup                    = settings::solver::prefer_eigs_over_bfgs;
    auto eigs_iter_max_backup                  = settings::solver::eigs_iter_max;
    auto eigs_tol_min_backup                   = settings::solver::eigs_tol_min;
    auto eigs_ncv_backup                       = settings::solver::eigs_ncv;
    auto multisite_mps_when_backup             = settings::strategy::multisite_mps_when;
    auto multisite_mps_site_max_backup         = settings::strategy::multisite_mps_site_max;
    auto multisite_mps_site_def_backup         = settings::strategy::multisite_mps_site_def;
    settings::solver::eigs_iter_max            = 100;
    settings::solver::eigs_tol_min             = std::min(1e-14, settings::solver::eigs_tol_min);
    settings::solver::eigs_ncv                 = std::max(35ul, settings::solver::eigs_ncv);
    settings::solver::prefer_eigs_over_bfgs    = OptEigs::ALWAYS;
    settings::strategy::multisite_mps_when     = MultisiteWhen::ALWAYS;
    settings::strategy::multisite_mps_site_max = 2;
    settings::strategy::multisite_mps_site_def = 2;

    auto len      = tensors.get_length<long>();
    auto pos      = tensors.get_position<long>();
    auto dir      = pos > len / 2 ? -1l : 1l;
    auto site_seq = dir > 0 ? num::range<long>(0, len - 1, 1) : num::range<long>(1, len, -1);

    tensors.activate_sites({tensors.get_position()});
    tensors.rebuild_edges();
    tools::log->info("Trying to move sites | pos {} | dir {}", pos, dir);
    std::vector<std::string> report;
    auto                     ene_old = tools::finite::measure::energy(tensors);
    auto                     var_old = tools::finite::measure::energy_variance(tensors);
    for(const auto &site : site_seq) {
        std::vector<long> tgt_pos_seq = dir > 0 ? num::range<long>(site + 1, len - 1, 1) : num::range<long>(0, site, -1);
        std::vector<long> tgt_pos_req = dir > 0 ? num::range<long>(site, len - 2, -1) : num::range<long>(1, site + 1, 1);
        tgt_pos_seq.insert(tgt_pos_seq.end(), tgt_pos_req.begin(), tgt_pos_req.end());
        tools::log->debug("Moving site {} dir {} | seq: {}", site, dir, tgt_pos_seq);
        for(const auto &tgt_pos : tgt_pos_seq) {
            update_state();
            print_status();
            if(sites_mps)
                report.emplace_back(fmt::format("sites [{:2}, {:2}] @ pos [{:2}, {:2}] | variance {:.5e}", sites_mps->at(tensors.active_sites.front()),
                                                sites_mps->at(tensors.active_sites.back()), tensors.active_sites.front(), tensors.active_sites.back(),
                                                tools::finite::measure::energy_variance(tensors)));
            status.step += 1;
            tensors.move_site_to_pos(static_cast<size_t>(site), tgt_pos, sites_mps, sites_mpo, tgt_pos);
            tools::log->debug("Labels    : {}", tensors.state->get_labels());
            tools::log->debug("Sites mps : {}", sites_mps.value());
            tools::log->debug("Sites mpo : {}", sites_mpo.value());
        }
        tools::log->info("Resetting MPO's");
        tensors.rebuild_mpo();
        tensors.rebuild_mpo_squared();
        tensors.move_center_point_to_inward_edge();
        tensors.activate_sites();
        tensors.rebuild_edges();
        sites_mpo = std::nullopt;
        sites_mps = std::nullopt;
        status.iter += 1;
        check_convergence();
        if(status.variance_mpo_converged_for > 0) break;
    }
    clear_convergence_status();
    tools::log->info("Finished moving sites");
    for(const auto &r : report) tools::log->info("{}", r);
    tools::log->info("Energy    {:.16f} --> {:.16f}", ene_old, tools::finite::measure::energy(tensors));
    tools::log->info("Variance  {:9.3e} --> {:9.3e}", var_old, tools::finite::measure::energy_variance(tensors));

    if(not tensors.position_is_inward_edge())
        throw except::logic_error("Position {} and direction {} is not an inward edge", tensors.get_position(), tensors.state->get_direction());
    settings::solver::eigs_iter_max            = eigs_iter_max_backup;
    settings::solver::eigs_tol_min             = eigs_tol_min_backup;
    settings::solver::eigs_ncv                 = eigs_ncv_backup;
    settings::solver::prefer_eigs_over_bfgs    = prefer_eigs_backup;
    settings::strategy::multisite_mps_when     = multisite_mps_when_backup;
    settings::strategy::multisite_mps_site_max = multisite_mps_site_max_backup;
    settings::strategy::multisite_mps_site_def = multisite_mps_site_def_backup;
}

void AlgorithmFinite::update_variance_max_digits(std::optional<double> energy) {
    if(not tensors.position_is_inward_edge()) return;
    if(tensors.active_sites.empty()) return;
    if(not energy) energy = tools::finite::measure::energy(tensors);
    double energy_abs                 = std::abs(energy.value());
    double energy_pow                 = energy_abs * energy_abs;
    double digits10                   = std::numeric_limits<double>::digits10;
    double energy_top                 = settings::precision::use_mpo_energy_shift ? energy_abs : energy_pow;
    double energy_exp                 = std::ceil(std::max(0.0, std::log10(energy_top))) + 1;
    double max_digits                 = std::floor(std::max(0.0, digits10 - energy_exp));
    status.energy_variance_max_digits = static_cast<size_t>(max_digits);
    status.energy_variance_prec_limit = std::pow(10.0, -max_digits);
    tools::log->debug("Estimated limit on energy variance precision: {:.3e}", status.energy_variance_prec_limit);
}

void AlgorithmFinite::update_bond_dimension_limit() {
    if(not tensors.position_is_inward_edge()) return;
    status.bond_max                   = settings::get_bond_max(status.algo_type);
    status.bond_limit_has_reached_max = status.bond_lim >= status.bond_max;
    if(settings::strategy::bond_increase_when == UpdateWhen::NEVER) {
        status.bond_lim                   = status.bond_max;
        status.bond_limit_has_reached_max = true;
        return;
    }
    if(status.bond_limit_has_reached_max) return;
    auto tic = tid::tic_scope("bond_grow");

    if constexpr(settings::debug) {
        if(tools::log->level() == spdlog::level::trace) {
            double truncation_threshold = 2 * settings::solver::svd_truncation_lim;
            size_t trunc_bond_count     = tensors.state->num_sites_truncated(truncation_threshold);
            size_t bond_at_lim_count    = tensors.state->num_bonds_at_limit(status.bond_lim);
            tools::log->trace("Truncation threshold  : {:<.8e}", truncation_threshold);
            tools::log->trace("Truncation errors     : {}", tensors.state->get_truncation_errors());
            tools::log->trace("Bond dimensions       : {}", tools::finite::measure::bond_dimensions(*tensors.state));
            tools::log->trace("Truncated bond count  : {} ", trunc_bond_count);
            tools::log->trace("Bonds at limit  count : {} ", bond_at_lim_count);
            tools::log->trace("Entanglement entropies: {} ", tools::finite::measure::entanglement_entropies(*tensors.state));
        }
    }
    // If we got here we want to increase the bond dimension limit progressively during the simulation
    bool is_saturated      = status.algorithm_saturated_for > 0; // Allow one round while saturated so that extra efforts get a chance.
    bool is_has_stuck      = status.algorithm_has_stuck_for > 0;
    bool is_truncated      = tensors.state->is_limited_by_bond(status.bond_lim) or tensors.state->is_truncated(status.trnc_lim);
    bool grow_if_truncated = settings::strategy::bond_increase_when == UpdateWhen::TRUNCATED;
    bool grow_if_saturated = settings::strategy::bond_increase_when == UpdateWhen::SATURATED;
    bool grow_if_has_stuck = settings::strategy::bond_increase_when == UpdateWhen::STUCK;

    if(grow_if_truncated and not is_truncated) {
        tools::log->info("State is not limited by its bond dimension. Kept current bond limit {}", status.bond_lim);
        return;
    }
    if(grow_if_saturated and not is_saturated) {
        tools::log->info("Algorithm is not saturated. Kept current bond limit {}", status.bond_lim);
        return;
    }
    if(grow_if_has_stuck and not is_has_stuck) {
        tools::log->info("Algorithm is not stuck. Kept current bond limit {}", status.bond_lim);
        return;
    }

    // Do a projection to make sure the saved data is in the correct sector
    if(settings::strategy::project_on_bond_update) {
        tensors.project_to_nearest_axis(settings::strategy::target_axis, svd::config(status.bond_lim, status.trnc_lim));
        tensors.rebuild_edges();
    }
    // Write current results before updating bond dimension
    write_to_file(StorageEvent::BOND_INCREASE);
    if(settings::strategy::randomize_on_bond_update and status.bond_lim >= 32)
        randomize_state(ResetReason::BOND_UPDATE, StateInit::RANDOMIZE_PREVIOUS_STATE, std::nullopt, std::nullopt);

    // If we got to this point we will update the bond dimension by a factor
    auto rate = settings::strategy::bond_increase_rate;
    if(rate <= 1.0) throw except::runtime_error("Error: get_bond_grow_rate == {:.3f} | must be larger than one", rate);

    auto bond_new = static_cast<double>(status.bond_lim);
    if(rate <= 2.0 and rate > 1.0) {
        bond_new = std::ceil(bond_new * rate);
        bond_new = num::round_up_to_multiple_of<double>(bond_new, 4);
    } else if(rate > 2.0) {
        bond_new = bond_new + rate;
    } else
        throw except::logic_error("Expected grow_rate > 1.0. Got {}", rate);
    bond_new = std::min(bond_new, static_cast<double>(status.bond_max));

    tools::log->info("Updating bond dimension limit {} -> {} | truncated {} | saturated {}", status.bond_lim, bond_new, is_truncated, is_saturated);
    status.bond_lim                   = static_cast<long>(bond_new);
    status.bond_limit_has_reached_max = status.bond_lim == status.bond_max;
    status.algorithm_has_stuck_for    = 0;
    status.algorithm_saturated_for    = 0;

    // Last sanity check before leaving here
    if(status.bond_lim > status.bond_max) throw except::logic_error("bond_lim is larger than get_bond_max! {} > {}", status.bond_lim, status.bond_max);
}

void AlgorithmFinite::reduce_bond_dimension_limit(double rate, UpdateWhen when, StorageEvent storage_event) {
    // We reduce the bond dimension limit during FES whenever entanglement has stopped changing
    if(not tensors.position_is_inward_edge()) return;
    if(when == UpdateWhen::NEVER) return;
    if(iter_last_bond_reduce == 0) iter_last_bond_reduce = status.iter;
    size_t iter_since_reduce = std::max(status.iter, iter_last_bond_reduce) - std::min(status.iter, iter_last_bond_reduce);
    bool   reduce_saturated  = when == UpdateWhen::SATURATED and status.algorithm_saturated_for > 0;
    bool   reduce_truncated  = when == UpdateWhen::TRUNCATED and tensors.state->is_truncated(status.trnc_lim);
    bool   reduce_iteration  = when == UpdateWhen::ITERATION and iter_since_reduce >= 1;
    if(reduce_saturated or reduce_truncated or reduce_iteration or iter_since_reduce >= 20) {
        write_to_file(storage_event, CopyPolicy::OFF);

        auto bond_new = static_cast<double>(status.bond_lim);
        if(rate > 0.0 and rate < 1.0)
            bond_new *= rate;
        else if(rate >= 1.0)
            bond_new -= rate;
        else
            throw except::logic_error("invalid rate {}", rate);
        bond_new = std::floor(std::max(bond_new, 1.0));
        if(bond_new == static_cast<double>(status.bond_lim))
            status.algo_stop = AlgorithmStop::SUCCESS; // There would be no change in bond_lim
        else {
            if(storage_event != StorageEvent::NONE) tools::log->info("Updating bond dimension limit {} -> {}", status.bond_lim, bond_new);
            status.bond_lim = static_cast<long>(bond_new);
        }
        iter_last_bond_reduce = status.iter;
    }
}

void AlgorithmFinite::update_truncation_error_limit() {
    if(not tensors.position_is_inward_edge()) return;
    if(status.trnc_lim == 0.0) throw std::runtime_error("trnc_lim is zero!");
    status.trnc_min                   = settings::solver::svd_truncation_lim;
    status.trnc_limit_has_reached_min = status.trnc_lim <= status.trnc_min;
    if(settings::strategy::trnc_decrease_when == UpdateWhen::NEVER or settings::strategy::trnc_decrease_rate == 0.0) {
        status.trnc_lim                   = status.trnc_min;
        status.trnc_limit_has_reached_min = true;
        return;
    }
    if(status.trnc_limit_has_reached_min) return;
    auto tic = tid::tic_scope("trnc_down");

    if constexpr(settings::debug) {
        if(tools::log->level() == spdlog::level::trace) {
            double truncation_threshold = 2 * settings::solver::svd_truncation_lim;
            size_t trunc_bond_count     = tensors.state->num_sites_truncated(truncation_threshold);
            tools::log->trace("Truncation threshold  : {:<.8e}", truncation_threshold);
            tools::log->trace("Truncation errors     : {}", tensors.state->get_truncation_errors());
            tools::log->trace("Bond dimensions       : {}", tools::finite::measure::bond_dimensions(*tensors.state));
            tools::log->trace("Truncated bond count  : {} ", trunc_bond_count);
            tools::log->trace("Entanglement entropies: {} ", tools::finite::measure::entanglement_entropies(*tensors.state));
        }
    }

    // If we got here we want to decrease the truncation error limit progressively during the simulation
    bool is_saturated      = status.algorithm_saturated_for > 0; // Allow one round while saturated so that extra efforts get a chance.
    bool is_has_stuck      = status.algorithm_has_stuck_for > 0; // Allow one round while saturated so that extra efforts get a chance.
    bool is_truncated      = tensors.state->is_limited_by_bond(status.bond_lim) or tensors.state->is_truncated(status.trnc_lim);
    bool drop_if_truncated = settings::strategy::trnc_decrease_when == UpdateWhen::TRUNCATED;
    bool drop_if_saturated = settings::strategy::trnc_decrease_when == UpdateWhen::SATURATED;
    bool drop_if_has_stuck = settings::strategy::trnc_decrease_when == UpdateWhen::STUCK;

    if(drop_if_truncated and not is_truncated) {
        tools::log->info("State is not truncated. Kept current truncation error limit {:8.2e}", status.trnc_lim);
        return;
    }
    if(drop_if_saturated and not is_saturated) {
        tools::log->info("Algorithm is not saturated. Kept current truncation error limit {:8.2e}", status.trnc_lim);
        return;
    }
    if(drop_if_has_stuck and not is_has_stuck) {
        tools::log->info("Algorithm is not stuck. Kept current truncation error limit {:8.2e}", status.trnc_lim);
        return;
    }

    // Do a projection to make sure the saved data is in the correct sector
    if(settings::strategy::project_on_bond_update) {
        tensors.project_to_nearest_axis(settings::strategy::target_axis, svd::config(status.bond_lim, status.trnc_lim));
        tensors.rebuild_edges();
    }
    // Write current results before updating the truncation error limit
    write_to_file(StorageEvent::TRNC_DECREASE);

    // If we got to this point we will update the truncation error limit by a factor
    auto rate = settings::strategy::trnc_decrease_rate;
    if(rate > 1.0 or rate < 0) throw except::runtime_error("Error: trnc_decrease_rate == {:8.2e} | must be in [0, 1]");

    auto trnc_new = std::max(status.trnc_min, status.trnc_lim * rate);

    tools::log->info("Updating truncation error limit {:8.2e} -> {:8.2e} | truncated {} | saturated {}", status.trnc_lim, trnc_new, is_truncated, is_saturated);
    status.trnc_lim                   = trnc_new;
    status.trnc_limit_has_reached_min = status.trnc_lim == status.trnc_min;

    // Last sanity check before leaving here
    if(status.trnc_lim < status.trnc_min) throw except::logic_error("trnc_lim is smaller than trnc_min ! {:8.2e} > {:8.2e}", status.trnc_lim, status.trnc_min);
}

void AlgorithmFinite::update_expansion_factor_alpha() {
    if(settings::strategy::max_env_expansion_alpha > 0) {
        // Set a good initial value to start with
        if(status.env_expansion_alpha == 0) status.env_expansion_alpha = std::min(status.energy_variance_lowest, settings::strategy::max_env_expansion_alpha);

        // Update alpha
        double old_expansion_alpha = status.env_expansion_alpha;
        double factor_up           = std::pow(1e+1, 1.0 / static_cast<double>(settings::model::model_size)); // A sweep increases alpha by x10
        double factor_dn           = std::pow(1e-3, 1.0 / static_cast<double>(settings::model::model_size)); // A sweep decreases alpha by x1000

        bool   var_has_improved  = status.energy_variance_lowest / status.env_expansion_variance < 0.9;
        bool   var_has_converged = status.variance_mpo_converged_for > 0 or status.energy_variance_lowest < settings::precision::variance_convergence_threshold;
        double alpha_upper_limit = settings::strategy::max_env_expansion_alpha;
        double alpha_lower_limit = std::min(alpha_upper_limit, status.energy_variance_lowest);

        if(status.algorithm_has_stuck_for == 0 or var_has_improved or var_has_converged)
            status.env_expansion_alpha *= factor_dn;
        else
            status.env_expansion_alpha *= factor_up;

        status.env_expansion_alpha = std::clamp(status.env_expansion_alpha, alpha_lower_limit, alpha_upper_limit);
        if(status.env_expansion_alpha < old_expansion_alpha) {
            status.env_expansion_step     = status.step;
            status.env_expansion_variance = status.energy_variance_lowest;
            tools::log->trace("Decreased alpha {:8.2e} -> {:8.2e}", old_expansion_alpha, status.env_expansion_alpha);
        } else if(status.env_expansion_alpha > old_expansion_alpha) {
            tools::log->trace("Increased alpha {:8.2e} -> {:8.2e}", old_expansion_alpha, status.env_expansion_alpha);
        }
    }
}

void AlgorithmFinite::randomize_model() {
    tools::log->info("Randomizing model");
    auto tic = tid::tic_scope("rnd_model");
    tensors.randomize_model();
    clear_convergence_status();
}

void AlgorithmFinite::randomize_state(ResetReason reason, StateInit state_init, std::optional<StateInitType> state_type, std::optional<std::string> sector,
                                      std::optional<bool> use_eigenspinors, std::optional<size_t> bitfield, std::optional<long> bond_lim,
                                      std::optional<double> trnc_lim) {
    auto t_rnd = tid::tic_scope("rnd_state");
    if(reason == ResetReason::SATURATED) {
        if(status.num_resets >= settings::strategy::max_resets)
            return tools::log->warn("Skipped reset: num resets {} >= max resets {}", status.num_resets, settings::strategy::max_resets);
        else
            status.num_resets++; // Only increment if doing it for saturation reasons
    }
    if(not state_type) state_type = tensors.state->is_real() ? StateInitType::REAL : StateInitType::CPLX;
    if(not sector) sector = settings::strategy::initial_axis;
    if(not use_eigenspinors) use_eigenspinors = settings::strategy::use_eigenspinors;
    if(not bitfield) bitfield = settings::input::bitfield;
    if(not bond_lim) {
        bond_lim = settings::get_bond_init(status.algo_type);
        if(settings::strategy::bond_increase_when == UpdateWhen::NEVER and state_init == StateInit::RANDOMIZE_PREVIOUS_STATE)
            bond_lim = static_cast<long>(std::pow(2, std::floor(std::log2(tensors.state->find_largest_bond())))); // Nearest power of two from below
    }
    if(not trnc_lim) {
        trnc_lim = settings::solver::svd_truncation_init;
        if(state_init == StateInit::RANDOMIZE_PREVIOUS_STATE) trnc_lim = 1e-2;
    }

    tensors.activate_sites(settings::solver::max_size_shift_invert, 2); // Activate a pair of sites so that asserts and measurements work
    tensors.rebuild_edges();
    tensors.randomize_state(reason, state_init, state_type.value(), sector.value(), use_eigenspinors.value(), bitfield.value(), bond_lim.value());

    if(settings::strategy::project_initial_state and qm::spin::half::is_valid_axis(sector.value())) {
        tools::log->info("Projecting state | target sector {} | norm {:.16f} | spin components: {:+.16f}", sector.value(),
                         tools::finite::measure::norm(*tensors.state), fmt::join(tools::finite::measure::spin_components(*tensors.state), ", "));
        tensors.project_to_nearest_axis(sector.value(), svd::config(bond_lim, trnc_lim));
        tensors.rebuild_edges();
        // Note! After running this function we should rebuild edges! However, there are usually no sites active at this point, so we do it further down.
    }

    clear_convergence_status();
    status.reset();
    status.iter      = 0;
    status.step      = 0;
    status.position  = tensors.state->get_position<long>();
    status.direction = tensors.state->get_direction();
    status.algo_stop = AlgorithmStop::NONE;
    if(settings::strategy::bond_increase_when != UpdateWhen::NEVER) status.bond_lim = bond_lim.value();
    if(reason == ResetReason::NEW_STATE) excited_state_number++;
    if(tensors.state->find_largest_bond() > bond_lim.value())
        //        tools::log->warn("Faulty truncation after randomize. Max found bond is {}, but bond limit is {}", tensors.state->find_largest_bond(),
        //        bond_lim.value());
        throw except::runtime_error("Faulty truncation after randomize. Max found bond dimension is {}, but bond limit is {}",
                                    tensors.state->find_largest_bond(), bond_lim.value());

    tensors.rebuild_edges();
    tools::log->info("Randomization successful:");
    tools::log->info("-- State labels             : {}", tensors.state->get_labels());
    tools::log->info("-- Normalization            : {:.16f}", tools::finite::measure::norm(*tensors.state));
    tools::log->info("-- Spin components (X,Y,Z)  : {:.16f}", fmt::join(tools::finite::measure::spin_components(*tensors.state), ", "));
    tools::log->info("-- Bond dimensions          : {}", tools::finite::measure::bond_dimensions(*tensors.state));

    if(status.algo_type != AlgorithmType::fLBIT) {
        tools::log->info("-- Energy per site          : {}", tools::finite::measure::energy_per_site(tensors));
        tools::log->info("-- Energy density           : {}", tools::finite::measure::energy_normalized(tensors, status.energy_min, status.energy_max));
        tools::log->info("-- Energy variance          : {:8.2e}", tools::finite::measure::energy_variance(tensors));
    }
}

void AlgorithmFinite::try_projection(std::optional<std::string> target_sector) {
    if(not tensors.position_is_inward_edge()) return;
    if(not target_sector and projected_iter == status.iter) return;
    if(status.variance_mpo_converged_for > 0 and status.spin_parity_has_converged) return; // No need
    size_t iter_since_last_projection = std::max(projected_iter, status.iter) - projected_iter;

    bool project_on_spin_saturation = settings::strategy::project_on_saturation > 0 and not status.spin_parity_has_converged and
                                      iter_since_last_projection >= settings::strategy::project_on_saturation;

    bool project_on_var_saturation = settings::strategy::project_on_saturation > 0 and status.algorithm_saturated_for > 0 and
                                     iter_since_last_projection >= settings::strategy::project_on_saturation;

    bool project_on_every_iter = settings::strategy::project_on_every_iter > 0 and iter_since_last_projection >= settings::strategy::project_on_every_iter;

    bool project_to_given_sector = target_sector.has_value();

    if(project_on_every_iter or project_on_var_saturation or project_to_given_sector or project_on_spin_saturation) {
        if(not target_sector) target_sector = settings::strategy::target_axis;
        if(not qm::spin::half::is_valid_axis(target_sector.value())) return; // Do not project unless the target sector is one of +- xyz
        std::string msg;
        if(project_on_spin_saturation) msg += " | reason: spin component has not converged";
        if(project_on_var_saturation) msg += fmt::format(" | reason: run every {} iter on variance saturation", settings::strategy::project_on_saturation);
        if(project_on_every_iter) msg += fmt::format(" | reason: run every {} iter", settings::strategy::project_on_every_iter);
        tools::log->info("Trying projection to {}{}", target_sector.value(), msg);

        auto sector_sign   = qm::spin::half::get_sign(target_sector.value());
        auto variance_old  = tools::finite::measure::energy_variance(tensors);
        auto spincomp_old  = tools::finite::measure::spin_components(*tensors.state);
        auto entropies_old = tools::finite::measure::entanglement_entropies(*tensors.state);
        if(sector_sign != 0) {
            tensors.project_to_nearest_axis(target_sector.value(), svd::config(status.bond_lim, status.trnc_lim));
            tensors.rebuild_edges();
        } else {
            // We have a choice here.
            // If no sector sign has been given, and the spin component along the requested axis is near zero,
            // then we may inadvertently project to a sector opposite to the target state.
            // If that happened, we would get stuck in a local minima.
            // One reasonable thing to do here is to compare the variance of both projections,
            // and keep the one with the lowest variance.
            // Of course, one problem is that if the spin component is already in one sector,
            // projecting to the other sector will zero the norm. So we can only make this
            // decision if the |spin component| << 1. Maybe < 0.5 is enough?
            auto spin_component_along_requested_axis = tools::finite::measure::spin_component(*tensors.state, target_sector.value());
            tools::log->debug("Spin component along {} = {:.16f}", target_sector.value(), spin_component_along_requested_axis);
            if(std::abs(spin_component_along_requested_axis) < 0.5) {
                // Here we deem the spin component undecided enough to make a safe projection to both sides for comparison
                auto tensors_neg  = tensors;
                auto tensors_pos  = tensors;
                auto variance_neg = std::numeric_limits<double>::quiet_NaN();
                auto variance_pos = std::numeric_limits<double>::quiet_NaN();
                try {
                    tools::log->debug("Trying projection to -{}", target_sector.value());
                    tensors_neg.project_to_nearest_axis(fmt::format("-{}", target_sector.value()), svd::config(status.bond_lim, status.trnc_lim));
                    variance_neg = tools::finite::measure::energy_variance(tensors_neg);
                } catch(const std::exception &ex) { throw except::runtime_error("Projection to -{} failed: {}", target_sector.value(), ex.what()); }

                try {
                    tools::log->debug("Trying projection to +{}", target_sector.value());
                    tensors_pos.project_to_nearest_axis(fmt::format("+{}", target_sector.value()), svd::config(status.bond_lim, status.trnc_lim));
                    variance_pos = tools::finite::measure::energy_variance(tensors_pos);
                } catch(const std::exception &ex) { throw except::runtime_error("Projection to +{} failed: {}", target_sector.value(), ex.what()); }

                tools::log->debug("Variance after projection to -{} = {:8.2e}", target_sector.value(), variance_neg);
                tools::log->debug("Variance after projection to +{} = {:8.2e}", target_sector.value(), variance_pos);
                if(std::isnan(variance_neg) and std::isnan(variance_pos))
                    tools::log->warn("Both -{0} and +{0} projections failed to yield a valid variance", target_sector.value());

                if(not std::isnan(variance_neg) and variance_neg < variance_pos)
                    tensors = tensors_neg;
                else if(not std::isnan(variance_pos))
                    tensors = tensors_pos;
            } else {
                // Here the spin component is close to one sector. We just project to the nearest sector
                // It may turn out that the spin component is almost exactly +-1 already, then no projection happens, but other
                // routines may go through, such as sign selection on MPO² projection.

                tensors.project_to_nearest_axis(target_sector.value(), svd::config(status.bond_lim, status.trnc_lim));
            }
            tensors.rebuild_edges();
            auto variance_new  = tools::finite::measure::energy_variance(tensors);
            auto spincomp_new  = tools::finite::measure::spin_components(*tensors.state);
            auto entropies_new = tools::finite::measure::entanglement_entropies(*tensors.state);
            tools::log->info("Projection result: variance {:8.2e} -> {:8.2e}  | spin components {:.16f} -> {:.16f}", variance_old, variance_new,
                             fmt::join(spincomp_old, ", "), fmt::join(spincomp_new, ", "));
            if(tools::log->level() <= spdlog::level::debug)
                for(const auto &[i, e] : iter::enumerate(entropies_old)) {
                    tools::log->debug("entropy [{:>2}] = {:>8.6f} --> {:>8.6f} | change {:8.5e}", i, e, entropies_new[i], entropies_new[i] - e);
                }
        }
        if(target_sector.value() == settings::strategy::target_axis) projected_iter = status.iter;
        write_to_file(StorageEvent::PROJ_STATE, CopyPolicy::OFF);
    }
}

void AlgorithmFinite::try_parity_shift() {
    if(not settings::precision::use_mpo_parity_shift) return;
    if(not tensors.position_is_inward_edge()) return;
    if(not qm::spin::half::is_valid_axis(settings::strategy::target_axis)) return;
    tensors.set_parity_shift_mpo_squared(settings::strategy::target_axis);
}

void AlgorithmFinite::try_parity_sep() {
    if(not tensors.position_is_inward_edge()) return;
    size_t iter_since_last_projection = std::max(projected_iter, status.iter) - projected_iter;
    bool   project_on_every_iter = settings::strategy::project_on_every_iter > 0 and iter_since_last_projection >= settings::strategy::project_on_every_iter;

    if(project_on_every_iter and status.iter >= 10) { tensors.set_psfactor(1.0); }
}

AlgorithmFinite::log_entry::log_entry(const AlgorithmStatus &s, const TensorsFinite &t)
    : status(s), variance(status.algo_type == AlgorithmType::fLBIT ? 0.0 : tools::finite::measure::energy_variance(t)),
      entropies(tools::finite::measure::entanglement_entropies(*t.state)) {}

void AlgorithmFinite::check_convergence_variance(std::optional<double> threshold, std::optional<double> saturation_sensitivity) {
    if(not tensors.position_is_inward_edge()) return;
    if(not threshold) threshold = std::max(status.energy_variance_prec_limit, settings::precision::variance_convergence_threshold);
    if(not saturation_sensitivity) saturation_sensitivity = settings::precision::variance_saturation_sensitivity;
    tools::log->trace("Checking convergence of variance mpo | convergence threshold {:.2e} | sensitivity {:.2e}", threshold.value(),
                      saturation_sensitivity.value());

    if(algorithm_history.empty() or algorithm_history.back().status.step < status.step)
        algorithm_history.emplace_back(status, tensors);
    else
        algorithm_history.back() = log_entry(status, tensors);

    // Gather the variance history
    std::vector<double> var_mpo_iter;
    std::transform(algorithm_history.begin(), algorithm_history.end(), std::back_inserter(var_mpo_iter),
                   [](const log_entry &h) -> double { return h.variance; });

    //    var_mpo_iter.emplace_back(tools::finite::measure::energy_variance(tensors));
    auto report = check_saturation(var_mpo_iter, saturation_sensitivity.value(), SaturationScale::log);
    if(report.has_computed) {
        status.variance_mpo_converged_for                          = count_convergence(var_mpo_iter, threshold.value(), report.saturated_point);
        status.variance_mpo_saturated_for                          = report.saturated_count;
        algorithm_history.back().status.variance_mpo_converged_for = status.variance_mpo_converged_for;
        algorithm_history.back().status.variance_mpo_saturated_for = status.variance_mpo_saturated_for;

        if(tools::log->level() >= spdlog::level::debug)
            tools::log->debug("Energy variance convergence: converged {} | saturated {} (since {})", status.variance_mpo_converged_for, report.saturated_count,
                              report.saturated_point);
        if(tools::log->level() <= spdlog::level::trace) {
            tools::log->trace("Energy variance slope details:");
            tools::log->trace(" -- sensitivity        = {:7.4e}", saturation_sensitivity.value());
            tools::log->trace(" -- threshold          = {:7.4e}", threshold.value());
            tools::log->trace(" -- saturated point    = {} ", report.saturated_point);
            tools::log->trace(" -- saturated count    = {} ", report.saturated_count);
            tools::log->trace(" -- converged count    = {} ", status.variance_mpo_converged_for);
            tools::log->trace(" -- var history        = {:7.4e}", fmt::join(report.Y_vec, ", "));
            tools::log->trace(" -- avg history        = {:7.4e}", fmt::join(report.Y_avg, ", "));
            tools::log->trace(" -- std history        = {:7.4e}", fmt::join(report.Y_std, ", "));
            tools::log->trace(" -- sat history        = {}", report.Y_sat, ", ");
        }
    }
}

void AlgorithmFinite::check_convergence_entg_entropy(std::optional<double> saturation_sensitivity) {
    if(not tensors.position_is_inward_edge()) return;
    tools::log->trace("Checking convergence of entanglement");
    if(not saturation_sensitivity) saturation_sensitivity = settings::precision::entropy_saturation_sensitivity;

    if(algorithm_history.empty() or algorithm_history.back().status.step < status.step)
        algorithm_history.emplace_back(status, tensors);
    else
        algorithm_history.back() = log_entry(status, tensors);

    // Gather the entropy history
    size_t                           entropies_size = tensors.get_length() + 1;
    std::vector<SaturationReport>    reports(entropies_size);
    std::vector<std::vector<double>> entropy_iter(entropies_size);

    for(size_t site = 0; site < entropies_size; site++) {
        std::transform(algorithm_history.begin(), algorithm_history.end(), std::back_inserter(entropy_iter[site]),
                       [entropies_size, site](const log_entry &h) -> double {
                           if(h.entropies.empty()) throw except::runtime_error("Entanglement entropies are missing from algorithm history entry");
                           if(h.entropies.size() != entropies_size)
                               throw except::runtime_error("Entanglement entropies have the wrong size {} != {}", h.entropies.size(), entropies_size);
                           return h.entropies[site];
                       });
        reports[site] = check_saturation(entropy_iter[site], saturation_sensitivity.value(), SaturationScale::lin);
    }

    bool all_computed = std::all_of(reports.begin(), reports.end(), [](const SaturationReport &r) { return r.has_computed; });
    if(all_computed) {
        // Find the report which saturated last
        auto last_saturated_itr = std::min_element(reports.begin(), reports.end(), [](const SaturationReport &r1, const SaturationReport &r2) -> bool {
            return r1.saturated_count < r2.saturated_count;
        });
        if(last_saturated_itr != reports.end()) {
            auto  last_saturated_site         = static_cast<size_t>(std::distance(reports.begin(), last_saturated_itr));
            auto &report                      = reports[last_saturated_site];
            status.entanglement_saturated_for = report.saturated_count;
            if(tools::log->level() >= spdlog::level::debug)
                tools::log->debug("Entanglement ent. convergence at site {}: converged {} | saturated {} iters (since {})", last_saturated_site,
                                  status.entanglement_converged_for, report.saturated_count, report.saturated_point);
            if(tools::log->level() <= spdlog::level::trace) {
                tools::log->trace("Entanglement slope details:");
                tools::log->trace(" -- site               = {}", last_saturated_site);
                tools::log->trace(" -- sensitivity        = {:7.4e}", saturation_sensitivity.value());
                tools::log->trace(" -- saturated point    = {} ", report.saturated_point);
                tools::log->trace(" -- saturated count    = {} ", report.saturated_count);
                tools::log->trace(" -- ent history        = {:7.4e}", fmt::join(report.Y_vec, ", "));
                tools::log->trace(" -- avg history        = {:7.4e}", fmt::join(report.Y_avg, ", "));
                tools::log->trace(" -- std history        = {:7.4e}", fmt::join(report.Y_std, ", "));
                tools::log->trace(" -- sat history        = {}", report.Y_sat);
            }
        }
    }
    status.entanglement_converged_for = status.entanglement_saturated_for;

    algorithm_history.back().status.entanglement_converged_for = status.entanglement_converged_for;
    algorithm_history.back().status.entanglement_saturated_for = status.entanglement_saturated_for;
}

void AlgorithmFinite::check_convergence_spin_parity_sector(std::string_view target_sector, double threshold) {
    static constexpr std::array<std::string_view, 9> valid_sectors = {"x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z"};
    bool sector_is_valid = std::find(valid_sectors.begin(), valid_sectors.end(), target_sector) != valid_sectors.end();
    if(sector_is_valid) {
        auto axis                        = qm::spin::half::get_axis_unsigned(settings::strategy::target_axis);
        auto sign                        = qm::spin::half::get_sign(settings::strategy::target_axis);
        auto spin_components             = tools::finite::measure::spin_components(*tensors.state);
        auto spin_component_along_axis   = tools::finite::measure::spin_component(*tensors.state, settings::strategy::target_axis);
        status.spin_parity_has_converged = std::abs(std::abs(spin_component_along_axis) - 1) <= threshold;
        if(status.spin_parity_has_converged and spin_component_along_axis * sign < 0)
            tools::log->warn("Spin component {} has converged: {:.16f} but requested sector was {}", axis, fmt::join(spin_components, ", "), target_sector);
        if(not status.spin_parity_has_converged) {
            tools::log->info("Spin component {} not converged: {:.16f} | threshold {:8.2e}", target_sector, fmt::join(spin_components, ", "), threshold);
        } else {
            tools::log->debug("Spin component {} has converged: {:.16f} | threshold {:8.2e}", target_sector, fmt::join(spin_components, ", "), threshold);
        }
    } else
        status.spin_parity_has_converged = true; // Probably no sector was specified
}

void AlgorithmFinite::clear_convergence_status() {
    tools::log->trace("Clearing convergence status");
    algorithm_history.clear();
    status.algo_stop                  = AlgorithmStop::NONE;
    status.algorithm_has_finished     = false;
    status.algorithm_has_succeeded    = false;
    status.algorithm_has_to_stop      = false;
    status.algorithm_has_stuck_for    = 0;
    status.algorithm_saturated_for    = 0;
    status.algorithm_converged_for    = 0;
    status.entanglement_converged_for = 0;
    status.entanglement_saturated_for = 0;
    status.variance_mpo_converged_for = 0;
    status.variance_mpo_saturated_for = 0;
    status.bond_limit_has_reached_max = false;
    status.trnc_limit_has_reached_min = false;
    status.spin_parity_has_converged  = false;
}

void AlgorithmFinite::write_to_file(StorageEvent storage_event, CopyPolicy copy_policy) {
    if(not write_enabled) return;
    tools::finite::h5::save::simulation(*h5file, tensors, status, storage_event, copy_policy);
}

void AlgorithmFinite::write_to_file(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges, StorageEvent storage_event,
                                    CopyPolicy copy_policy) {
    if(not write_enabled) return;
    tools::finite::h5::save::simulation(*h5file, state, model, edges, status, storage_event, copy_policy);
}

template<typename T>
void AlgorithmFinite::write_to_file(const T &data, std::string_view name, StorageEvent storage_event, CopyPolicy copy_policy) {
    if(not write_enabled) return;
    auto sinfo = StorageInfo(status, tensors.state->get_name(), storage_event);
    tools::finite::h5::save::data(*h5file, sinfo, data, name, copy_policy);
}

template void AlgorithmFinite::write_to_file(const Eigen::Tensor<std::complex<double>, 2> &data, std::string_view name, StorageEvent storage_event,
                                             CopyPolicy copy_policy);

void AlgorithmFinite::print_status() {
    if(num::mod(status.step, settings::print_freq(status.algo_type)) != 0) return;
    if(settings::print_freq(status.algo_type) == 0) return;

    std::string report;
    //    report += fmt::format("{:<} ", status.algo_type_sv());
    report += fmt::format(FMT_STRING("{:<} "), tensors.state->get_name());
    report += fmt::format(FMT_STRING("iter:{:<4} "), status.iter);
    report += fmt::format(FMT_STRING("step:{:<5} "), status.step);
    report += fmt::format(FMT_STRING("L:{} "), tensors.get_length());
    std::string site_str;
    if(tensors.active_sites.empty()) site_str = fmt::format(FMT_STRING("{:^6}"), tensors.state->get_position<long>());
    if(tensors.active_sites.size() == 1) site_str = fmt::format(FMT_STRING("{:^6}"), tensors.active_sites.front());
    if(tensors.active_sites.size() >= 2) {
        auto frnt = sites_mps.has_value() ? sites_mps->at(tensors.active_sites.front()) : tensors.active_sites.front();
        auto back = sites_mps.has_value() ? sites_mps->at(tensors.active_sites.back()) : tensors.active_sites.back();
        if(tensors.position_is_at(static_cast<long>(tensors.active_sites.front())))
            site_str = fmt::format(FMT_STRING("{:>2}.{:>2} "), frnt, back);
        else if(tensors.position_is_at(static_cast<long>(tensors.active_sites.back())))
            site_str = fmt::format(FMT_STRING("{:>2} {:>2}."), frnt, back);
        else
            site_str = fmt::format(FMT_STRING("{:>2} {:>2} "), frnt, back);
    }
    if(tensors.state->get_direction() > 0) {
        report += fmt::format("l:|{}⟩ ", site_str);
    } else if(tensors.state->get_direction() < 0) {
        report += fmt::format("l:⟨{}| ", site_str);
    }

    double epsite = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy(tensors);
    report += fmt::format(FMT_STRING("E:{:<20.16f} "), epsite);

    if(status.algo_type == AlgorithmType::xDMRG) { report += fmt::format(FMT_STRING("e:{:<5.3f} "), status.energy_dens); }
    report += fmt::format(FMT_STRING("Sₑ({:>2}):{:<10.8f} "), tensors.state->get_position<long>(),
                          tools::finite::measure::entanglement_entropy_current(*tensors.state));

    double variance = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy_variance(tensors);
    report += fmt::format(FMT_STRING("σ²H:{:<8.2e} [{:<8.2e}] "), variance, status.energy_variance_lowest);
    report += fmt::format(FMT_STRING("ε:{:<8.2e} "), tensors.state->get_truncation_error_active_max());
    if(settings::strategy::multisite_mps_site_def == 1) report += fmt::format(FMT_STRING("α:{:<8.2e} "), status.env_expansion_alpha);
    report += fmt::format(FMT_STRING("χ:{:<3}|{:<3}|"), settings::get_bond_max(status.algo_type), status.bond_lim);
    size_t comma_width       = settings::strategy::multisite_mps_site_max <= 2 ? 0 : 2; // ", "
    size_t bracket_width     = 2;                                                       // The {} edges
    size_t bond_single_width = static_cast<size_t>(std::log10(settings::get_bond_max(status.algo_type))) + 1;
    size_t bond_num_elements = settings::strategy::multisite_mps_site_max == 1 ? 1 : settings::strategy::multisite_mps_site_max - 1;
    size_t bond_string_width = bracket_width + (bond_single_width + comma_width) // The width of an element like " 54,"
                                                   * bond_num_elements;          // Number of bonds
    std::vector<long> bonds_merged = tools::finite::measure::bond_dimensions_merged(*tensors.state);
    if(bonds_merged.empty())
        report += fmt::format(FMT_STRING("{0:<{1}} "), " ", bond_string_width);
    else {
        std::string bonds_string = fmt::format("{}", bonds_merged);
        report += fmt::format(FMT_STRING("{0:<{1}} "), bonds_string, bond_string_width);
    }

    if(last_optmode and last_optspace)
        report += fmt::format(FMT_STRING("opt:[{}|{}] "), enum2sv(last_optmode.value()).substr(0, 3), enum2sv(last_optspace.value()).substr(0, 3));
    report += fmt::format(FMT_STRING("stk:{:<1} "), status.algorithm_has_stuck_for);
    report += fmt::format(FMT_STRING("sat:{:<1}[σ² {:<1} Sₑ {:<1}] "), status.algorithm_saturated_for, status.variance_mpo_saturated_for,
                          status.entanglement_saturated_for);
    report += fmt::format(FMT_STRING("time:{:<9} "), fmt::format("{:>7.1f}s", tid::get_unscoped("t_tot").get_time()));
    report += fmt::format(FMT_STRING("mem[rss {:<.1f}|peak {:<.1f}|vm {:<.1f}]MB "), debug::mem_rss_in_mb(), debug::mem_hwm_in_mb(), debug::mem_vm_in_mb());
    tools::log->info(report);
}

void AlgorithmFinite::print_status_full() {
    tensors.redo_all_measurements();
    tools::log->info("{:=^60}", "");
    tools::log->info("= {: ^56} =", fmt::format("Completed [{}][{}]", status.algo_type_sv(), tensors.state->get_name()));
    tools::log->info("{:=^60}", "");
    tools::log->info("Stop reason                        = {}", status.algo_stop_sv());
    tools::log->info("Sites                              = {}", tensors.get_length());
    tools::log->info("Position                           = {}", status.position);
    tools::log->info("Direction                          = {}", status.direction);
    tools::log->info("Iterations (full chain sweeps)     = {}", status.iter);
    tools::log->info("Steps (moves along the chain)      = {}", status.step);
    tools::log->info("Total time                         = {:<.1f} s = {:<.2f} min", tid::get_unscoped("t_tot").get_time(),
                     tid::get_unscoped("t_tot").get_time() / 60);

    if(status.algo_type != AlgorithmType::fLBIT) {
        double energy = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy(tensors);
        double epsite = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy_per_site(tensors);
        tools::log->info("Energy          E                  = {:<.16f}", energy);
        tools::log->info("Energy per site E/L                = {:<.16f}", epsite);
        if(status.algo_type == AlgorithmType::xDMRG)
            tools::log->info("Energy density (rescaled 0 to 1) ε = {:<6.4f}",
                             tools::finite::measure::energy_normalized(tensors, status.energy_min, status.energy_max));
        double variance = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy_variance(tensors);
        tools::log->info("Energy variance σ²(H)              = {:<8.2e}", variance);
    }
    tools::log->info("Bond dimension maximum χmax        = {}", settings::get_bond_max(status.algo_type));
    tools::log->info("Bond dimensions χ                  = {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->info("Bond dimension  χ (mid)            = {}", tools::finite::measure::bond_dimension_midchain(*tensors.state));
    tools::log->info("Entanglement entropies Sₑ          = {:8.2e}", fmt::join(tools::finite::measure::entanglement_entropies(*tensors.state), ", "));
    tools::log->info("Entanglement entropy   Sₑ (mid)    = {:8.2e}", tools::finite::measure::entanglement_entropy_midchain(*tensors.state), ", ");
    if(status.algo_type == AlgorithmType::fLBIT) {
        tools::log->info("Number entropies Sₙ                = {:8.2e}", fmt::join(tools::finite::measure::number_entropies(*tensors.state), ", "));
        tools::log->info("Number entropy   Sₙ (mid)          = {:8.2e}", tools::finite::measure::number_entropy_midchain(*tensors.state), ", ");
    }
    tools::log->info("Spin components (global X,Y,Z)     = {:.16f}", fmt::join(tools::finite::measure::spin_components(*tensors.state), ", "));

    if(status.algo_type == AlgorithmType::xDMRG) {
        tools::finite::measure::expectation_values_xyz(*tensors.state);
        tools::finite::measure::correlation_matrix_xyz(*tensors.state);
        tools::finite::measure::structure_factors_xyz(*tensors.state);
        tools::finite::measure::kvornings_marker(*tensors.state);
        tools::log->info("Expectation values ⟨σx⟩            = {:+9.6f}",
                         fmt::join(tenx::span(tensors.state->measurements.expectation_values_sx.value()), ", "));
        tools::log->info("Expectation values ⟨σy⟩            = {:+9.6f}",
                         fmt::join(tenx::span(tensors.state->measurements.expectation_values_sy.value()), ", "));
        tools::log->info("Expectation values ⟨σz⟩            = {:+9.6f}",
                         fmt::join(tenx::span(tensors.state->measurements.expectation_values_sz.value()), ", "));
        tools::log->info("Structure f. L⁻¹ ∑_ij ⟨σx_i σx_j⟩² = {:+.16f}", tensors.state->measurements.structure_factor_x.value());
        tools::log->info("Structure f. L⁻¹ ∑_ij ⟨σy_i σy_j⟩² = {:+.16f}", tensors.state->measurements.structure_factor_y.value());
        tools::log->info("Structure f. L⁻¹ ∑_ij ⟨σz_i σz_j⟩² = {:+.16f}", tensors.state->measurements.structure_factor_z.value());
        tools::log->info("Kvornings marker ⟨σ+..σz..σ-⟩      = {:.8f}", fmt::join(tenx::span(tensors.state->measurements.kvornings_marker.value()), ", "));
    }

    tools::log->info("Truncation Error limit             = {:8.2e}", status.trnc_lim);
    tools::log->info("Truncation Errors ε                = {:8.2e}", fmt::join(tensors.state->get_truncation_errors(), ", "));
    tools::log->info("Algorithm has succeeded            = {:<}", status.algorithm_has_succeeded);
    tools::log->info("Algorithm has saturated for        = {:<}", status.algorithm_saturated_for);
    tools::log->info("Algorithm has got stuck for        = {:<}", status.algorithm_has_stuck_for);
    tools::log->info("Algorithm has converged for        = {:<}", status.algorithm_converged_for);

    if(status.algo_type != AlgorithmType::fLBIT) {
        tools::log->info("σ²                                 = Converged : {:<4}  Saturated: {:<4}", status.variance_mpo_converged_for,
                         status.variance_mpo_saturated_for);
    }
    tools::log->info("Sₑ                                 = Converged : {:<4}  Saturated: {:<4}", status.entanglement_converged_for,
                     status.entanglement_saturated_for);
    tools::log->info("Mem RSS                            = {:<.1f} MB", debug::mem_rss_in_mb());
    tools::log->info("Mem Peak                           = {:<.1f} MB", debug::mem_hwm_in_mb());
    tools::log->info("Mem VM                             = {:<.1f} MB", debug::mem_vm_in_mb());
    tools::log->info("{:=^60}", "");
}
