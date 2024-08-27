#include "AlgorithmFinite.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "debug/info.h"
#include "general/iter.h"
#include "math/cast.h"
#include "math/num.h"
#include "math/stat.h"
#include "math/tenx/span.h"
#include "qm/spin.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/env/EnvEne.h"
#include "tensors/site/env/EnvVar.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/h5/storage_info.h"
#include "tools/common/log.h"
#include "tools/finite/env.h"
#include "tools/finite/env/EnvExpansionResult.h"
#include "tools/finite/h5.h"
#include "tools/finite/measure.h"
#include "tools/finite/mps.h"
#include "tools/finite/multisite.h"
#include "tools/finite/ops.h"
#include "tools/finite/print.h"
#include <h5pp/h5pp.h>

AlgorithmFinite::AlgorithmFinite(OptRitz opt_ritz, AlgorithmType algo_type) : AlgorithmBase(opt_ritz, algo_type) {
    tools::log->trace("Constructing class_algorithm_finite");
    tensors.initialize(algo_type, settings::model::model_type, settings::model::model_size, 0);
    dmrg_degeneracy_score.resize(tensors.get_length<size_t>(), std::numeric_limits<double>::quiet_NaN());
    dmrg_eigs_tol = settings::precision::eigs_tol_max;
}

AlgorithmFinite::AlgorithmFinite(std::shared_ptr<h5pp::File> h5ppFile_, OptRitz opt_ritz, AlgorithmType algo_type)
    : AlgorithmBase(std::move(h5ppFile_), opt_ritz, algo_type) {
    tools::log->trace("Constructing class_algorithm_finite");
    tensors.initialize(algo_type, settings::model::model_type, settings::model::model_size, 0);
    dmrg_degeneracy_score.resize(tensors.get_length<size_t>(), std::numeric_limits<double>::quiet_NaN());
    dmrg_eigs_tol = settings::precision::eigs_tol_max;
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
    tools::log->info("Starting {} simulation", status.algo_type_sv());
    auto t_tot  = tid::get_unscoped("t_tot", tid::level::normal).tic_token();
    auto t_algo = tid::tic_scope(status.algo_type_sv());
    tid::set_level(settings::timer::level);
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

void AlgorithmFinite::run_rbds_analysis() {
    if(settings::strategy::rbds_rate == 0) return;
    last_optcost   = std::nullopt;
    last_optsolver = std::nullopt;
    tools::log     = tools::Logger::setLogger(fmt::format("{}-rbds", status.algo_type_sv()), settings::console::loglevel, settings::console::timestamp);
    tools::log->info("Starting {} reverse bond dimension scaling with rate bond rate [{}] of model [{}] for state [{}]", status.algo_type_sv(),
                     settings::strategy::rbds_rate, enum2sv(settings::model::model_type), tensors.state->get_name());
    auto t_rbds         = tid::tic_scope("rbds");
    auto tensors_backup = tensors;
    auto status_backup  = status;
    tensors.move_center_point_to_inward_edge();
    tensors.activate_sites({tensors.get_position()});
    tensors.rebuild_edges();
    // Generate a list of bond dimension limits
    auto bond_max    = std::min(settings::get_bond_max(status.algo_type), safe_cast<long>(std::pow(2.0, tensors.get_length<double>() / 2)));
    auto bond_limits = std::vector<long>();

    if(settings::strategy::rbds_rate < 1)
        for(long b = bond_max; b > 0l; b = safe_cast<long>(settings::strategy::rbds_rate * safe_cast<double>(b))) bond_limits.emplace_back(b);
    else {
        for(long b = 1l; b <= bond_max; b = safe_cast<long>(settings::strategy::rbds_rate + safe_cast<double>(b))) bond_limits.emplace_back(b);
        std::reverse(bond_limits.begin(), bond_limits.end()); // Needs to go from high to low
    }

    for(const auto &bond_lim : bond_limits) {
        status.bond_lim = bond_lim;
        if(bond_lim < tensors.state->get_largest_bond()) {
            // Cut down the bond dimension with SVDs only
            tools::finite::mps::normalize_state(*tensors.state, svd::config(bond_lim, status.trnc_min), NormPolicy::ALWAYS);
            tensors.clear_cache();
            tensors.clear_measurements();
            tensors.rebuild_edges();

            // Retain the max truncation error, otherwise it is lost in the next pass
            // auto truncation_errors = tensors.state->get_truncation_errors();
            // tensors.state->keep_max_truncation_errors(truncation_errors);
        }
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_rbds->get_time();
        write_to_file(StorageEvent::RBDS_STEP, CopyPolicy::OFF);
        print_status();
    }
    // Reset our logger
    tools::log = tools::Logger::getLogger(std::string(status.algo_type_sv()));
    tensors    = tensors_backup;
    status     = status_backup;
}

void AlgorithmFinite::run_rtes_analysis() {
    if(settings::strategy::rtes_rate <= 1.0) return;
    last_optcost   = std::nullopt;
    last_optsolver = std::nullopt;
    tools::log     = tools::Logger::setLogger(fmt::format("{}-rtes", status.algo_type_sv()), settings::console::loglevel, settings::console::timestamp);
    tools::log->info("Starting {} reverse truncation error scaling with growth rate [{}] of model [{}] for state [{}]", status.algo_type_sv(),
                     settings::strategy::rtes_rate, enum2sv(settings::model::model_type), tensors.state->get_name());
    auto t_rbds         = tid::tic_scope("rtes");
    auto tensors_backup = tensors;
    auto status_backup  = status;
    tensors.move_center_point_to_inward_edge();
    tensors.activate_sites({tensors.get_position()});
    tensors.rebuild_edges();
    // Generate a list of truncation error limits
    auto trnc_min    = std::min(1e-1, settings::precision::svd_truncation_min);
    auto trnc_limits = std::vector<double>();
    for(double t = trnc_min; t < 1e-1 + 1e-8; t *= settings::strategy::rtes_rate) trnc_limits.emplace_back(t);

    for(const auto &trnc_lim : trnc_limits) {
        status.trnc_lim = trnc_lim;
        if(trnc_lim > tensors.state->get_smallest_schmidt_value()) {
            // Cut down the bond dimension with SVDs only
            tools::finite::mps::normalize_state(*tensors.state, svd::config(status.bond_max, trnc_lim), NormPolicy::ALWAYS);
            tensors.clear_cache();
            tensors.clear_measurements();
            tensors.rebuild_edges();
            // Retain the max truncation error, otherwise it is lost in the next pass
            // auto truncation_errors = tensors.state->get_truncation_errors();
            // tensors.state->keep_max_truncation_errors(truncation_errors);
        }
        status.wall_time = tid::get_unscoped("t_tot").get_time();
        status.algo_time = t_rbds->get_time();
        write_to_file(StorageEvent::RTES_STEP, CopyPolicy::OFF);
        print_status();
    }

    // Reset our logger
    tools::log = tools::Logger::getLogger(std::string(status.algo_type_sv()));
    tensors    = tensors_backup;
    status     = status_backup;
}

void AlgorithmFinite::run_postprocessing() {
    tools::log->info("Running default postprocessing for {}", status.algo_type_sv());
    auto tic = tid::tic_scope("post");
    if(has_flag(settings::strategy::projection_policy, ProjectionPolicy::FINISHED)) {
        tensors.project_to_nearest_axis(settings::strategy::target_axis, svd::config(status.bond_lim, status.trnc_lim));
        tensors.rebuild_edges();
    }
    write_to_file(StorageEvent::FINISHED, CopyPolicy::FORCE); // For final mps
    print_status_full();
    run_rbds_analysis();
    run_rtes_analysis();
    tools::log->info("Finished default postprocessing for {}", status.algo_type_sv());
}

int AlgorithmFinite::get_eigs_iter_max() const {
    using namespace settings::precision;
    if(settings::strategy::etol_decrease_when != UpdatePolicy::NEVER) return safe_cast<int>(eigs_iter_max);
    auto iter_min = std::min(safe_cast<double>(eigs_iter_min), safe_cast<double>(eigs_iter_max));
    auto iter_max = std::max(safe_cast<double>(eigs_iter_min), safe_cast<double>(eigs_iter_max));
    if(eigs_iter_gain_policy == EigsIterGainPolicy::NEVER) return safe_cast<int>(iter_min);
    if(has_flag(eigs_iter_gain_policy, EigsIterGainPolicy::FIN_BOND) and !status.bond_limit_has_reached_max) return safe_cast<int>(iter_min);
    if(has_flag(eigs_iter_gain_policy, EigsIterGainPolicy::FIN_TRNC) and !status.trnc_limit_has_reached_min) return safe_cast<int>(iter_min);
    bool   gain_iteration   = has_flag(eigs_iter_gain_policy, EigsIterGainPolicy::ITERATION) and status.iter > 0;
    bool   gain_varsat      = has_flag(eigs_iter_gain_policy, EigsIterGainPolicy::VARSAT) and status.variance_mpo_saturated_for > 0;
    bool   gain_saturated   = has_flag(eigs_iter_gain_policy, EigsIterGainPolicy::SATURATED) and status.algorithm_saturated_for > 0;
    bool   gain_stuck       = has_flag(eigs_iter_gain_policy, EigsIterGainPolicy::STUCK) and status.algorithm_has_stuck_for > 0;
    size_t iteration        = gain_iteration ? status.iter : 0;
    size_t iter_varsat      = gain_varsat ? status.variance_mpo_saturated_for : 0;
    size_t iter_saturated   = gain_saturated ? status.algorithm_saturated_for : 0;
    size_t iter_stuck       = gain_stuck ? status.algorithm_has_stuck_for : 0;
    size_t iter_no_progress = std::max({iteration, iter_varsat, iter_saturated, iter_stuck});
    double iter_multiplier  = std::max(1.0, safe_cast<double>(std::pow(eigs_iter_gain, iter_no_progress)));
    double iter_extra       = iter_min * iter_multiplier; // Run more when stuck
    auto   result           = std::min(iter_extra, iter_max);
    if(result != std::clamp(result, iter_min, iter_max))
        throw except::runtime_error("eigs iter out of range: {} | min {} | max {}", result, eigs_iter_min, eigs_iter_max);
    return safe_cast<int>(result);
}

AlgorithmFinite::OptMeta AlgorithmFinite::get_opt_meta() {
    tools::log->trace("get_opt_meta: configuring optimization step");
    OptMeta m1;
    m1.label = "opt meta";

    // The first decision is easy. Real or complex optimization
    if(tensors.is_real()) m1.optType = OptType::REAL;
    // Set the target eigenvalue
    m1.optRitz = status.opt_ritz;

    // Set the default svd limits
    m1.svd_cfg = svd::config(status.bond_lim, status.trnc_lim);
    // Set up handy quantities
    auto var_conv_thresh = std::max(status.energy_variance_prec_limit, settings::precision::variance_convergence_threshold);
    bool var_isconverged = status.variance_mpo_converged_for > 0 or var_latest < var_conv_thresh;

    // When nearly degenerate we prefer many iterations over many sites
    auto position       = static_cast<size_t>(std::max<long>(0l, tensors.get_position<long>()));
    bool has_degeneracy = dmrg_degeneracy_score[position] >= 0.5 * settings::xdmrg::try_shifting_when_degen;

    // Copy settings
    m1.min_sites         = std::min(tensors.get_length<size_t>(), settings::strategy::dmrg_min_blocksize);
    m1.max_sites         = m1.min_sites;
    m1.subspace_tol      = settings::precision::target_subspace_error;
    m1.primme_projection = "primme_proj_refined"; // converges as [refined < harmonic < RR] (in iterations) (sometimes by a lot) with ~1-10% more time
    m1.primme_method     = "PRIMME_DYNAMIC";
    // m1.primme_method     = "PRIMME_GD_Olsen_plusK";
    m1.eigs_nev          = settings::precision::eigs_nev_min;
    if(status.algorithm_has_stuck_for > 0) {
        m1.eigs_nev = std::clamp(static_cast<int>(1 + status.algorithm_has_stuck_for), settings::precision::eigs_nev_min, settings::precision::eigs_nev_max);
    }
    m1.eigv_target = 0.0;                           // We always target 0 when OptCost::VARIANCE
    m1.eigs_ncv    = settings::precision::eigs_ncv; // Sets to log2(problem_size) if a nonpositive value is given here
    // #pragma message "Revert the envexpandmode to var"
    // if(status.algo_type == AlgorithmType::xDMRG and status.algorithm_has_stuck_for > 0) { m1.expand_mode = EnvExpandMode::ENE | EnvExpandMode::VAR; }
    // Set up a multiplier for number of iterations
    m1.eigs_iter_max = get_eigs_iter_max() * m1.eigs_nev.value();
    m1.eigs_tol      = dmrg_eigs_tol;
    if(settings::strategy::etol_decrease_when == UpdatePolicy::NEVER) {
        // Increase the precision of the eigenvalue solver as variance decreases, and when stuck
        m1.eigs_tol = std::clamp(status.energy_variance_lowest * 1e-1, settings::precision::eigs_tol_min, settings::precision::eigs_tol_max);
        if(status.algorithm_has_stuck_for > 1 and status.bond_limit_has_reached_max and status.trnc_limit_has_reached_min)
            m1.eigs_tol = settings::precision::eigs_tol_min;
    }

    // Set up the subspace (environment) expansion
    m1.expand_side = EnvExpandSide::FORWARD;
    m1.expand_mode = status.algo_type == AlgorithmType::xDMRG ? EnvExpandMode::VAR : EnvExpandMode::ENE;
    // m1.expand_mode = status.algo_type == AlgorithmType::xDMRG ? EnvExpandMode::ENE : EnvExpandMode::ENE;

    m1.optCost = status.algo_type == AlgorithmType::xDMRG ? OptCost::VARIANCE : OptCost::ENERGY;
    m1.optAlgo = OptAlgo::DIRECT;

    m1.optSolver = OptSolver::EIGS;

    if(status.iter < settings::strategy::iter_max_warmup) {
        // If early in the simulation we can use more sites with lower bond dimension o find a good starting point
        m1.max_problem_size = settings::precision::eig_max_size; // Try to use full diagonalization instead (may resolve level spacing issues early on)
        m1.max_sites        = settings::strategy::dmrg_max_blocksize;
    } else {
        m1.max_problem_size = settings::strategy::dmrg_max_prob_size;
        m1.max_sites        = has_degeneracy ? 1 : dmrg_blocksize; // When nearly degenerate we prefer many iterations over many sites
    }

    if(status.algo_type == AlgorithmType::xDMRG and settings::xdmrg::try_directx2_when_stuck and tensors.get_position<long>() >= 0) {
        if(has_degeneracy and not var_isconverged and var_latest < 1e-7) {
            // Try this desperate measure for one iteration...
            // m1.max_sites        = 1; // tensors.get_length<size_t>();
            // m1.max_problem_size = settings::strategy::dmrg_max_prob_size;
            ; // 2 << (m1.max_sites - 1);
            m1.optAlgo = OptAlgo::DIRECTX2;
            // m1.eigs_iter_max = 5000000;
        }
    }

    // Set up the problem size here
    m1.chosen_sites = tools::finite::multisite::generate_site_list(*tensors.state, m1.max_problem_size, m1.max_sites, m1.min_sites, m1.label);
    m1.problem_dims = tools::finite::multisite::get_dimensions(*tensors.state, m1.chosen_sites);
    m1.problem_size = tools::finite::multisite::get_problem_size(*tensors.state, m1.chosen_sites);

    // Do eig instead of eigs when it's cheap (e.g. near the edges or early in the simulation)
    m1.optSolver = m1.problem_size <= settings::precision::eig_max_size ? OptSolver::EIG : OptSolver::EIGS;

    bool full_system_optimization = m1.chosen_sites.size() == tensors.get_length<size_t>();
    if(status.algo_type == AlgorithmType::xDMRG and                                   //
       status.algorithm_has_stuck_for > 1 and                                         //
       status.variance_mpo_saturated_for > settings::strategy::iter_max_saturated and //
       status.trnc_limit_has_reached_min and                                          //
       status.bond_limit_has_reached_max and                                          //
       settings::xdmrg::try_shifting_when_degen > 0 and                               //
       m1.problem_size > settings::precision::eig_max_size and                        //
       not full_system_optimization)                                                  //
    {
        m1.eigs_nev = std::max(m1.eigs_nev.value_or(1), 2);
    }

    if(full_system_optimization and status.algo_type == AlgorithmType::xDMRG and var_latest < 1e-8) {
        m1.optCost = OptCost::ENERGY;
        m1.optRitz = OptRitz::SM;
    }

    // if(m1.optSolver == OptSolver::EIGS and m1.chosen_sites.size() == 1 and tensors.state->get_direction() < 0) {
    // #pragma message "testing DIRECTZ"
    //         m1.eigs_ncv = std::max(m1.eigs_ncv.value_or(32), 32);
    //         m1.eigs_tol = 1e-12;
    //         m1.optAlgo  = OptAlgo::DIRECTZ;
    //     }

    // m1.eigs_nev = settings::precision::eigs_nev_max;
    // m1.eigs_ncv = std::max(m1.eigs_ncv.value_or(32), 32);
    // m1.eigs_tol = 1e-12;
    // m1.optAlgo = OptAlgo::DIRECTX2;
    // m1.optRitz = OptRitz::SM;

    if(var_isconverged) {
        m1.expand_mode = EnvExpandMode::NONE; // No need to expand
    }
    m1.validate();

    return m1;
}

void AlgorithmFinite::expand_environment(EnvExpandMode envexpMode, EnvExpandSide envexpSide, std::optional<svd::config> svd_cfg) {
    // Set a good initial value to start with
    // if(status.env_expansion_alpha == 0.0) status.env_expansion_alpha = std::min(settings::strategy::env_expansion_alpha_max, status.energy_variance_lowest);
    // alpha = alpha.value_or(status.env_expansion_alpha);
    // if(alpha.value() <= 0) return;
    if(not svd_cfg.has_value()) {
        #pragma message "Restore the bond and trnc lims for expansion?"
        auto bond_lim = envexpSide == EnvExpandSide::FORWARD ? status.bond_lim : status.bond_lim;
        auto trnc_lim = envexpSide == EnvExpandSide::FORWARD ? status.trnc_lim : status.trnc_lim;
        svd_cfg       = svd::config(bond_lim, trnc_lim);
    }

    auto res = tensors.expand_environment(envexpMode, envexpSide, svd_cfg.value());
    if(not res.ok) return;
    tools::log->debug("Expanded environment {} {} | bond {} {} | αₑ:{:.2e} αᵥ:{:.2e} | var {:.3e} -> {:.3e} | ene {:.16f} -> {:.16f}", flag2str(envexpMode),
                      enum2sv(envexpSide), res.posL, res.posR, res.alpha_ene, res.alpha_var, res.var_old, res.var_new, res.ene_old, res.ene_new);
    using namespace settings::strategy;
    status.env_expansion_alpha = !std::isnan(res.alpha_ene) ? res.alpha_ene : res.alpha_var;
}

void AlgorithmFinite::move_center_point(std::optional<long> num_moves) {
    auto var_before_move = 0.0;
    if(tensors.state->get_position<long>() >= 0) { var_before_move = tools::finite::measure::energy_variance(tensors); }
    auto old_pos = status.position;
    if(not num_moves.has_value()) {
        if(tensors.active_sites.empty())
            num_moves = 1;
        else {
            long posL_edge     = 0;
            long posR_edge     = tensors.get_length<long>() - 1;
            long posL_active   = safe_cast<long>(tensors.active_sites.front());
            long posR_active   = safe_cast<long>(tensors.active_sites.back());
            long num_active    = safe_cast<long>(tensors.active_sites.size());
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
            } else {
                num_moves = 1;
                if(num_active >= 2 and status.iter < settings::strategy::iter_max_warmup) {
                    num_moves = num_active - 1; // Move so that the center point moves out of the active region
                }
            }
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
            status.step += tensors.move_center_point(svd::config(status.bond_lim, status.trnc_lim)); // Shouldn't truncate.
            // Do not go past the edge if you aren't there already!
            // It's important to stay at the inward edge, so we can do convergence checks and so on
            if(tensors.position_is_inward_edge()) break;
        }
        // if(tensors.position_is_outward_edge())
        // throw except::logic_error("Invalid position after moving {} steps: the position is outward edge: pos {} | dir {}", num_moves,
        // tensors.get_position<long>(), tensors.state->get_direction());
        tensors.clear_active_sites();
    } catch(std::exception &e) {
        tools::finite::print::dimensions(tensors);
        throw except::runtime_error("Failed to move center point: {}", e.what());
    }
    status.position  = tensors.state->get_position<long>();
    status.direction = tensors.state->get_direction();
    if(status.position >= 0) {
        tensors.rebuild_edges();
        tensors.activate_sites({tensors.get_position<size_t>()});
        auto var_after_move = tools::finite::measure::energy_variance(tensors);
        tools::log->debug("Moved center position {} -> {} | var {:.2e} -> {:.2e}", old_pos, status.position, var_before_move, var_after_move);
    }
}

void AlgorithmFinite::set_energy_shift_mpo() {
    if(not settings::precision::use_energy_shifted_mpo) return;
    if(not tensors.position_is_inward_edge()) return;
    auto energy_shift = tools::finite::measure::energy(tensors);
    tensors.set_energy_shift_mpo(energy_shift); // Avoid catastrophic cancellation by shifting energy on each mpo by E/L
}

void AlgorithmFinite::rebuild_tensors() {
    if(not tensors.position_is_inward_edge()) return;
    tensors.rebuild_mpo();          // The shift clears our squared mpo's. So we have to rebuild them.
    tensors.rebuild_mpo_squared();  // The shift clears our squared mpo's. So we have to rebuild them.
    tensors.compress_mpo_squared(); // Compress the mpo's if compression is enabled
    tensors.rebuild_edges();        // The shift modified all our mpo's. So we have to rebuild all the edges.
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
    auto eigs_iter_max_backup          = settings::precision::eigs_iter_max;
    auto eigs_tol_min_backup           = settings::precision::eigs_tol_min;
    auto eigs_ncv_backup               = settings::precision::eigs_ncv;
    auto dmrg_blocksize_backup         = dmrg_blocksize;
    settings::precision::eigs_iter_max = 100;
    settings::precision::eigs_tol_min  = std::min(1e-14, settings::precision::eigs_tol_min);
    settings::precision::eigs_ncv      = std::max(35, settings::precision::eigs_ncv);

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
            tensors.move_site_to_pos(safe_cast<size_t>(site), tgt_pos, sites_mps, sites_mpo, tgt_pos);
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
    settings::precision::eigs_iter_max = eigs_iter_max_backup;
    settings::precision::eigs_tol_min  = eigs_tol_min_backup;
    settings::precision::eigs_ncv      = eigs_ncv_backup;
    dmrg_blocksize                     = dmrg_blocksize_backup;
}

void AlgorithmFinite::update_precision_limit(std::optional<double> energy_upper_bound) {
    if(not tensors.position_is_inward_edge()) return;
    // The variance precision limit depends on the Hamiltonian operator norm ~ largest eigenvalue.
    // We can get a rough order of magnitude etimate the largest eigenvalue by adding the absolute value of all the
    // Hamiltonian couplings and fields.
    if(not energy_upper_bound) energy_upper_bound = tensors.model->get_energy_upper_bound();
    double energy_abs                 = std::abs(energy_upper_bound.value());
    double digits10                   = std::numeric_limits<double>::digits10;
    double energy_exp                 = std::ceil(std::max(0.0, std::log10(energy_abs)));
    double max_digits                 = std::floor(std::max(0.0, digits10 - energy_exp));
    status.energy_variance_max_digits = safe_cast<size_t>(max_digits);
    status.energy_variance_prec_limit = std::pow(10.0, -max_digits);
    tools::log->info("Estimated limit on energy variance precision: {:.3e} | energy upper bound {:.3f}", status.energy_variance_prec_limit,
                     energy_upper_bound.value());
}

void AlgorithmFinite::update_bond_dimension_limit() {
    if(not tensors.position_is_inward_edge()) return;
    status.bond_limit_has_reached_max = status.bond_lim >= status.bond_max;

    if(status.iter < settings::strategy::iter_max_warmup and not has_flag(settings::strategy::bond_increase_when, UpdatePolicy::WARMUP)) {
        tools::log->debug("Kept bond dimension limit: {} | warmup", status.bond_min);
        status.bond_lim = status.bond_min;
        return;
    }
    if(settings::strategy::bond_increase_when == UpdatePolicy::NEVER or settings::strategy::bond_increase_rate <= 1) {
        status.bond_lim                   = status.bond_max;
        status.bond_limit_has_reached_max = true;
        return;
    }
    if(status.bond_limit_has_reached_max) return;
    auto tic = tid::tic_scope("bond_grow", tid::level::higher);

    if constexpr(settings::debug) {
        if(tools::log->level() == spdlog::level::trace) {
            double truncation_threshold = 2 * settings::precision::svd_truncation_min;
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
    bool                          is_on_warmup      = status.iter < settings::strategy::iter_max_warmup;
    bool                          is_truncated      = tensors.state->is_limited_by_bond(status.bond_lim) or tensors.state->is_truncated(status.trnc_lim);
    bool                          is_saturated      = status.algorithm_saturated_for > 0; // Allow one round while saturated so that extra efforts get a chance.
    bool                          is_has_stuck      = status.algorithm_has_stuck_for > 0;
    bool                          is_next_iter      = true; // Should be checked
    bool                          is_even_iter      = status.iter % 2 == 0;
    bool                          grow_if_on_warmup = has_flag(settings::strategy::bond_increase_when, UpdatePolicy::WARMUP);
    bool                          grow_if_next_iter = has_flag(settings::strategy::bond_increase_when, UpdatePolicy::ITERATION);
    bool                          grow_if_even_iter = has_flag(settings::strategy::bond_increase_when, UpdatePolicy::FULLSWEEP);
    bool                          grow_if_truncated = has_flag(settings::strategy::bond_increase_when, UpdatePolicy::TRUNCATED);
    bool                          grow_if_saturated = has_flag(settings::strategy::bond_increase_when, UpdatePolicy::SATURATED);
    bool                          grow_if_has_stuck = has_flag(settings::strategy::bond_increase_when, UpdatePolicy::STUCK);
    std::vector<std::string_view> reason;
    if(grow_if_on_warmup and is_on_warmup) reason.emplace_back("warmup");
    if(grow_if_truncated and is_truncated) reason.emplace_back("truncated");
    if(grow_if_saturated and is_saturated) reason.emplace_back("saturated");
    if(grow_if_has_stuck and is_has_stuck) reason.emplace_back("stuck");
    if(grow_if_next_iter and is_next_iter) reason.emplace_back("next iter");
    if(grow_if_even_iter and is_even_iter) reason.emplace_back("even iter");
    if(reason.empty()) {
        tools::log->trace("Kept bond dimension limit: {} | no reason to update", status.bond_lim);
        return;
    }
    // Write current results before updating bond dimension
    write_to_file(StorageEvent::BOND_UPDATE);

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

    tools::log->info("Updated bond dimension limit: {} -> {} | reason: {}", status.bond_lim, bond_new, fmt::join(reason, " | "));
    status.bond_lim                   = safe_cast<long>(bond_new);
    status.bond_limit_has_reached_max = status.bond_lim == status.bond_max;
    // status.algorithm_has_stuck_for    = 0;
    // status.algorithm_saturated_for    = 0;

    // Last sanity check before leaving here
    if(status.bond_lim > status.bond_max) throw except::logic_error("bond_lim is larger than get_bond_max! {} > {}", status.bond_lim, status.bond_max);
}

void AlgorithmFinite::reduce_bond_dimension_limit(double rate, UpdatePolicy when, StorageEvent storage_event) {
    // We reduce the bond dimension limit during RBDS whenever entanglement has stopped changing
    if(not tensors.position_is_inward_edge()) return;
    if(when == UpdatePolicy::NEVER) return;
    if(iter_last_bond_reduce == 0) iter_last_bond_reduce = status.iter;
    size_t iter_since_reduce = std::max(status.iter, iter_last_bond_reduce) - std::min(status.iter, iter_last_bond_reduce);
    bool   reduce_saturated  = when == UpdatePolicy::SATURATED and status.algorithm_saturated_for > 0;
    bool   reduce_truncated  = when == UpdatePolicy::TRUNCATED and tensors.state->is_truncated(status.trnc_lim);
    bool   reduce_iteration  = when == UpdatePolicy::ITERATION and iter_since_reduce >= 1;
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
        if(bond_new == static_cast<double>(status.bond_lim)) {
            // There would be no change in bond_lim
            status.algo_stop = AlgorithmStop::SUCCESS;
        } else {
            if(storage_event != StorageEvent::NONE) tools::log->info("Updated bond dimension limit {} -> {}", status.bond_lim, bond_new);
            status.bond_lim = safe_cast<long>(bond_new);
        }
        iter_last_bond_reduce = status.iter;
    }
}

void AlgorithmFinite::update_truncation_error_limit() {
    if(not tensors.position_is_inward_edge()) return;
    if(status.trnc_lim == 0.0) throw std::runtime_error("trnc_lim is zero!");
    status.trnc_limit_has_reached_min = status.trnc_lim <= status.trnc_min;

    if(status.iter < settings::strategy::iter_max_warmup and not has_flag(settings::strategy::trnc_decrease_when, UpdatePolicy::WARMUP)) {
        status.trnc_lim = status.trnc_max;
        return;
    }

    if(settings::strategy::trnc_decrease_when == UpdatePolicy::NEVER or settings::strategy::trnc_decrease_rate == 0.0) {
        status.trnc_lim                   = status.trnc_min;
        status.trnc_limit_has_reached_min = true;
        return;
    }

    if(status.trnc_limit_has_reached_min) return;
    auto tic = tid::tic_scope("trnc_down", tid::level::higher);

    if constexpr(settings::debug) {
        if(tools::log->level() == spdlog::level::trace) {
            double truncation_threshold = 2 * settings::precision::svd_truncation_min;
            size_t trunc_bond_count     = tensors.state->num_sites_truncated(truncation_threshold);
            tools::log->trace("Truncation threshold  : {:<.8e}", truncation_threshold);
            tools::log->trace("Truncation errors     : {}", tensors.state->get_truncation_errors());
            tools::log->trace("Bond dimensions       : {}", tools::finite::measure::bond_dimensions(*tensors.state));
            tools::log->trace("Truncated bond count  : {} ", trunc_bond_count);
            tools::log->trace("Entanglement entropies: {} ", tools::finite::measure::entanglement_entropies(*tensors.state));
        }
    }

    // If we got here we want to decrease the truncation error limit progressively during the simulation
    bool                          is_truncated      = tensors.state->is_limited_by_bond(status.bond_lim) or tensors.state->is_truncated(status.trnc_lim);
    bool                          is_saturated      = status.algorithm_saturated_for > 0; // Allow one round while saturated so that extra efforts get a chance.
    bool                          is_has_stuck      = status.algorithm_has_stuck_for > 0;
    bool                          is_next_iter      = true; // Should be checked
    bool                          is_even_iter      = status.iter % 2 == 0;
    bool                          drop_if_truncated = has_flag(settings::strategy::trnc_decrease_when, UpdatePolicy::TRUNCATED);
    bool                          drop_if_saturated = has_flag(settings::strategy::trnc_decrease_when, UpdatePolicy::SATURATED);
    bool                          drop_if_has_stuck = has_flag(settings::strategy::trnc_decrease_when, UpdatePolicy::STUCK);
    bool                          drop_if_next_iter = has_flag(settings::strategy::trnc_decrease_when, UpdatePolicy::ITERATION);
    bool                          drop_if_even_iter = has_flag(settings::strategy::trnc_decrease_when, UpdatePolicy::FULLSWEEP);
    std::vector<std::string_view> reason;
    if(drop_if_truncated and is_truncated) reason.emplace_back("truncated");
    if(drop_if_saturated and is_saturated) reason.emplace_back("saturated");
    if(drop_if_has_stuck and is_has_stuck) reason.emplace_back("stuck");
    if(drop_if_next_iter and is_next_iter) reason.emplace_back("next iter");
    if(drop_if_even_iter and is_even_iter) reason.emplace_back("even iter");
    if(reason.empty()) {
        tools::log->trace("Kept truncation error limit:  {:8.2e} | no reason to update", status.trnc_lim);
        return;
    }

    // Write current results before updating the truncation error limit
    write_to_file(StorageEvent::TRNC_UPDATE);

    // If we got to this point we will update the truncation error limit by a factor
    auto rate = settings::strategy::trnc_decrease_rate;
    if(rate > 1.0 or rate < 0) throw except::runtime_error("Error: trnc_decrease_rate == {:8.2e} | must be in [0, 1]");

    auto trnc_new = std::max(status.trnc_min, status.trnc_lim * rate);
    if(trnc_new < status.trnc_min / rate) trnc_new = status.trnc_min; // If the truncation error is close enough to reaching min, just set min.

    tools::log->info("Updated truncation error limit: {:8.2e} -> {:8.2e} | reasons: {}", status.trnc_lim, trnc_new, fmt::join(reason, " | "));
    status.trnc_lim                   = trnc_new;
    status.trnc_limit_has_reached_min = std::abs(status.trnc_lim - status.trnc_min) < std::numeric_limits<double>::epsilon();
    // status.algorithm_has_stuck_for    = 0;
    // status.algorithm_saturated_for    = 0;
    // Last sanity check before leaving here
    if(status.trnc_lim < status.trnc_min) throw except::logic_error("trnc_lim is smaller than trnc_min ! {:8.2e} > {:8.2e}", status.trnc_lim, status.trnc_min);
}

void AlgorithmFinite::update_dmrg_blocksize() {
    auto old_bsize = dmrg_blocksize;
    if(settings::strategy::dmrg_blocksize_policy == BlockSizePolicy::MIN) {
        dmrg_blocksize = settings::strategy::dmrg_min_blocksize;
        if(old_bsize != dmrg_blocksize) {
            tools::log->info("Updated blocksize {} -> {} | policy: {}", old_bsize, dmrg_blocksize, flag2str(settings::strategy::dmrg_blocksize_policy));
        }
        return;
    }
    if(settings::strategy::dmrg_blocksize_policy == BlockSizePolicy::MAX) {
        dmrg_blocksize = settings::strategy::dmrg_max_blocksize;
        if(old_bsize != dmrg_blocksize) {
            tools::log->info("Updated blocksize {} -> {} | policy: {}", old_bsize, dmrg_blocksize, flag2str(settings::strategy::dmrg_blocksize_policy));
        }
        return;
    }
    if(has_flag(settings::strategy::dmrg_blocksize_policy, BlockSizePolicy::FIN_BOND) and !status.bond_limit_has_reached_max) {
        dmrg_blocksize = settings::strategy::dmrg_min_blocksize;
        if(old_bsize != dmrg_blocksize) {
            tools::log->info("Updated blocksize {} -> {} | policy: {}", old_bsize, dmrg_blocksize, flag2str(settings::strategy::dmrg_blocksize_policy));
        }
        return;
    }
    if(has_flag(settings::strategy::dmrg_blocksize_policy, BlockSizePolicy::FIN_TRNC) and !status.trnc_limit_has_reached_min) {
        dmrg_blocksize = settings::strategy::dmrg_min_blocksize;
        if(old_bsize != dmrg_blocksize) {
            tools::log->info("Updated blocksize {} -> {} | policy: {}", old_bsize, dmrg_blocksize, flag2str(settings::strategy::dmrg_blocksize_policy));
        }
        return;
    }

    if(not tensors.position_is_inward_edge()) return;
    using namespace settings::strategy;
    tools::log->trace("Updating blocksize | policy: {}", flag2str(dmrg_blocksize_policy));
    bool has_converged = status.algorithm_converged_for > 0 or var_latest <= settings::precision::variance_convergence_threshold;
    bool has_sat_ent   = has_flag(dmrg_blocksize_policy, BlockSizePolicy::SAT_ENT) and status.entanglement_saturated_for > 0 and not has_converged;
    bool has_sat_icom  = has_flag(dmrg_blocksize_policy, BlockSizePolicy::SAT_ICOM) and infocom_saturated_for > 0 and not has_converged;
    bool has_sat_var   = has_flag(dmrg_blocksize_policy, BlockSizePolicy::SAT_VAR) and status.variance_mpo_saturated_for > 0 and not has_converged;
    bool has_sat_algo  = has_flag(dmrg_blocksize_policy, BlockSizePolicy::SAT_ALGO) and status.algorithm_saturated_for > 0 and not has_converged;
    bool has_stk_algo  = has_flag(dmrg_blocksize_policy, BlockSizePolicy::STK_ALGO) and status.algorithm_has_stuck_for > 0 and not has_converged;
    // bool has_saturated = has_flag(dmrg_blocksize_policy, BlockSizePolicy::SATURATED) and status.algorithm_saturated_for > 0 and not has_converged;
    // bool has_got_stuck = has_flag(dmrg_blocksize_policy, BlockSizePolicy::STUCK) and status.algorithm_has_stuck_for > 0 and not has_converged;
    bool has_no_status = !has_flag(dmrg_blocksize_policy, BlockSizePolicy::SAT_ENT) and  //
                         !has_flag(dmrg_blocksize_policy, BlockSizePolicy::SAT_ICOM) and //
                         !has_flag(dmrg_blocksize_policy, BlockSizePolicy::SAT_VAR) and  //
                         !has_flag(dmrg_blocksize_policy, BlockSizePolicy::SAT_ALGO) and //
                         !has_flag(dmrg_blocksize_policy, BlockSizePolicy::STK_ALGO);
    bool has_to_adjust = has_sat_ent or has_sat_icom or has_sat_var or has_sat_algo or has_stk_algo or has_no_status;

    BlockSizePolicy choice = BlockSizePolicy::MIN;
    if(has_flag(dmrg_blocksize_policy, BlockSizePolicy::MAX)) choice = BlockSizePolicy::MAX;
    if(has_flag(dmrg_blocksize_policy, BlockSizePolicy::ICOM)) choice = BlockSizePolicy::ICOM;

    BlockSizePolicy reasons = BlockSizePolicy::MIN;
    if(has_sat_ent) reasons |= BlockSizePolicy::SAT_ENT;
    if(has_sat_icom) reasons |= BlockSizePolicy::SAT_ICOM;
    if(has_sat_var) reasons |= BlockSizePolicy::SAT_VAR;
    if(has_sat_algo) reasons |= BlockSizePolicy::SAT_ALGO;
    if(has_stk_algo) reasons |= BlockSizePolicy::STK_ALGO;

    std::string msg;
    if(choice == BlockSizePolicy::MAX and has_to_adjust) {
        dmrg_blocksize = dmrg_max_blocksize;
        tools::log->info("Set MAX blocksize {} -> {}: {}", old_bsize, dmrg_blocksize, has_no_status ? "" : flag2str(reasons));
        return;
    }

    bool has_icom      = has_flag(dmrg_blocksize_policy, BlockSizePolicy::ICOM);
    bool has_icomplus1 = has_flag(dmrg_blocksize_policy, BlockSizePolicy::ICOMPLUS1);
    if((has_icom or has_icomplus1) and has_to_adjust) {
        /*! Update the number of sites used in dmrg updates.
         *  We try to match the blocksize to the "information center of mass", obtained from the information lattice.
         *  This should be a good blocksize to the correlations that exist in the state.
         */
        auto   ip        = InfoPolicy{.bits_max_error = -0.5, .eig_max_size = 3200, .svd_max_size = 1024, .svd_trnc_lim = std::max(status.trnc_lim, 1e-6)};
        double icom      = tools::finite::measure::information_center_of_mass(*tensors.state, ip);
        double blocksize = has_icomplus1 ? std::ceil(icom + 1) : std::ceil(icom);
        dmrg_blocksize   = std::clamp<size_t>(static_cast<size_t>(blocksize), dmrg_min_blocksize, dmrg_max_blocksize);
        tools::log->info("Set ICOM{} blocksize {} -> {}: icom = {:.3f}: {}", has_icomplus1 ? "+1" : "", old_bsize, dmrg_blocksize, icom,
                         has_no_status ? "" : flag2str(reasons));
    } else {
        dmrg_blocksize = dmrg_min_blocksize;
        if(old_bsize != dmrg_blocksize) { tools::log->info("Set MIN blocksize {} -> {} (default)", old_bsize, dmrg_blocksize); }
    }
}

void AlgorithmFinite::update_eigs_tolerance() {
    using namespace settings::precision;
    if(not tensors.position_is_inward_edge()) return;
    dmrg_eigs_tol = std::clamp(dmrg_eigs_tol, eigs_tol_min, eigs_tol_max);

    if(dmrg_eigs_tol == 0.0) throw std::runtime_error("dmrg_eigs_tol is zero!");
    if(status.iter < settings::strategy::iter_max_warmup and not has_flag(settings::strategy::etol_decrease_when, UpdatePolicy::WARMUP)) {
        dmrg_eigs_tol = eigs_tol_max;
        return;
    }

    if(settings::strategy::etol_decrease_when == UpdatePolicy::NEVER or settings::strategy::etol_decrease_rate == 0.0) {
        dmrg_eigs_tol = eigs_tol_min;
        return;
    }
    auto tic = tid::tic_scope("etol_down", tid::level::higher);

    // If we got here we want to decrease the eigs tolerance progressively during the simulation
    bool is_saturated      = status.algorithm_saturated_for > 0; // Allow one round while saturated so that extra efforts get a chance.
    bool is_has_stuck      = status.algorithm_has_stuck_for > 0;
    bool is_next_iter      = true; // Should be checked
    bool is_even_iter      = status.iter % 2 == 0;
    bool drop_if_saturated = has_flag(settings::strategy::etol_decrease_when, UpdatePolicy::SATURATED);
    bool drop_if_has_stuck = has_flag(settings::strategy::etol_decrease_when, UpdatePolicy::STUCK);
    bool drop_if_next_iter = has_flag(settings::strategy::etol_decrease_when, UpdatePolicy::ITERATION);
    bool drop_if_even_iter = has_flag(settings::strategy::etol_decrease_when, UpdatePolicy::FULLSWEEP);

    std::vector<std::string_view> reason;
    if(drop_if_saturated and is_saturated) reason.emplace_back("saturated");
    if(drop_if_has_stuck and is_has_stuck) reason.emplace_back("stuck");
    if(drop_if_next_iter and is_next_iter) reason.emplace_back("next iter");
    if(drop_if_even_iter and is_even_iter) reason.emplace_back("even iter");
    if(reason.empty()) {
        tools::log->trace("Kept truncation error limit:  {:8.2e} | no reason to update", status.trnc_lim);
        return;
    }

    // Write current results before updating
    // write_to_file(StorageEvent::TRNC_UPDATE);

    // If we got to this point we will update by a factor
    auto rate = settings::strategy::etol_decrease_rate;
    if(rate > 1.0 or rate < 0) throw except::runtime_error("Error: etol_decrease_rate == {:8.2e} | must be in [0, 1]");

    auto etol_new = std::max(eigs_tol_min, dmrg_eigs_tol * rate);
    if(etol_new < eigs_tol_min / rate) etol_new = eigs_tol_min; // If the tolerance is close enough to reaching min, just set min.

    tools::log->info("Updated eigs tolerance: {:8.2e} -> {:8.2e} | reasons: {}", dmrg_eigs_tol, etol_new, fmt::join(reason, " | "));
    dmrg_eigs_tol = etol_new;
    // status.trnc_lim                   = trnc_new;
    // status.trnc_limit_has_reached_min = std::abs(status.trnc_lim - status.trnc_min) < std::numeric_limits<double>::epsilon();
    // status.algorithm_has_stuck_for = 0;
    // status.algorithm_saturated_for = 0;
    // Last sanity check before leaving here
    if(dmrg_eigs_tol < eigs_tol_min) throw except::logic_error("dmrg_eigs_tol is smaller than eigs_tol_min ! {:8.2e} > {:8.2e}", dmrg_eigs_tol, eigs_tol_min);
}

void AlgorithmFinite::initialize_model() {
    tools::log->info("Initializing model");
    tensors.initialize_model();
    clear_convergence_status();
    tools::finite::print::model(*tensors.model);
}

void AlgorithmFinite::initialize_state(ResetReason reason, StateInit state_init, std::optional<StateInitType> state_type, std::optional<std::string> axis,
                                       std::optional<bool> use_eigenspinors, std::optional<std::string> pattern, std::optional<long> bond_lim,
                                       std::optional<double> trnc_lim) {
    auto t_rnd = tid::tic_scope("rnd_state", tid::level::higher);
    if(not state_type) state_type = tensors.state->is_real() ? StateInitType::REAL : StateInitType::CPLX;
    if(not axis) axis = settings::strategy::initial_axis;
    if(not use_eigenspinors) use_eigenspinors = settings::strategy::use_eigenspinors;
    if(not pattern) pattern = settings::strategy::initial_pattern;
    if(not bond_lim) {
        bond_lim = settings::get_bond_min(status.algo_type);
        if(settings::strategy::bond_increase_when == UpdatePolicy::NEVER and state_init == StateInit::RANDOMIZE_PREVIOUS_STATE)
            bond_lim = safe_cast<long>(std::pow(2, std::floor(std::log2(tensors.state->get_largest_bond())))); // Nearest power of two from below
    }
    if(not trnc_lim) {
        trnc_lim = settings::precision::svd_truncation_max;
        if(state_init == StateInit::RANDOMIZE_PREVIOUS_STATE) trnc_lim = 1e-2;
    }

    tensors.activate_sites(settings::precision::eigs_max_size_shift_invert, 2); // Activate a pair of sites so that asserts and measurements work
    tensors.rebuild_edges();
    tensors.initialize_state(reason, state_init, state_type.value(), axis.value(), use_eigenspinors.value(), bond_lim.value(), pattern.value());
    if(has_flag(settings::strategy::projection_policy, ProjectionPolicy::INIT) and qm::spin::half::is_valid_axis(axis.value())) {
        tools::log->info("Projecting state | target sector {} | norm {:.16f} | spin components: {::+.16f}", axis.value(),
                         tools::finite::measure::norm(*tensors.state), tools::finite::measure::spin_components(*tensors.state));
        tensors.project_to_nearest_axis(axis.value(), svd::config(bond_lim, trnc_lim));
        tensors.rebuild_edges();
        // Note! After running this function we should rebuild edges! However, there are usually no sites active at this point, so we do it further down.
    }
    settings::strategy::initial_pattern = pattern.value();
    clear_convergence_status();
    status.reset();
    status.iter      = 0;
    status.step      = 0;
    status.position  = tensors.state->get_position<long>();
    status.direction = tensors.state->get_direction();
    status.algo_stop = AlgorithmStop::NONE;
    if(settings::strategy::bond_increase_when != UpdatePolicy::NEVER) status.bond_lim = bond_lim.value();
    if(tensors.state->get_largest_bond() > bond_lim.value())
        //        tools::log->warn("Faulty truncation after randomize. Max found bond is {}, but bond limit is {}", tensors.state->get_largest_bond(),
        //        bond_lim.value());
        throw except::runtime_error("Faulty truncation after randomize. Max found bond dimension is {}, but bond limit is {}",
                                    tensors.state->get_largest_bond(), bond_lim.value());

    tensors.rebuild_edges();
    tools::log->info("State initialization successful:");
    tools::log->info("-- name          : {}", tensors.state->get_name());
    tools::log->info("-- type          : {}", enum2sv(state_init));
    tools::log->info("-- value         : {}", enum2sv(state_type.value()));
    tools::log->info("-- axis          : {}", axis.value());
    tools::log->info("-- pattern       : {}", settings::strategy::initial_pattern);
    tools::log->info("-- labels        : {}", tensors.state->get_labels());
    tools::log->info("-- norm          : {:.16f}", tools::finite::measure::norm(*tensors.state));
    tools::log->info("-- spin (X,Y,Z)  : {::.16f}", tools::finite::measure::spin_components(*tensors.state));
    tools::log->info("-- bond dimensions          : {}", tools::finite::measure::bond_dimensions(*tensors.state));

    if(status.algo_type != AlgorithmType::fLBIT) {
        tools::log->info("-- energy                   : {}", tools::finite::measure::energy(tensors));
        if(!std::isnan(status.energy_min + status.energy_max))
            tools::log->info("-- energy density           : {}", tools::finite::measure::energy_normalized(tensors, status.energy_min, status.energy_max));
        tools::log->info("-- energy variance          : {:8.2e}", tools::finite::measure::energy_variance(tensors));
    }
    write_to_file(StorageEvent::INIT);
}

void AlgorithmFinite::try_projection(std::optional<std::string> target_axis) {
    if(not tensors.position_is_inward_edge()) return;
    if(settings::strategy::projection_policy == ProjectionPolicy::NEVER) return;
    if(status.spin_parity_has_converged) return; // No need
    size_t iter_since_last_projection = std::max(projected_iter, status.iter) - projected_iter;
    if(iter_since_last_projection == 0) return;

    using namespace settings;
    using enum ProjectionPolicy;
    bool projection_if_stuck = has_flag(strategy::projection_policy, STUCK) and status.algorithm_has_stuck_for > 0;
    bool projection_if_conv  = has_flag(strategy::projection_policy, CONVERGED) and status.algorithm_converged_for > 0;
    bool projection_if_iter  = has_flag(strategy::projection_policy, ITER) and iter_since_last_projection > 0;

    if(projection_if_stuck or projection_if_conv or projection_if_iter or target_axis.has_value()) {
        if(not target_axis) target_axis = settings::strategy::target_axis;
        if(not qm::spin::half::is_valid_axis(target_axis.value())) return; // Do not project unless the target sector is one of +- xyz
        std::vector<std::string> msg;
        if(target_axis.has_value()) msg.emplace_back(fmt::format("given sector {}", target_axis.value()));
        if(projection_if_stuck) msg.emplace_back("stuck");
        if(projection_if_conv) msg.emplace_back("converged");
        if(projection_if_iter) msg.emplace_back("iter");
        tools::log->info("Trying projection to {} | reasons: {}", target_axis.value(), fmt::join(msg, " | "));

        auto sector_sign   = qm::spin::half::get_sign(target_axis.value());
        auto energy_old    = tools::finite::measure::energy(tensors);
        auto variance_old  = tools::finite::measure::energy_variance(tensors);
        auto spincomp_old  = tools::finite::measure::spin_components(*tensors.state);
        auto entropies_old = tools::finite::measure::entanglement_entropies(*tensors.state);
        auto svd_cfg       = svd::config(status.bond_lim, status.trnc_lim);
        if(sector_sign != 0) {
            tensors.project_to_nearest_axis(target_axis.value(), svd_cfg);
            tensors.rebuild_edges();
            auto spincomp_new = tools::finite::measure::spin_components(*tensors.state);
            if(spincomp_new != spincomp_old) {
                if(target_axis.value() == settings::strategy::target_axis) projected_iter = status.iter;
                write_to_file(StorageEvent::PROJECTION, CopyPolicy::OFF);
            }
        } else {
            // We have to make a choice if no sector sign has been given, and the spin component along the requested axis is << 1.
            // The simplest thing is to compare the resuls of both projections, e.g. calculate P(+z)|psi> and P(-z)|psi>,
            // and then keep either the one with the lowest variance, or the one with the smallest energy, depending on the selected ritz.
            // Of course, one problem is that if the spin component is already in one sector,
            // projecting to the other sector will zero the norm. So we can only make this
            // decision if the |spin spin_component_along_requested_axis| < 1, otherwise we loose all precision.
            // We choose |spin_component_along_requested_axis| < 0.7 here, but this choice is arbitrary.
            auto spin_component_along_requested_axis = tools::finite::measure::spin_component(*tensors.state, target_axis.value());
            tools::log->debug("Spin component along {} = {:.16f}", target_axis.value(), spin_component_along_requested_axis);
            if(std::abs(spin_component_along_requested_axis) < 0.7) {
                // Here we deem the spin component undecided enough to make a safe projection to both sides for comparison
                auto tensors_neg  = tensors;
                auto tensors_pos  = tensors;
                auto energy_neg   = std::numeric_limits<double>::quiet_NaN();
                auto energy_pos   = std::numeric_limits<double>::quiet_NaN();
                auto variance_neg = std::numeric_limits<double>::quiet_NaN();
                auto variance_pos = std::numeric_limits<double>::quiet_NaN();
                try {
                    tools::log->debug("Trying projection to -{}", target_axis.value());
                    auto target_neg = fmt::format("-{}", target_axis.value());
                    tensors_neg.project_to_nearest_axis(target_neg, svd_cfg);
                    tensors_neg.rebuild_edges();
                    energy_neg   = tools::finite::measure::energy(tensors_neg);
                    variance_neg = tools::finite::measure::energy_variance(tensors_neg);

                } catch(const std::exception &ex) { throw except::runtime_error("Projection to -{} failed: {}", target_axis.value(), ex.what()); }

                try {
                    tools::log->debug("Trying projection to +{}", target_axis.value());
                    auto target_pos = fmt::format("+{}", target_axis.value());
                    tensors_pos.project_to_nearest_axis(target_pos, svd_cfg);
                    tensors_pos.rebuild_edges();
                    energy_pos   = tools::finite::measure::energy(tensors_pos);
                    variance_pos = tools::finite::measure::energy_variance(tensors_pos);
                } catch(const std::exception &ex) { throw except::runtime_error("Projection to +{} failed: {}", target_axis.value(), ex.what()); }

                tools::log->debug("Projection to -{}: Energy: {:.16f} | Variance {:.3e}", target_axis.value(), energy_neg, variance_neg);
                tools::log->debug("Projection to +{}: Energy: {:.16f} | Variance {:.3e}", target_axis.value(), energy_pos, variance_pos);
                if(std::isnan(variance_neg) and std::isnan(variance_pos))
                    tools::log->warn("Both -{0} and +{0} projections failed to yield a valid variance", target_axis.value());

                switch(status.opt_ritz) {
                    case OptRitz::SM: { // Take the smallest |energy|
                        if(std::abs(energy_neg) < std::abs(energy_pos)) {
                            tensors = tensors_neg;
                        } else if(not std::isnan(energy_pos)) {
                            tensors = tensors_pos;
                        }
                        break;
                    }
                    case OptRitz::SR: { // Take the smallest energy
                        if(energy_neg < energy_pos) {
                            tensors = tensors_neg;
                        } else if(not std::isnan(energy_pos)) {
                            tensors = tensors_pos;
                        }
                        break;
                    }
                    case OptRitz::LR: { // Take the largest energy
                        if(energy_neg > energy_pos) {
                            tensors = tensors_neg;
                        } else if(not std::isnan(energy_pos)) {
                            tensors = tensors_pos;
                        }
                        break;
                    }
                    case OptRitz::IS:
                    case OptRitz::TE:
                    case OptRitz::NONE: { // Take the smallest variance
                        if(variance_neg < variance_pos) {
                            tensors = tensors_neg;
                        } else if(not std::isnan(variance_neg)) {
                            tensors = tensors_pos;
                        }
                        break;
                    }
                }

            } else {
                // Here the spin component is close to one sector. We just project to the nearest sector
                // It may turn out that the spin component is almost exactly +-1 already, then no projection happens, but other
                // routines may go through, such as sign selection on MPO² projection.
                tensors.project_to_nearest_axis(target_axis.value(), svd::config(status.bond_lim, status.trnc_lim));
            }
            tensors.rebuild_edges();
            auto energy_new    = tools::finite::measure::energy(tensors);
            auto variance_new  = tools::finite::measure::energy_variance(tensors);
            auto spincomp_new  = tools::finite::measure::spin_components(*tensors.state);
            auto entropies_new = tools::finite::measure::entanglement_entropies(*tensors.state);
            if(spincomp_new != spincomp_old) {
                tools::log->info("Projection result: energy {:.16f} -> {:.16f} variance {:.4e} -> {:.4e}  | spin components {:.16f} -> {:.16f}", energy_old,
                                 energy_new, variance_old, variance_new, fmt::join(spincomp_old, ", "), fmt::join(spincomp_new, ", "));
                if(tools::log->level() <= spdlog::level::debug)
                    for(const auto &[i, e] : iter::enumerate(entropies_old)) {
                        tools::log->debug("entropy [{:>2}] = {:>8.6f} --> {:>8.6f} | change {:8.5e}", i, e, entropies_new[i], entropies_new[i] - e);
                    }
                if(target_axis.value() == settings::strategy::target_axis) projected_iter = status.iter;
                write_to_file(StorageEvent::PROJECTION, CopyPolicy::OFF);
            }
        }
    }
}

void AlgorithmFinite::set_parity_shift_mpo(std::optional<std::string> target_axis) {
    if(not settings::precision::use_parity_shifted_mpo) return;
    if(not tensors.position_is_inward_edge()) return;
    target_axis = target_axis.value_or(settings::strategy::target_axis);
    // If ritz == SR we shift the spectrum of the non-targeted sector UP in energy.
    // If ritz == LR we shift the spectrum of the non-targetd sector DOWN in energy
    tensors.set_parity_shift_mpo(settings::get_ritz(status.algo_type), target_axis.value());
}

void AlgorithmFinite::set_parity_shift_mpo_squared(std::optional<std::string> target_axis) {
    if(not settings::precision::use_parity_shifted_mpo_squared) return;
    if(not tensors.position_is_inward_edge()) return;
    target_axis = target_axis.value_or(settings::strategy::target_axis);
    tensors.set_parity_shift_mpo_squared(target_axis.value());
}

void AlgorithmFinite::check_convergence() {
    if(not tensors.position_is_inward_edge()) return;
    auto t_con = tid::tic_scope("conv");

    check_convergence_variance();
    check_convergence_entg_entropy();
    check_convergence_spin_parity_sector(settings::strategy::target_axis);
    check_convergence_icom();

    // Determine if the algorithm has saturated
    size_t sum_saturated_for = status.variance_mpo_saturated_for + status.entanglement_saturated_for + infocom_saturated_for;
    int    num_enabled       = 0;
    bool   variance_enabled  = settings::precision::variance_saturation_sensitivity > 0;
    bool   entropy_enabled   = settings::precision::entropy_saturation_sensitivity > 0;
    bool   infocom_enabled   = settings::precision::infocom_saturation_sensitivity > 0;
    if(variance_enabled) num_enabled++;
    if(entropy_enabled) num_enabled++;
    if(infocom_enabled) num_enabled++;
    int  variance_saturated = variance_enabled and status.variance_mpo_saturated_for > 1;
    int  entropy_saturated  = entropy_enabled and status.entanglement_saturated_for > 1;
    int  infocom_saturated  = infocom_enabled and infocom_saturated_for > 1;
    int  all_saturated      = variance_saturated + entropy_saturated + infocom_saturated >= num_enabled;
    bool max_saturated      = sum_saturated_for > settings::strategy::iter_max_saturated * num_enabled;
    if(all_saturated or max_saturated) {
        status.algorithm_saturated_for++;
    } else {
        status.algorithm_saturated_for = 0;
    }

    // Determine if the algorithm has converged

    bool variance_mpo_converged = settings::precision::variance_saturation_sensitivity > 0 ? status.variance_mpo_converged_for > 0 : true;
    bool entanglement_converged = settings::precision::entropy_saturation_sensitivity > 0 ? status.entanglement_converged_for > 0 : true;
    if(variance_mpo_converged and entanglement_converged and status.spin_parity_has_converged)
        status.algorithm_converged_for++;
    else
        status.algorithm_converged_for = 0;

    // Determine if the algorithm is stuck
    if(all_saturated and status.algorithm_saturated_for > 1 and status.algorithm_converged_for == 0)
        status.algorithm_has_stuck_for++;
    else
        status.algorithm_has_stuck_for = 0;

    if(status.iter < settings::strategy::iter_max_warmup) {
        status.algorithm_saturated_for = 0;
        status.algorithm_has_stuck_for = 0;
    }

    status.algorithm_has_succeeded = status.algorithm_converged_for > settings::strategy::iter_min_converged;
    status.algorithm_has_to_stop =
        status.bond_limit_has_reached_max and status.trnc_limit_has_reached_min and status.algorithm_has_stuck_for >= settings::strategy::iter_max_stuck;

    tools::log->info("Algorithm report: converged {} (σ² {} Sₑ {} spin {}) | saturated {} (σ² {} Sₑ {} iₗ {}) | stuck {} | succeeded {} | has to stop {} | var "
                     "prec limit {:8.2e}",
                     status.algorithm_converged_for, status.variance_mpo_converged_for, status.entanglement_converged_for, status.spin_parity_has_converged,
                     status.algorithm_saturated_for, status.variance_mpo_saturated_for, status.entanglement_saturated_for, infocom_saturated_for,
                     status.algorithm_has_stuck_for, status.algorithm_has_succeeded, status.algorithm_has_to_stop, status.energy_variance_prec_limit);

    status.algo_stop = AlgorithmStop::NONE;
    if(status.iter >= settings::get_iter_min(status.algo_type)) {
        if(status.iter >= settings::get_iter_max(status.algo_type)) status.algo_stop = AlgorithmStop::MAX_ITERS;
        if(status.algorithm_has_succeeded) status.algo_stop = AlgorithmStop::SUCCESS;
        if(status.algorithm_has_to_stop) status.algo_stop = AlgorithmStop::SATURATED;
    }
}

AlgorithmFinite::log_entry::log_entry(const AlgorithmStatus &s, const TensorsFinite &t) : entropies(tools::finite::measure::entanglement_entropies(*t.state)) {
    status    = s;
    energy    = status.algo_type == AlgorithmType::fLBIT ? 0.0 : tools::finite::measure::energy(t);
    variance  = status.algo_type == AlgorithmType::fLBIT ? 0.0 : tools::finite::measure::energy_variance(t);
    entropies = tools::finite::measure::entanglement_entropies(*t.state);
    auto ip   = InfoPolicy{.bits_max_error = -0.5,
                           .eig_max_size   = 3200,
                           .svd_max_size   = 1024,
                           .svd_trnc_lim   = std::max(status.trnc_lim, 1e-6),
                           .useCache       = UseCache::TRUE};
    icom      = status.algo_type == AlgorithmType::fLBIT ? 0.0 : tools::finite::measure::information_center_of_mass(*t.state, ip);
    time      = status.wall_time;
}

void AlgorithmFinite::check_convergence_variance(std::optional<double> threshold, std::optional<double> saturation_sensitivity) {
    if(not tensors.position_is_inward_edge()) return;
    if(not saturation_sensitivity) saturation_sensitivity = settings::precision::variance_saturation_sensitivity;
    if(saturation_sensitivity <= 0) return;
    if(not threshold) {
        if(status.iter < settings::strategy::iter_max_warmup) {
            threshold = settings::precision::variance_convergence_threshold;
        } else {
            threshold = std::max(status.energy_variance_prec_limit, settings::precision::variance_convergence_threshold);
        }
    }
    saturation_sensitivity = std::max(saturation_sensitivity.value(), std::sqrt(status.trnc_lim)); // Large trnc causes noise that never saturates
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
    auto report = check_saturation(var_mpo_iter, saturation_sensitivity.value(), SaturationPolicy::val | SaturationPolicy::mid | SaturationPolicy::log);
    if(report.has_computed) {
        status.variance_mpo_converged_for                          = count_convergence(var_mpo_iter, threshold.value(), report.saturated_point);
        status.variance_mpo_saturated_for                          = report.saturated_count;
        algorithm_history.back().status.variance_mpo_converged_for = status.variance_mpo_converged_for;
        algorithm_history.back().status.variance_mpo_saturated_for = status.variance_mpo_saturated_for;

        if(tools::log->level() >= spdlog::level::debug)
            tools::log->debug("Energy variance convergence: converged {} | saturated {} (since {})", status.variance_mpo_converged_for, report.saturated_count,
                              report.saturated_point);
        if(tools::log->level() <= spdlog::level::trace) {
            std::vector<double> times;
            std::vector<double> energies;
            std::transform(algorithm_history.begin(), algorithm_history.end(), std::back_inserter(times), [](const log_entry &h) -> double { return h.time; });
            std::transform(algorithm_history.begin(), algorithm_history.end(), std::back_inserter(energies),
                           [](const log_entry &h) -> double { return h.energy; });
            tools::log->trace("Energy variance convergence details:");
            tools::log->trace(" -- sensitivity     = {:7.4e}", saturation_sensitivity.value());
            tools::log->trace(" -- threshold       = {:7.4e}", threshold.value());
            tools::log->trace(" -- saturated point = {} ", report.saturated_point);
            tools::log->trace(" -- saturated count = {} ", report.saturated_count);
            tools::log->trace(" -- converged count = {} ", status.variance_mpo_converged_for);
            tools::log->trace(" -- sat             = {}", report.Y_sat);
            tools::log->trace(" -- val             = {::7.4e}", report.Y_vec);
            tools::log->trace(" -- ene             = {::7.4e}", energies);
            tools::log->trace(" -- time            = {::7.4e}", times);
            tools::log->trace(" -- avg             = {::7.4e}", report.Y_avg);
            tools::log->trace(" -- med             = {::7.4e}", report.Y_med);
            tools::log->trace(" -- min             = {::7.4e}", report.Y_min);
            tools::log->trace(" -- max             = {::7.4e}", report.Y_max);
            tools::log->trace(" -- mid             = {::7.4e}", report.Y_mid);
            tools::log->trace(" -- dif             = {::7.4e}", report.Y_dif);
            tools::log->trace(" -- std_val         = {::7.4e}", report.Y_vec_std);
            tools::log->trace(" -- std_avg         = {::7.4e}", report.Y_avg_std);
            tools::log->trace(" -- std_med         = {::7.4e}", report.Y_med_std);
            tools::log->trace(" -- std_min         = {::7.4e}", report.Y_min_std);
            tools::log->trace(" -- std_max         = {::7.4e}", report.Y_max_std);
            tools::log->trace(" -- std_mid         = {::7.4e}", report.Y_mid_std);
            tools::log->trace(" -- dif_avg         = {::7.4e}", report.Y_dif_avg);
            tools::log->trace(" -- std_mov         = {::7.4e}", report.Y_mov_std);
        }
    }
}

void AlgorithmFinite::check_convergence_icom(std::optional<double> saturation_sensitivity) {
    if(not tensors.position_is_inward_edge()) return;
    if(not saturation_sensitivity) saturation_sensitivity = settings::precision::infocom_saturation_sensitivity;
    if(saturation_sensitivity <= 0) return;
    tools::log->trace("Checking convergence of icom | sensitivity {:.2e}", saturation_sensitivity.value());

    if(algorithm_history.empty() or algorithm_history.back().status.step < status.step)
        algorithm_history.emplace_back(status, tensors);
    else
        algorithm_history.back() = log_entry(status, tensors);

    // Gather the information center of mass history
    std::vector<double> infocom_iter;
    std::transform(algorithm_history.begin(), algorithm_history.end(), std::back_inserter(infocom_iter), [](const log_entry &h) -> double { return h.icom; });

    //    var_mpo_iter.emplace_back(tools::finite::measure::energy_variance(tensors));
    auto report = check_saturation(infocom_iter, saturation_sensitivity.value(), SaturationPolicy::val | SaturationPolicy::mid | SaturationPolicy::mov);
    if(report.has_computed) {
        infocom_saturated_for = report.saturated_count;

        if(tools::log->level() >= spdlog::level::debug)
            tools::log->debug("icom saturation: saturated {} (since {})", infocom_saturated_for, report.saturated_count, report.saturated_point);
        if(tools::log->level() <= spdlog::level::trace) {
            std::vector<double> times;
            std::transform(algorithm_history.begin(), algorithm_history.end(), std::back_inserter(times), [](const log_entry &h) -> double { return h.time; });
            tools::log->trace("ICOM saturation details:");
            tools::log->trace(" -- sensitivity     = {:7.4e}", saturation_sensitivity.value());
            tools::log->trace(" -- saturated point = {} ", report.saturated_point);
            tools::log->trace(" -- saturated count = {} ", report.saturated_count);
            tools::log->trace(" -- converged count = {} ", status.variance_mpo_converged_for);
            tools::log->trace(" -- sat             = {}", report.Y_sat);
            tools::log->trace(" -- val             = {::7.4e}", report.Y_vec);
            tools::log->trace(" -- min             = {::7.4e}", report.Y_min);
            tools::log->trace(" -- max             = {::7.4e}", report.Y_max);
            tools::log->trace(" -- mid             = {::7.4e}", report.Y_mid);
            tools::log->trace(" -- std_var         = {::7.4e}", report.Y_vec_std);
            tools::log->trace(" -- std_min         = {::7.4e}", report.Y_min_std);
            tools::log->trace(" -- std_max         = {::7.4e}", report.Y_max_std);
            tools::log->trace(" -- std_mid         = {::7.4e}", report.Y_mid_std);
            tools::log->trace(" -- ste_mov         = {::7.4e}", report.Y_mov_std);
            tools::log->trace(" -- time            = {::7.4e}", times);
        }
    }
}

void AlgorithmFinite::check_convergence_entg_entropy(std::optional<double> saturation_sensitivity) {
    if(not tensors.position_is_inward_edge()) return;
    if(not saturation_sensitivity) saturation_sensitivity = settings::precision::entropy_saturation_sensitivity;
    if(saturation_sensitivity <= 0) return;
    tools::log->trace("Checking convergence of entanglement");
    if(algorithm_history.empty() or algorithm_history.back().status.step < status.step)
        algorithm_history.emplace_back(status, tensors);
    else
        algorithm_history.back() = log_entry(status, tensors);

    if(status.algo_type == AlgorithmType::fLBIT) {
        status.entanglement_saturated_for                          = 0;
        status.entanglement_converged_for                          = 0;
        algorithm_history.back().status.entanglement_converged_for = 0;
        algorithm_history.back().status.entanglement_saturated_for = 0;
        return;
    }

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
        reports[site] = check_saturation(entropy_iter[site], saturation_sensitivity.value(), SaturationPolicy::val | SaturationPolicy::mov);
    }

    bool all_computed = std::all_of(reports.begin(), reports.end(), [](const SaturationReport &r) { return r.has_computed; });
    if(all_computed) {
        // Find the report which saturated last
        auto last_saturated_itr = std::min_element(reports.begin(), reports.end(), [](const SaturationReport &r1, const SaturationReport &r2) -> bool {
            return r1.saturated_count < r2.saturated_count;
        });
        if(last_saturated_itr != reports.end()) {
            auto  last_saturated_site         = safe_cast<size_t>(std::distance(reports.begin(), last_saturated_itr));
            auto &report                      = reports[last_saturated_site];
            status.entanglement_saturated_for = report.saturated_count;
            if(tools::log->level() >= spdlog::level::debug)
                tools::log->debug("Entanglement ent. convergence at site {}: converged {} | saturated {} iters (since {})", last_saturated_site,
                                  status.entanglement_converged_for, report.saturated_count, report.saturated_point);
            if(tools::log->level() <= spdlog::level::trace) {
                tools::log->trace("Entanglement convergence details:");
                tools::log->trace(" -- site            = {}", last_saturated_site);
                tools::log->trace(" -- sensitivity     = {:7.4e}", saturation_sensitivity.value());
                tools::log->trace(" -- saturated point = {} ", report.saturated_point);
                tools::log->trace(" -- saturated count = {} ", report.saturated_count);
                tools::log->trace(" -- val             = {::7.4e}", report.Y_vec);
                tools::log->trace(" -- min             = {::7.4e}", report.Y_min);
                tools::log->trace(" -- max             = {::7.4e}", report.Y_max);
                tools::log->trace(" -- mid             = {::7.4e}", report.Y_mid);
                tools::log->trace(" -- std_var         = {::7.4e}", report.Y_vec_std);
                tools::log->trace(" -- std_min         = {::7.4e}", report.Y_min_std);
                tools::log->trace(" -- std_max         = {::7.4e}", report.Y_max_std);
                tools::log->trace(" -- std_mid         = {::7.4e}", report.Y_mid_std);
                tools::log->trace(" -- ste_mov         = {::7.4e}", report.Y_mov_std);
            }
        }
    }
    status.entanglement_converged_for                          = status.entanglement_saturated_for;
    algorithm_history.back().status.entanglement_converged_for = status.entanglement_converged_for;
    algorithm_history.back().status.entanglement_saturated_for = status.entanglement_saturated_for;
}

void AlgorithmFinite::check_convergence_spin_parity_sector(std::string_view target_axis, double threshold) {
    static constexpr std::array<std::string_view, 9> valid_sectors = {"x", "+x", "-x", "y", "+y", "-y", "z", "+z", "-z"};
    bool sector_is_valid                                           = std::find(valid_sectors.begin(), valid_sectors.end(), target_axis) != valid_sectors.end();
    if(sector_is_valid) {
        tools::log->trace("Checking convergence of spin components");
        auto target_axus               = qm::spin::half::get_axis_unsigned(target_axis);
        auto target_sign               = qm::spin::half::get_sign(target_axis);
        auto spin_components           = tools::finite::measure::spin_components(*tensors.state);
        auto spin_component_along_axus = tools::finite::measure::spin_component(*tensors.state, target_axus);
        bool spin_along_axus_near_abs1 = std::abs(std::abs(spin_component_along_axus) - 1) <= threshold;
        auto spin_sign_along_axus      = tools::finite::measure::spin_sign(*tensors.state, target_axus);
        // We may have shifted the spin parity sector in the MPO or MPO².
        // If we shifted for +z, that means we are targeting +z and therefore shifted -z out of the way.
        // So we would need to make sure we are actually converging to +z and not -z.
        auto [target_ritz_mpo, target_sign_mpo, target_axus_mpo] = tensors.get_parity_shift_mpo();
        auto [target_sign_mpo2, target_axus_mpo2]                = tensors.get_parity_shift_mpo_squared();

        bool target_axus_mpo_ok  = target_axus_mpo.empty() or target_axus_mpo == target_axus;   // MPO parity shift disabled or same axus
        bool target_sign_mpo_ok  = spin_sign_along_axus * target_sign_mpo >= 0;                 // MPO parity shift disabled or same sign
        bool target_axus_mpo2_ok = target_axus_mpo2.empty() or target_axus_mpo2 == target_axus; // MPO² parity shift disabled or same axus
        bool target_sign_mpo2_ok = spin_sign_along_axus * target_sign_mpo2 >= 0;                // MPO² parity shift disabled or same sign

        bool target_axus_ok = target_axus_mpo_ok and target_axus_mpo2_ok;
        bool target_sign_ok = target_sign_mpo_ok and target_sign_mpo2_ok;

        std::string msg;
        if(!target_axus_mpo_ok or !target_sign_mpo_ok) msg += fmt::format(" mpo {}{}", fmt::format("{:+}", target_sign_mpo).front(), target_axus_mpo);
        if(!target_axus_mpo2_ok or !target_sign_mpo2_ok) msg += fmt::format(" mpo² {}{}", fmt::format("{:+}", target_sign_mpo2).front(), target_axus_mpo2);

        if(!msg.empty()) {
            tools::log->warn("check_convergence_spin_parity_sector: mismatch: target {} | spin parity shift:{} | spin components: {::.16f}", target_axis, msg,
                             spin_components);
            if(spin_along_axus_near_abs1) {
                // We seem to have converged in the opposite parity sector.
                // The algorithm will likely not be able to escape the current sector, so we might as well set it as the target.
                auto target_axis_opposite = fmt::format("{}{}", fmt::format("{:+}", spin_sign_along_axus).front(), target_axus);
                tools::log->warn("check_convergence_spin_parity_sector: resetting spin parity target: {}", target_axis_opposite);
                set_parity_shift_mpo(target_axis_opposite);
                set_parity_shift_mpo_squared(target_axis_opposite);
                rebuild_tensors(); // Rebuilds and compresses mpos, then rebuilds the environments
            }
        }

        status.spin_parity_has_converged = target_axus_ok and target_sign_ok and spin_along_axus_near_abs1;

        if(status.spin_parity_has_converged and spin_component_along_axus * target_sign < 0)
            tools::log->warn("Spin components: {::.16f} | {} converged ({} requested) | threshold {:8.2}", spin_components, target_axus, target_axis,
                             threshold);
        if(not status.spin_parity_has_converged) {
            tools::log->info("Spin components: {::.16f} | {} not converged ({} requested) | threshold {:8.2e}", spin_components, target_axus, target_axis,
                             threshold);
        } else {
            tools::log->debug("Spin components: {::.16f} | {} converged ({} requested) | threshold {:8.2e}", spin_components, target_axus, target_axis,
                              threshold);
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
    infocom_saturated_for             = 0;
    status.bond_limit_has_reached_max = false;
    status.trnc_limit_has_reached_min = false;
    status.spin_parity_has_converged  = false;
    // status.energy_variance_lowest     = 1.0;
    // status.env_expansion_alpha        = 0.0;
}

void AlgorithmFinite::write_to_file(StorageEvent storage_event, CopyPolicy copy_policy) {
    if(not h5file) return;
    status.event = storage_event;
    tools::finite::h5::save::simulation(*h5file, tensors, status, copy_policy);
    status.event = StorageEvent::NONE;
}

void AlgorithmFinite::write_to_file(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges, StorageEvent storage_event,
                                    CopyPolicy copy_policy) {
    if(not h5file) return;
    status.event = storage_event;
    tools::finite::h5::save::simulation(*h5file, state, model, edges, status, copy_policy);
    status.event = StorageEvent::NONE;
}

template<typename T>
void AlgorithmFinite::write_to_file(const T &data, std::string_view name, StorageEvent storage_event, CopyPolicy copy_policy) {
    if(not h5file) return;
    status.event = storage_event;
    auto sinfo   = StorageInfo(status, tensors.state->get_name());
    tools::finite::h5::save::data(*h5file, sinfo, data, name, copy_policy);
    status.event = StorageEvent::NONE;
}

template void AlgorithmFinite::write_to_file(const Eigen::Tensor<cplx, 2> &data, std::string_view name, StorageEvent storage_event, CopyPolicy copy_policy);
template void AlgorithmFinite::write_to_file(const Eigen::Tensor<cplx, 1> &data, std::string_view name, StorageEvent storage_event, CopyPolicy copy_policy);
template void AlgorithmFinite::write_to_file(const Eigen::Tensor<real, 1> &data, std::string_view name, StorageEvent storage_event, CopyPolicy copy_policy);
void          AlgorithmFinite::print_status() {
    if(num::mod(status.step, settings::print_freq(status.algo_type)) != 0) return;
    if(settings::print_freq(status.algo_type) == 0) return;
    std::string report;
    //    report += fmt::format("{:<} ", status.algo_type_sv());
    report += fmt::format("{:<} ", tensors.state->get_name());
    report += fmt::format("iter:{:<4} ", status.iter);
    report += fmt::format("step:{:<5} ", status.step);
    report += fmt::format("L:{} ", tensors.get_length());
    std::string site_str;
    if(tensors.active_sites.empty()) site_str = fmt::format("{:^6}", tensors.state->get_position<long>());
    if(tensors.active_sites.size() == 1) {
        if(tensors.state->get_direction() > 0) {
            site_str = fmt::format("{:^6}", tensors.active_sites.front());
        } else {
            site_str = fmt::format("{:^6}", tensors.active_sites.front());
        }
    }
    if(tensors.active_sites.size() >= 2) {
        auto frnt = sites_mps.has_value() ? sites_mps->at(tensors.active_sites.front()) : tensors.active_sites.front();
        auto back = sites_mps.has_value() ? sites_mps->at(tensors.active_sites.back()) : tensors.active_sites.back();
        if(tensors.position_is_at(safe_cast<long>(tensors.active_sites.front())))
            site_str = fmt::format("{:>2}.{:>2} ", frnt, back);
        else if(tensors.position_is_at(safe_cast<long>(tensors.active_sites.back())))
            site_str = fmt::format("{:>2} {:>2}.", frnt, back);
        else
            site_str = fmt::format("{:>2} {:>2} ", frnt, back);
    }
    if(tensors.state->get_direction() > 0) {
        report += fmt::format("l:|{}⟩ ", site_str);
    } else if(tensors.state->get_direction() < 0) {
        report += fmt::format("l:⟨{}| ", site_str);
    }

    if(status.algo_type == AlgorithmType::xDMRG and !std::isnan(status.energy_dens)) { report += fmt::format("e:{:<5.3f} ", status.energy_dens); }
    double ene = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy(tensors);
    double var = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy_variance(tensors);
    report += fmt::format("E:{:<20.16f} ", ene);
    report += fmt::format("σ²H:{:<8.2e} [{:<8.2e}] ", var, status.energy_variance_lowest);
    report += fmt::format("Sₑ({:>2}):{:<10.8f} ", tensors.state->get_position<long>(), tools::finite::measure::entanglement_entropy_current(*tensors.state));

    report += fmt::format("ε:{:<8.2e} ", tensors.state->get_truncation_error_active_max());
    // report += fmt::format("α:{:<8.2e} ", status.env_expansion_alpha);
    if(status.bond_max == status.bond_lim) {
        report += fmt::format("χ:{:<3}|", status.bond_max);
    } else {
        report += fmt::format("χ:{:<3}|{:<3}|", status.bond_max, status.bond_lim);
    }
    auto bonds_msites = std::clamp(settings::strategy::dmrg_min_blocksize - 1, 1ul, tensors.get_length<size_t>());
    auto bonds_maxims = std::vector<long>(bonds_msites, status.bond_max);
    auto bonds_merged = tools::finite::measure::bond_dimensions_active(*tensors.state);
    auto bonds_padlen = fmt::format("{}", bonds_maxims).size();
    auto bonds_string = fmt::format("{}", bonds_merged);
    report += fmt::format("{0:<{1}} ", bonds_string, bonds_padlen);

    if(last_optcost and last_optsolver and last_optalgo) {
        std::string short_optalgo, short_optcost;
        switch(last_optalgo.value()) {
            case OptAlgo::DIRECT: short_optalgo = "DIR"; break;
            case OptAlgo::DIRECTZ: short_optalgo = "DRZ"; break;
            case OptAlgo::DIRECTX2: short_optalgo = "DX2"; break;
            case OptAlgo::MPSEIGS: short_optalgo = "MPS"; break;
            case OptAlgo::SHIFTINV: short_optalgo = "SHI"; break;
            case OptAlgo::SUBSPACE: short_optalgo = "SUB"; break;
            default: short_optalgo = "???";
        }
        switch(last_optcost.value()) {
            case OptCost::OVERLAP: short_optcost = "OVE"; break;
            case OptCost::ENERGY: short_optcost = "ENE"; break;
            case OptCost::VARIANCE: short_optcost = "VAR"; break;
            default: short_optcost = "???";
        }
        report += fmt::format("opt:[{}|{}|{}] ", short_optcost, short_optalgo, enum2sv(last_optsolver.value()));
    }
    std::string stat = fmt::format("sat:{:<1}", status.algorithm_saturated_for);
    if(status.algorithm_has_stuck_for > 0) { stat = fmt::format("stk:{:<1} ", status.algorithm_has_stuck_for); }
    if(status.algorithm_converged_for > 0) { stat = fmt::format("con:{:<1} ", status.algorithm_converged_for); }
    report += stat;
    if(settings::precision::variance_saturation_sensitivity > 0 or settings::precision::entropy_saturation_sensitivity > 0) {
        std::vector<std::string> satstr;
        if(settings::precision::variance_saturation_sensitivity > 0) satstr.emplace_back(fmt::format("σ²H:{:<}", status.variance_mpo_saturated_for));
        if(settings::precision::entropy_saturation_sensitivity > 0) satstr.emplace_back(fmt::format("Sₑ:{:<1}", status.entanglement_saturated_for));
        if(settings::precision::infocom_saturation_sensitivity > 0) satstr.emplace_back(fmt::format("iₗ:{:<1}", infocom_saturated_for));
        report += fmt::format("[{}]", fmt::join(satstr, " "));
    } else {
        report += " ";
    }

    report += fmt::format("time:{:<9} ", fmt::format("{:>7.1f}s", tid::get_unscoped("t_tot").get_time()));
    report += fmt::format("mem[rss {:<.1f}|peak {:<.1f}|vm {:<.1f}]MB ", debug::mem_rss_in_mb(), debug::mem_hwm_in_mb(), debug::mem_vm_in_mb());
    tools::log->info(report);
}

void AlgorithmFinite::print_status_full() {
    tensors.clear_cache();
    tensors.clear_measurements();
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
        tools::log->info("Energy          E                  = {:<.16f}", energy);
        if(status.algo_type == AlgorithmType::xDMRG)
            tools::log->info("Energy density (rescaled 0 to 1) ε = {:<6.4f}",
                             tools::finite::measure::energy_normalized(tensors, status.energy_min, status.energy_max));
        double variance = tensors.active_sites.empty() ? std::numeric_limits<double>::quiet_NaN() : tools::finite::measure::energy_variance(tensors);
        tools::log->info("Energy variance σ²(H)              = {:<8.2e}", variance);
    }
    tools::log->info("Bond dimension maximum χmax        = {}", status.bond_max);
    tools::log->info("Bond dimensions χ                  = {}", tools::finite::measure::bond_dimensions(*tensors.state));
    tools::log->info("Bond dimension  χ (mid)            = {}", tools::finite::measure::bond_dimension_midchain(*tensors.state));
    tools::log->info("Entanglement entropies Sₑ          = {::8.2e}", tools::finite::measure::entanglement_entropies(*tensors.state));
    tools::log->info("Entanglement entropy   Sₑ (mid)    = {:8.2e}", tools::finite::measure::entanglement_entropy_midchain(*tensors.state));
    if(status.algo_type == AlgorithmType::fLBIT) {
        tools::log->info("Number entropies Sₙ                = {::8.2e}", tools::finite::measure::number_entropies(*tensors.state));
        tools::log->info("Number entropy   Sₙ (mid)          = {:8.2e}", tools::finite::measure::number_entropy_midchain(*tensors.state));
    }
    tools::log->info("Spin components (global X,Y,Z)     = {::.16f}", tools::finite::measure::spin_components(*tensors.state));

    if(status.algo_type == AlgorithmType::xDMRG) {
        auto expectation_values_xyz = tools::finite::measure::expectation_values_xyz(*tensors.state);
        auto structure_factor_xyz   = tools::finite::measure::structure_factor_xyz(*tensors.state);
        auto opdm_spectrum          = tools::finite::measure::opdm_spectrum(*tensors.state);
        tools::log->info("Expectation values ⟨σx⟩            = {::+9.6f}", tenx::span(expectation_values_xyz[0]));
        tools::log->info("Expectation values ⟨σy⟩            = {::+9.6f}", tenx::span(expectation_values_xyz[1]));
        tools::log->info("Expectation values ⟨σz⟩            = {::+9.6f}", tenx::span(expectation_values_xyz[2]));
        tools::log->info("Structure f. L⁻¹ ∑_ij ⟨σx_i σx_j⟩² = {:+.16f}", structure_factor_xyz[0]);
        tools::log->info("Structure f. L⁻¹ ∑_ij ⟨σy_i σy_j⟩² = {:+.16f}", structure_factor_xyz[1]);
        tools::log->info("Structure f. L⁻¹ ∑_ij ⟨σz_i σz_j⟩² = {:+.16f}", structure_factor_xyz[2]);
        tools::log->info("OPDM spectrum ⟨σ+..σz..σ-⟩         = {:.8f}", fmt::join(tenx::span(opdm_spectrum), ", "));
    }

    tools::log->info("Truncation Error limit             = {:8.2e}", status.trnc_lim);
    tools::log->info("Truncation Errors ε                = {::8.2e}", tensors.state->get_truncation_errors());
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
