#include "settings.h"
#include "debug/exceptions.h"
#include "loader.h"

bool settings::algorithm_is_on(AlgorithmType algo_type) {
    switch(algo_type) {
        case AlgorithmType::iDMRG: return settings::idmrg::on;
        case AlgorithmType::fDMRG: return settings::fdmrg::on;
        case AlgorithmType::xDMRG: return settings::xdmrg::on;
        case AlgorithmType::iTEBD: return settings::itebd::on;
        case AlgorithmType::fLBIT: return settings::flbit::on;
        default: return false;
    }
}
long settings::get_bond_max(AlgorithmType algo_type) {
    switch(algo_type) {
        case AlgorithmType::iDMRG: return settings::idmrg::bond_max;
        case AlgorithmType::fDMRG: return settings::fdmrg::bond_max;
        case AlgorithmType::xDMRG: return settings::xdmrg::bond_max;
        case AlgorithmType::iTEBD: return settings::itebd::bond_max;
        case AlgorithmType::fLBIT: return settings::flbit::bond_max;
        default: return 64;
    }
}
size_t settings::print_freq(AlgorithmType algo_type) {
    switch(algo_type) {
        case AlgorithmType::iDMRG: return settings::idmrg::print_freq;
        case AlgorithmType::fDMRG: return settings::fdmrg::print_freq;
        case AlgorithmType::xDMRG: return settings::xdmrg::print_freq;
        case AlgorithmType::iTEBD: return settings::itebd::print_freq;
        case AlgorithmType::fLBIT: return settings::flbit::print_freq;
        default: return 1;
    }
}

long settings::get_bond_init(AlgorithmType algo_type) {
    switch(algo_type) {
        case AlgorithmType::iDMRG: return settings::idmrg::bond_init;
        case AlgorithmType::fDMRG: return settings::fdmrg::bond_init;
        case AlgorithmType::xDMRG: return settings::xdmrg::bond_init;
        case AlgorithmType::iTEBD: return settings::itebd::bond_init;
        case AlgorithmType::fLBIT: return settings::flbit::bond_init;
        default: return 8;
    }
}
bool settings::store_wave_function(AlgorithmType algo_type) {
    switch(algo_type) {
        case AlgorithmType::iDMRG: return false;
        case AlgorithmType::fDMRG: return settings::fdmrg::store_wavefn;
        case AlgorithmType::xDMRG: return settings::xdmrg::store_wavefn;
        case AlgorithmType::iTEBD: return false;
        case AlgorithmType::fLBIT: return settings::flbit::store_wavefn;
        default: return false;
    }
}

/*
 * The section below attempts to find and read the parameters from the given inputfile.
 *
 */

void settings::load(Loader &dmrg_config) {
    if(not dmrg_config.file_exists) throw except::runtime_error("Could not load config [{}]: File does not exist", dmrg_config.file_path.string());
    dmrg_config.load();
    input::config_filename      = dmrg_config.file_path.string();
    input::config_file_contents = dmrg_config.get_config_file_as_string();
    /* clang-format off */
    dmrg_config.load_parameter("threading::omp_threads"                       , threading::omp_threads);
    dmrg_config.load_parameter("threading::stl_threads"                       , threading::stl_threads);

    dmrg_config.load_parameter("input::seed"                                  , input::seed);
    dmrg_config.load_parameter("input::bitfield"                              , input::bitfield);

    dmrg_config.load_parameter("storage::output_filepath"                     , storage::output_filepath);
    dmrg_config.load_parameter("storage::save_timers"                         , storage::save_timers);
    dmrg_config.load_parameter("storage::savepoint_keep_newest_only"          , storage::savepoint_keep_newest_only);
    dmrg_config.load_parameter("storage::savepoint_frequency"                 , storage::savepoint_frequency);
    dmrg_config.load_parameter("storage::checkpoint_keep_newest_only"         , storage::checkpoint_keep_newest_only);
    dmrg_config.load_parameter("storage::checkpoint_frequency"                , storage::checkpoint_frequency);
    dmrg_config.load_parameter("storage::use_temp_dir"                        , storage::use_temp_dir);
    dmrg_config.load_parameter("storage::copy_from_temp_freq"                 , storage::copy_from_temp_freq);
    dmrg_config.load_parameter("storage::temp_dir"                            , storage::temp_dir);
    dmrg_config.load_parameter("storage::compression_level"                   , storage::compression_level);
    dmrg_config.load_parameter("storage::file_collision_policy"               , storage::file_collision_policy);
    dmrg_config.load_parameter("storage::file_resume_policy"                  , storage::file_resume_policy);
    dmrg_config.load_parameter("storage::file_resume_name"                    , storage::file_resume_name);
    dmrg_config.load_parameter("storage::file_resume_iter"                    , storage::file_resume_iter);
    dmrg_config.load_parameter("storage::storage_level_model"                 , storage::storage_level_model);
    dmrg_config.load_parameter("storage::storage_level_savepoint"             , storage::storage_level_savepoint);
    dmrg_config.load_parameter("storage::storage_level_checkpoint"            , storage::storage_level_checkpoint);
    dmrg_config.load_parameter("storage::storage_level_finished"              , storage::storage_level_finished);
    dmrg_config.load_parameter("storage::storage_level_proj_state"            , storage::storage_level_proj_state);
    dmrg_config.load_parameter("storage::storage_level_init_state"            , storage::storage_level_init_state);
    dmrg_config.load_parameter("storage::storage_level_emin_state"            , storage::storage_level_emin_state);
    dmrg_config.load_parameter("storage::storage_level_emax_state"            , storage::storage_level_emax_state);
    dmrg_config.load_parameter("storage::storage_level_bond_state"            , storage::storage_level_bond_state);
    dmrg_config.load_parameter("storage::storage_level_trnc_state"            , storage::storage_level_trnc_state);
    dmrg_config.load_parameter("storage::storage_level_fes_state"             , storage::storage_level_fes_state);

    dmrg_config.load_parameter("model::model_type"                            , model::model_type);
    dmrg_config.load_parameter("model::model_size"                            , model::model_size);
    dmrg_config.load_parameter("model::ising_tf_rf::J1"                       , model::ising_tf_rf::J1);
    dmrg_config.load_parameter("model::ising_tf_rf::J2"                       , model::ising_tf_rf::J2);
    dmrg_config.load_parameter("model::ising_tf_rf::h_tran"                   , model::ising_tf_rf::h_tran);
    dmrg_config.load_parameter("model::ising_tf_rf::h_mean"                   , model::ising_tf_rf::h_mean);
    dmrg_config.load_parameter("model::ising_tf_rf::h_wdth"                   , model::ising_tf_rf::h_wdth);
    dmrg_config.load_parameter("model::ising_tf_rf::spin_dim"                 , model::ising_tf_rf::spin_dim);
    dmrg_config.load_parameter("model::ising_tf_rf::distribution"             , model::ising_tf_rf::distribution);
    dmrg_config.load_parameter("model::ising_sdual::delta"                    , model::ising_sdual::delta);
    dmrg_config.load_parameter("model::ising_sdual::lambda"                   , model::ising_sdual::lambda);
    dmrg_config.load_parameter("model::ising_sdual::parity_sep"               , model::ising_sdual::parity_sep);
    dmrg_config.load_parameter("model::ising_sdual::distribution"             , model::ising_sdual::distribution);
    dmrg_config.load_parameter("model::ising_sdual::spin_dim"                 , model::ising_sdual::spin_dim);
    dmrg_config.load_parameter("model::ising_majorana::delta"                 , model::ising_majorana::delta);
    dmrg_config.load_parameter("model::ising_majorana::g"                     , model::ising_majorana::g);
    dmrg_config.load_parameter("model::ising_majorana::parity_sep"            , model::ising_majorana::parity_sep);
    dmrg_config.load_parameter("model::ising_majorana::distribution"          , model::ising_majorana::distribution);
    dmrg_config.load_parameter("model::ising_majorana::spin_dim"              , model::ising_majorana::spin_dim);
    dmrg_config.load_parameter("model::lbit::J1_mean"                         , model::lbit::J1_mean);
    dmrg_config.load_parameter("model::lbit::J2_mean"                         , model::lbit::J2_mean);
    dmrg_config.load_parameter("model::lbit::J3_mean"                         , model::lbit::J3_mean);
    dmrg_config.load_parameter("model::lbit::J1_wdth"                         , model::lbit::J1_wdth);
    dmrg_config.load_parameter("model::lbit::J2_wdth"                         , model::lbit::J2_wdth);
    dmrg_config.load_parameter("model::lbit::J3_wdth"                         , model::lbit::J3_wdth);
    dmrg_config.load_parameter("model::lbit::J2_xcls"                         , model::lbit::J2_xcls);
    dmrg_config.load_parameter("model::lbit::J2_span"                         , model::lbit::J2_span);

    dmrg_config.load_parameter("model::lbit::f_mixer"                         , model::lbit::f_mixer);
    dmrg_config.load_parameter("model::lbit::u_layer"                         , model::lbit::u_layer);
    dmrg_config.load_parameter("model::lbit::spin_dim"                        , model::lbit::spin_dim);
    dmrg_config.load_parameter("model::lbit::distribution"                    , model::lbit::distribution);

    dmrg_config.load_parameter("strategy::move_sites_when_stuck"              , strategy::move_sites_when_stuck);
    dmrg_config.load_parameter("strategy::project_on_saturation"              , strategy::project_on_saturation);
    dmrg_config.load_parameter("strategy::project_on_every_iter"              , strategy::project_on_every_iter);
    dmrg_config.load_parameter("strategy::project_on_bond_update"             , strategy::project_on_bond_update);
    dmrg_config.load_parameter("strategy::project_initial_state"              , strategy::project_initial_state);
    dmrg_config.load_parameter("strategy::project_final_state"                , strategy::project_final_state);
    dmrg_config.load_parameter("strategy::randomize_on_bond_update"           , strategy::randomize_on_bond_update);
    dmrg_config.load_parameter("strategy::randomize_early"                    , strategy::randomize_early);
    dmrg_config.load_parameter("strategy::use_eigenspinors"                   , strategy::use_eigenspinors);
    dmrg_config.load_parameter("strategy::max_resets"                         , strategy::max_resets);
    dmrg_config.load_parameter("strategy::max_stuck_iters"                    , strategy::max_stuck_iters);
    dmrg_config.load_parameter("strategy::max_saturation_iters"               , strategy::max_saturation_iters);
    dmrg_config.load_parameter("strategy::min_saturation_iters"               , strategy::min_saturation_iters);
    dmrg_config.load_parameter("strategy::min_converged_iters"                , strategy::min_converged_iters);
    dmrg_config.load_parameter("strategy::max_env_expansion_alpha"            , strategy::max_env_expansion_alpha);
    dmrg_config.load_parameter("strategy::multisite_mps_site_def"             , strategy::multisite_mps_site_def);
    dmrg_config.load_parameter("strategy::multisite_mps_site_max"             , strategy::multisite_mps_site_max);
    dmrg_config.load_parameter("strategy::multisite_mps_move"                 , strategy::multisite_mps_move);
    dmrg_config.load_parameter("strategy::multisite_mps_when"                 , strategy::multisite_mps_when);
    dmrg_config.load_parameter("strategy::target_axis"                        , strategy::target_axis);
    dmrg_config.load_parameter("strategy::initial_axis"                       , strategy::initial_axis);
    dmrg_config.load_parameter("strategy::initial_type"                       , strategy::initial_type);
    dmrg_config.load_parameter("strategy::initial_state"                      , strategy::initial_state);
    dmrg_config.load_parameter("strategy::secondary_states"                   , strategy::secondary_states);
    dmrg_config.load_parameter("strategy::fes_rate"                           , strategy::fes_rate);
    dmrg_config.load_parameter("strategy::bond_increase_when"                 , strategy::bond_increase_when);
    dmrg_config.load_parameter("strategy::bond_increase_rate"                 , strategy::bond_increase_rate);
    dmrg_config.load_parameter("strategy::trnc_decrease_when"                 , strategy::trnc_decrease_when);
    dmrg_config.load_parameter("strategy::trnc_decrease_rate"                 , strategy::trnc_decrease_rate);

    dmrg_config.load_parameter("solver::eigs_max_iter"                        , solver::eigs_max_iter);
    dmrg_config.load_parameter("solver::eigs_tolerance"                       , solver::eigs_tolerance);
    dmrg_config.load_parameter("solver::eigs_default_ncv"                     , solver::eigs_default_ncv);
    dmrg_config.load_parameter("solver::bfgs_max_iter"                        , solver::bfgs_max_iter);
    dmrg_config.load_parameter("solver::iter_stuck_multiplier"                , solver::iter_stuck_multiplier);
    dmrg_config.load_parameter("solver::max_size_full_eigs"                   , solver::max_size_full_eigs);
    dmrg_config.load_parameter("solver::max_size_shift_invert"                , solver::max_size_shift_invert);
    dmrg_config.load_parameter("solver::bfgs_fix_rnorm_w_eigs"                , solver::bfgs_fix_rnorm_w_eigs);
    dmrg_config.load_parameter("solver::prefer_eigs_over_bfgs"                , solver::prefer_eigs_over_bfgs);
    dmrg_config.load_parameter("solver::max_grad_tolerance"                   , solver::max_grad_tolerance);
    dmrg_config.load_parameter("solver::svd_truncation_lim"                   , solver::svd_truncation_lim);
    dmrg_config.load_parameter("solver::svd_truncation_init"                  , solver::svd_truncation_init);
    dmrg_config.load_parameter("solver::svd_switchsize_bdc"                   , solver::svd_switchsize_bdc);

    dmrg_config.load_parameter("precision::use_compressed_mpo_squared_all"    , precision::use_compressed_mpo_squared_all);
    dmrg_config.load_parameter("precision::use_compressed_mpo_squared_otf"    , precision::use_compressed_mpo_squared_otf);
    dmrg_config.load_parameter("precision::use_mpo_energy_shift"              , precision::use_mpo_energy_shift);
    dmrg_config.load_parameter("precision::use_projection_on_mpo_squared"     , precision::use_projection_on_mpo_squared);
    dmrg_config.load_parameter("precision::variance_convergence_threshold"    , precision::variance_convergence_threshold);
    dmrg_config.load_parameter("precision::variance_saturation_sensitivity"   , precision::variance_saturation_sensitivity);
    dmrg_config.load_parameter("precision::entropy_saturation_sensitivity"    , precision::entropy_saturation_sensitivity);
    dmrg_config.load_parameter("precision::target_subspace_error"             , precision::target_subspace_error);
    dmrg_config.load_parameter("precision::max_subspace_size"                 , precision::max_subspace_size);
    dmrg_config.load_parameter("precision::max_size_multisite"                , precision::max_size_multisite);
    dmrg_config.load_parameter("precision::max_norm_error"                    , precision::max_norm_error);

    //Parameters controlling finite-DMRG
    dmrg_config.load_parameter("fdmrg::on"                      , fdmrg::on);
    dmrg_config.load_parameter("fdmrg::max_iters"               , fdmrg::max_iters);
    dmrg_config.load_parameter("fdmrg::min_iters"               , fdmrg::min_iters);
    dmrg_config.load_parameter("fdmrg::bond_max"                , fdmrg::bond_max);
    dmrg_config.load_parameter("fdmrg::bond_init"               , fdmrg::bond_init);
    dmrg_config.load_parameter("fdmrg::print_freq "             , fdmrg::print_freq);
    dmrg_config.load_parameter("fdmrg::store_wavefn"            , fdmrg::store_wavefn);

    //Parameters controlling finite-LBIT
    dmrg_config.load_parameter("flbit::on"                      , flbit::on);
    dmrg_config.load_parameter("flbit::max_iters"               , flbit::max_iters);
    dmrg_config.load_parameter("flbit::min_iters"               , flbit::min_iters);
    dmrg_config.load_parameter("flbit::use_swap_gates"          , flbit::use_swap_gates);
    dmrg_config.load_parameter("flbit::bond_max"                , flbit::bond_max);
    dmrg_config.load_parameter("flbit::bond_init"               , flbit::bond_init);
    dmrg_config.load_parameter("flbit::time_start_real"         , flbit::time_start_real);
    dmrg_config.load_parameter("flbit::time_start_imag"         , flbit::time_start_imag);
    dmrg_config.load_parameter("flbit::time_final_real"         , flbit::time_final_real);
    dmrg_config.load_parameter("flbit::time_final_imag"         , flbit::time_final_imag);
    dmrg_config.load_parameter("flbit::time_num_steps"          , flbit::time_num_steps);
    dmrg_config.load_parameter("flbit::time_gate_id_threshold"  , flbit::time_gate_id_threshold);
    dmrg_config.load_parameter("flbit::print_freq"              , flbit::print_freq);
    dmrg_config.load_parameter("flbit::compute_lbit_length"     , flbit::compute_lbit_length);
    dmrg_config.load_parameter("flbit::compute_lbit_stats"      , flbit::compute_lbit_stats);
    dmrg_config.load_parameter("flbit::store_wavefn"            , flbit::store_wavefn);
    dmrg_config.load_parameter("flbit::save_swap_gates"         , flbit::save_swap_gates);

    //Parameters controlling excited state DMRG
    dmrg_config.load_parameter("xdmrg::on"                           , xdmrg::on);
    dmrg_config.load_parameter("xdmrg::max_iters"                    , xdmrg::max_iters);
    dmrg_config.load_parameter("xdmrg::min_iters"                    , xdmrg::min_iters);
    dmrg_config.load_parameter("xdmrg::opt_overlap_iters"            , xdmrg::opt_overlap_iters);
    dmrg_config.load_parameter("xdmrg::opt_overlap_bond_lim"         , xdmrg::opt_overlap_bond_lim);
    dmrg_config.load_parameter("xdmrg::opt_subspace_iters"           , xdmrg::opt_subspace_iters);
    dmrg_config.load_parameter("xdmrg::opt_subspace_bond_lim"        , xdmrg::opt_subspace_bond_lim);
    dmrg_config.load_parameter("xdmrg::bond_max"                     , xdmrg::bond_max);
    dmrg_config.load_parameter("xdmrg::bond_init"                    , xdmrg::bond_init);

    dmrg_config.load_parameter("xdmrg::print_freq "                  , xdmrg::print_freq);
    dmrg_config.load_parameter("xdmrg::store_wavefn"                 , xdmrg::store_wavefn);
    dmrg_config.load_parameter("xdmrg::energy_density_target"        , xdmrg::energy_density_target);
    dmrg_config.load_parameter("xdmrg::energy_density_window"        , xdmrg::energy_density_window);
    dmrg_config.load_parameter("xdmrg::max_states"                   , xdmrg::max_states);
    dmrg_config.load_parameter("xdmrg::finish_if_entanglm_saturated" , xdmrg::finish_if_entanglm_saturated);
    dmrg_config.load_parameter("xdmrg::finish_if_variance_saturated" , xdmrg::finish_if_variance_saturated);
    //Parameters controlling infinite-DMRG
    dmrg_config.load_parameter("idmrg::on"                           , idmrg::on);
    dmrg_config.load_parameter("idmrg::max_iters"                    , idmrg::max_iters);
    dmrg_config.load_parameter("idmrg::bond_max"                     , idmrg::bond_max);
    dmrg_config.load_parameter("idmrg::bond_init"                    , idmrg::bond_init);
    dmrg_config.load_parameter("idmrg::print_freq"                   , idmrg::print_freq);


    //Parameters controlling imaginary TEBD (Zero temperature)
    dmrg_config.load_parameter("itebd::on"                   , itebd::on       );
    dmrg_config.load_parameter("itebd::max_iters"            , itebd::max_iters);
    dmrg_config.load_parameter("itebd::time_step_init_real"  , itebd::time_step_init_real  );
    dmrg_config.load_parameter("itebd::time_step_init_imag"  , itebd::time_step_init_imag  );
    dmrg_config.load_parameter("itebd::time_step_min"        , itebd::time_step_min);
    dmrg_config.load_parameter("itebd::suzuki_order"         , itebd::suzuki_order);
    dmrg_config.load_parameter("itebd::bond_max"             , itebd::bond_max  );
    dmrg_config.load_parameter("itebd::bond_init"            , itebd::bond_init);
    dmrg_config.load_parameter("itebd::print_freq"           , itebd::print_freq);

    //Timers
    dmrg_config.load_parameter("timer::on"         , timer::on        );
    dmrg_config.load_parameter("timer::precision"  , timer::precision );
    dmrg_config.load_parameter("timer::level"      , timer::level     );
    //Console settings
    dmrg_config.load_parameter("console::loglevel"    , console::loglevel);
    dmrg_config.load_parameter("console::logh5pp"     , console::logh5pp);
    dmrg_config.load_parameter("console::timestamp"   , console::timestamp);
    /* clang-format off */
}



void settings::load(std::string_view  config_filename){
    Loader indata(config_filename);
    settings::load(indata);
}
