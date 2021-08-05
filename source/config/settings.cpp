#include "settings.h"
#include <config/loader.h>

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
long settings::chi_lim_max(AlgorithmType algo_type) {
    switch(algo_type) {
        case AlgorithmType::iDMRG: return settings::idmrg::chi_lim_max;
        case AlgorithmType::fDMRG: return settings::fdmrg::chi_lim_max;
        case AlgorithmType::xDMRG: return settings::xdmrg::chi_lim_max;
        case AlgorithmType::iTEBD: return settings::itebd::chi_lim_max;
        case AlgorithmType::fLBIT: return settings::flbit::chi_lim_max;
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
ChiGrow settings::chi_lim_grow(AlgorithmType algo_type) {
    switch(algo_type) {
        case AlgorithmType::iDMRG: return settings::idmrg::chi_lim_grow;
        case AlgorithmType::fDMRG: return settings::fdmrg::chi_lim_grow;
        case AlgorithmType::xDMRG: return settings::xdmrg::chi_lim_grow;
        case AlgorithmType::iTEBD: return settings::itebd::chi_lim_grow;
        case AlgorithmType::fLBIT: return settings::flbit::chi_lim_grow;
        default: return ChiGrow::OFF;
    }
}
double settings::chi_lim_grow_factor(AlgorithmType algo_type) {
    switch(algo_type) {
        case AlgorithmType::iDMRG: return settings::idmrg::chi_lim_grow_factor;
        case AlgorithmType::fDMRG: return settings::fdmrg::chi_lim_grow_factor;
        case AlgorithmType::xDMRG: return settings::xdmrg::chi_lim_grow_factor;
        case AlgorithmType::iTEBD: return settings::itebd::chi_lim_grow_factor;
        case AlgorithmType::fLBIT: return settings::flbit::chi_lim_grow_factor;
        default: return 1.5;
    }
}

long settings::chi_lim_init(AlgorithmType algo_type) {
    switch(algo_type) {
        case AlgorithmType::iDMRG: return settings::idmrg::chi_lim_init;
        case AlgorithmType::fDMRG: return settings::fdmrg::chi_lim_init;
        case AlgorithmType::xDMRG: return settings::xdmrg::chi_lim_init;
        case AlgorithmType::iTEBD: return settings::itebd::chi_lim_init;
        case AlgorithmType::fLBIT: return settings::flbit::chi_lim_init;
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
    if(not dmrg_config.file_exists) throw std::runtime_error(fmt::format("Could not load config [{}]: File does not exist", dmrg_config.file_path.string()));
    dmrg_config.load();
    input::config_filename      = dmrg_config.file_path.string();
    input::config_file_contents = dmrg_config.get_config_file_as_string();
    /* clang-format off */
    dmrg_config.load_parameter("input::seed"                                  , input::seed);
    dmrg_config.load_parameter("input::bitfield"                              , input::bitfield);

    dmrg_config.load_parameter("output::output_filepath"                      , output::output_filepath);
    dmrg_config.load_parameter("output::save_profiling"                       , output::save_profiling);
    dmrg_config.load_parameter("output::savepoint_keep_newest_only"           , output::savepoint_keep_newest_only);
    dmrg_config.load_parameter("output::savepoint_frequency"                  , output::savepoint_frequency);
    dmrg_config.load_parameter("output::checkpoint_keep_newest_only"          , output::checkpoint_keep_newest_only);
    dmrg_config.load_parameter("output::checkpoint_keep_chi_updates"          , output::checkpoint_keep_chi_updates);
    dmrg_config.load_parameter("output::checkpoint_frequency"                 , output::checkpoint_frequency);
    dmrg_config.load_parameter("output::use_temp_dir"                         , output::use_temp_dir);
    dmrg_config.load_parameter("output::copy_from_temp_freq"                  , output::copy_from_temp_freq);
    dmrg_config.load_parameter("output::temp_dir"                             , output::temp_dir);
    dmrg_config.load_parameter("output::compression_level"                    , output::compression_level);
    dmrg_config.load_parameter("output::file_collision_policy"                , output::file_collision_policy);
    dmrg_config.load_parameter("output::file_resume_policy"                   , output::file_resume_policy);
    dmrg_config.load_parameter("output::storage_level_model"                  , output::storage_level_model);
    dmrg_config.load_parameter("output::storage_level_savepoint"              , output::storage_level_savepoint);
    dmrg_config.load_parameter("output::storage_level_checkpoint"             , output::storage_level_checkpoint);
    dmrg_config.load_parameter("output::storage_level_good_state"             , output::storage_level_good_state);
    dmrg_config.load_parameter("output::storage_level_fail_state"             , output::storage_level_fail_state);
    dmrg_config.load_parameter("output::storage_level_proj_state"             , output::storage_level_proj_state);
    dmrg_config.load_parameter("output::storage_level_init_state"             , output::storage_level_init_state);
    dmrg_config.load_parameter("output::storage_level_emin_state"             , output::storage_level_emin_state);
    dmrg_config.load_parameter("output::storage_level_emax_state"             , output::storage_level_emax_state);

    dmrg_config.load_parameter("model::model_type"                            , model::model_type);
    dmrg_config.load_parameter("model::model_size"                            , model::model_size);
    dmrg_config.load_parameter("model::ising_tf_rf::J1"                       , model::ising_tf_rf::J1);
    dmrg_config.load_parameter("model::ising_tf_rf::J2"                       , model::ising_tf_rf::J2);
    dmrg_config.load_parameter("model::ising_tf_rf::h_tran"                   , model::ising_tf_rf::h_tran);
    dmrg_config.load_parameter("model::ising_tf_rf::h_mean"                   , model::ising_tf_rf::h_mean);
    dmrg_config.load_parameter("model::ising_tf_rf::h_stdv"                   , model::ising_tf_rf::h_stdv);
    dmrg_config.load_parameter("model::ising_tf_rf::spin_dim"                 , model::ising_tf_rf::spin_dim);
    dmrg_config.load_parameter("model::ising_tf_rf::distribution"             , model::ising_tf_rf::distribution);
    dmrg_config.load_parameter("model::ising_sdual::delta"                    , model::ising_sdual::delta);
    dmrg_config.load_parameter("model::ising_sdual::lambda"                   , model::ising_sdual::lambda);
    dmrg_config.load_parameter("model::ising_sdual::J_stdv"                   , model::ising_sdual::J_stdv);
    dmrg_config.load_parameter("model::ising_sdual::h_stdv"                   , model::ising_sdual::h_stdv);
    dmrg_config.load_parameter("model::ising_sdual::parity_sep"               , model::ising_sdual::parity_sep);
    dmrg_config.load_parameter("model::ising_sdual::distribution"             , model::ising_sdual::distribution);
    dmrg_config.load_parameter("model::ising_sdual::spin_dim"                 , model::ising_sdual::spin_dim);
    dmrg_config.load_parameter("model::lbit::J1_mean"                         , model::lbit::J1_mean);
    dmrg_config.load_parameter("model::lbit::J2_mean"                         , model::lbit::J2_mean);
    dmrg_config.load_parameter("model::lbit::J3_mean"                         , model::lbit::J3_mean);
    dmrg_config.load_parameter("model::lbit::J1_wdth"                         , model::lbit::J1_wdth);
    dmrg_config.load_parameter("model::lbit::J2_wdth"                         , model::lbit::J2_wdth);
    dmrg_config.load_parameter("model::lbit::J3_wdth"                         , model::lbit::J3_wdth);
    dmrg_config.load_parameter("model::lbit::J2_base"                         , model::lbit::J2_base);
    dmrg_config.load_parameter("model::lbit::J2_span"                         , model::lbit::J2_span);

    dmrg_config.load_parameter("model::lbit::f_mixer"                         , model::lbit::f_mixer);
    dmrg_config.load_parameter("model::lbit::u_layer"                         , model::lbit::u_layer);
    dmrg_config.load_parameter("model::lbit::spin_dim"                        , model::lbit::spin_dim);
    dmrg_config.load_parameter("model::lbit::distribution"                    , model::lbit::distribution);

    dmrg_config.load_parameter("strategy::krylov_opt_when_stuck"              , strategy::krylov_opt_when_stuck);
    dmrg_config.load_parameter("strategy::chi_quench_when_stuck"              , strategy::chi_quench_when_stuck);
    dmrg_config.load_parameter("strategy::perturb_when_stuck"                 , strategy::perturb_when_stuck);
    dmrg_config.load_parameter("strategy::expand_subspace_when_stuck"         , strategy::expand_subspace_when_stuck);
    dmrg_config.load_parameter("strategy::expand_on_saturation"               , strategy::expand_on_saturation);
    dmrg_config.load_parameter("strategy::project_on_saturation"              , strategy::project_on_saturation);
    dmrg_config.load_parameter("strategy::project_on_every_iter"              , strategy::project_on_every_iter);
    dmrg_config.load_parameter("strategy::project_on_chi_update"              , strategy::project_on_chi_update);
    dmrg_config.load_parameter("strategy::project_initial_state"              , strategy::project_initial_state);
    dmrg_config.load_parameter("strategy::randomize_on_chi_update"            , strategy::randomize_on_chi_update);
    dmrg_config.load_parameter("strategy::randomize_early"                    , strategy::randomize_early);
    dmrg_config.load_parameter("strategy::use_eigenspinors"                   , strategy::use_eigenspinors);
    dmrg_config.load_parameter("strategy::max_resets"                         , strategy::max_resets);
    dmrg_config.load_parameter("strategy::max_stuck_iters"                    , strategy::max_stuck_iters);
    dmrg_config.load_parameter("strategy::max_saturation_iters"               , strategy::max_saturation_iters);
    dmrg_config.load_parameter("strategy::min_saturation_iters"               , strategy::min_saturation_iters);
    dmrg_config.load_parameter("strategy::min_converged_iters"                , strategy::min_converged_iters);
    dmrg_config.load_parameter("strategy::max_expansion_iters"                , strategy::max_expansion_iters);
    dmrg_config.load_parameter("strategy::multisite_mps_size_def"             , strategy::multisite_mps_size_def);
    dmrg_config.load_parameter("strategy::multisite_mps_size_max"             , strategy::multisite_mps_size_max);
    dmrg_config.load_parameter("strategy::multisite_mps_step"                 , strategy::multisite_mps_step);
    dmrg_config.load_parameter("strategy::target_sector"                      , strategy::target_sector);
    dmrg_config.load_parameter("strategy::initial_type"                       , strategy::initial_type);
    dmrg_config.load_parameter("strategy::initial_state"                      , strategy::initial_state);
    dmrg_config.load_parameter("strategy::secondary_states"                   , strategy::secondary_states);

    dmrg_config.load_parameter("precision::eig_max_iter"                      , precision::eig_max_iter);
    dmrg_config.load_parameter("precision::eig_tolerance"                     , precision::eig_tolerance);
    dmrg_config.load_parameter("precision::eig_default_ncv"                   , precision::eig_default_ncv);
    dmrg_config.load_parameter("precision::svd_threshold"                     , precision::svd_threshold);
    dmrg_config.load_parameter("precision::svd_switchsize"                    , precision::svd_switchsize);
    dmrg_config.load_parameter("precision::use_compressed_mpo_squared_all"    , precision::use_compressed_mpo_squared_all);
    dmrg_config.load_parameter("precision::use_compressed_mpo_squared_otf"    , precision::use_compressed_mpo_squared_otf);
    dmrg_config.load_parameter("precision::use_reduced_mpo_energy"            , precision::use_reduced_mpo_energy);
    dmrg_config.load_parameter("precision::variance_convergence_threshold"    , precision::variance_convergence_threshold);
    dmrg_config.load_parameter("precision::variance_saturation_sensitivity"   , precision::variance_saturation_sensitivity);
    dmrg_config.load_parameter("precision::entropy_saturation_sensitivity"    , precision::entropy_saturation_sensitivity);
    dmrg_config.load_parameter("precision::target_subspace_error"             , precision::target_subspace_error);
    dmrg_config.load_parameter("precision::max_subspace_size"                 , precision::max_subspace_size);
    dmrg_config.load_parameter("precision::max_size_full_diag"                , precision::max_size_full_diag);
    dmrg_config.load_parameter("precision::max_size_part_diag"                , precision::max_size_part_diag);
    dmrg_config.load_parameter("precision::max_size_direct"                   , precision::max_size_direct);
    dmrg_config.load_parameter("precision::max_norm_error"                    , precision::max_norm_error);

    dmrg_config.load_parameter("threading::omp_threads"                       , threading::omp_threads);
    dmrg_config.load_parameter("threading::stl_threads"                       , threading::stl_threads);

    //Parameters controlling finite-DMRG
    dmrg_config.load_parameter("fdmrg::on"                      , fdmrg::on);
    dmrg_config.load_parameter("fdmrg::max_iters"               , fdmrg::max_iters);
    dmrg_config.load_parameter("fdmrg::min_iters"               , fdmrg::min_iters);
    dmrg_config.load_parameter("fdmrg::chi_lim_max"             , fdmrg::chi_lim_max);
    dmrg_config.load_parameter("fdmrg::chi_lim_grow"            , fdmrg::chi_lim_grow);
    dmrg_config.load_parameter("fdmrg::chi_lim_grow_factor"     , fdmrg::chi_lim_grow_factor);
    dmrg_config.load_parameter("fdmrg::chi_lim_init"            , fdmrg::chi_lim_init);
    dmrg_config.load_parameter("fdmrg::print_freq "             , fdmrg::print_freq);
    dmrg_config.load_parameter("fdmrg::store_wavefn"            , fdmrg::store_wavefn);

    //Parameters controlling finite-LBIT
    dmrg_config.load_parameter("flbit::on"                      , flbit::on);
    dmrg_config.load_parameter("flbit::max_iters"               , flbit::max_iters);
    dmrg_config.load_parameter("flbit::min_iters"               , flbit::min_iters);
    dmrg_config.load_parameter("flbit::chi_lim_max"             , flbit::chi_lim_max);
    dmrg_config.load_parameter("flbit::chi_lim_grow"            , flbit::chi_lim_grow);
    dmrg_config.load_parameter("flbit::chi_lim_grow_factor"     , flbit::chi_lim_grow_factor);
    dmrg_config.load_parameter("flbit::chi_lim_init"            , flbit::chi_lim_init);
    dmrg_config.load_parameter("flbit::time_start_real"         , flbit::time_start_real);
    dmrg_config.load_parameter("flbit::time_start_imag"         , flbit::time_start_imag);
    dmrg_config.load_parameter("flbit::time_final_real"         , flbit::time_final_real);
    dmrg_config.load_parameter("flbit::time_final_imag"         , flbit::time_final_imag);
    dmrg_config.load_parameter("flbit::time_num_steps"          , flbit::time_num_steps);
    dmrg_config.load_parameter("flbit::print_freq"              , flbit::print_freq);
    dmrg_config.load_parameter("flbit::compute_lbit_length"     , flbit::compute_lbit_length);
    dmrg_config.load_parameter("flbit::compute_lbit_stats"      , flbit::compute_lbit_stats);
    dmrg_config.load_parameter("flbit::store_wavefn"            , flbit::store_wavefn);

    //Parameters controlling excited state DMRG
    dmrg_config.load_parameter("xdmrg::on"                           , xdmrg::on);
    dmrg_config.load_parameter("xdmrg::max_iters"                    , xdmrg::max_iters);
    dmrg_config.load_parameter("xdmrg::min_iters"                    , xdmrg::min_iters);
    dmrg_config.load_parameter("xdmrg::olap_iters"                   , xdmrg::olap_iters);
    dmrg_config.load_parameter("xdmrg::vsub_iters"                   , xdmrg::vsub_iters);
    dmrg_config.load_parameter("xdmrg::chi_lim_max"                  , xdmrg::chi_lim_max);
    dmrg_config.load_parameter("xdmrg::chi_lim_grow"                 , xdmrg::chi_lim_grow);
    dmrg_config.load_parameter("xdmrg::chi_lim_grow_factor"          , xdmrg::chi_lim_grow_factor);
    dmrg_config.load_parameter("xdmrg::chi_lim_init"                 , xdmrg::chi_lim_init);
    dmrg_config.load_parameter("xdmrg::chi_lim_olap"                 , xdmrg::chi_lim_olap);
    dmrg_config.load_parameter("xdmrg::chi_lim_vsub"                 , xdmrg::chi_lim_vsub);
    dmrg_config.load_parameter("xdmrg::print_freq "                  , xdmrg::print_freq);
    dmrg_config.load_parameter("xdmrg::store_wavefn"                 , xdmrg::store_wavefn);
    dmrg_config.load_parameter("xdmrg::energy_density_target"        , xdmrg::energy_density_target);
    dmrg_config.load_parameter("xdmrg::energy_density_window"        , xdmrg::energy_density_window);
    dmrg_config.load_parameter("xdmrg::max_states"                   , xdmrg::max_states);
    dmrg_config.load_parameter("xdmrg::finish_if_entanglm_saturated" , xdmrg::finish_if_entanglm_saturated);
    dmrg_config.load_parameter("xdmrg::finish_if_variance_saturated" , xdmrg::finish_if_variance_saturated);

    //Parameters controlling infinite-DMRG
    dmrg_config.load_parameter("idmrg::on"                  , idmrg::on);
    dmrg_config.load_parameter("idmrg::max_iters"           , idmrg::max_iters);
    dmrg_config.load_parameter("idmrg::chi_lim_max"         , idmrg::chi_lim_max);
    dmrg_config.load_parameter("idmrg::chi_lim_grow"        , idmrg::chi_lim_grow);
    dmrg_config.load_parameter("idmrg::chi_lim_grow_factor" , idmrg::chi_lim_grow_factor);
    dmrg_config.load_parameter("idmrg::chi_lim_init"        , idmrg::chi_lim_init);
    dmrg_config.load_parameter("idmrg::print_freq"          , idmrg::print_freq);


    //Parameters controlling imaginary TEBD (Zero temperature)
    dmrg_config.load_parameter("itebd::on"                  , itebd::on       );
    dmrg_config.load_parameter("itebd::max_iters"           , itebd::max_iters);
    dmrg_config.load_parameter("itebd::time_step_init_real" , itebd::time_step_init_real  );
    dmrg_config.load_parameter("itebd::time_step_init_imag" , itebd::time_step_init_imag  );
    dmrg_config.load_parameter("itebd::time_step_min"       , itebd::time_step_min);
    dmrg_config.load_parameter("itebd::suzuki_order"        , itebd::suzuki_order);
    dmrg_config.load_parameter("itebd::chi_lim_max"         , itebd::chi_lim_max  );
    dmrg_config.load_parameter("itebd::chi_lim_grow"        , itebd::chi_lim_grow);
    dmrg_config.load_parameter("itebd::chi_lim_grow_factor" , itebd::chi_lim_grow_factor);
    dmrg_config.load_parameter("itebd::chi_lim_init"        , itebd::chi_lim_init);
    dmrg_config.load_parameter("itebd::print_freq"          , itebd::print_freq);

    //Profiling
    dmrg_config.load_parameter("profiling::on"        , profiling::on        );
    dmrg_config.load_parameter("profiling::precision" , profiling::precision );
    dmrg_config.load_parameter("profiling::extra"     , profiling::extra     );
    //Console settings
    dmrg_config.load_parameter("console::verbosity"   , console::verbosity);
    dmrg_config.load_parameter("console::timestamp"   , console::timestamp);
    /* clang-format off */
}



void settings::load(std::string_view  config_filename){
    Loader indata(config_filename);
    settings::load(indata);
}
