//
// Created by david on 2018-01-14.
//

#include "nmspc_settings.h"
#include <h5pp/h5pp.h>
#include <io/class_config_reader.h>
#include <config/class_dmrg_config.h>
using namespace std;

/*
 * The section below attempts to find and read the parameters from the given inputfile.
 *
 */

void settings::load_config(class_dmrg_config &dmrg_config) {
    if(not dmrg_config.file_exists) throw std::runtime_error(fmt::format("Could not load config [{}]: File does not exist", dmrg_config.file_path.string()));
    dmrg_config.load();
    input::config_filename      = dmrg_config.file_path.string();
    input::config_file_contents = dmrg_config.get_config_file_as_string();
    /* clang-format off */
    dmrg_config.load_parameter("input::seed"                                  , input::seed);
    dmrg_config.load_parameter("input::bitfield"                              , input::bitfield);
    dmrg_config.load_parameter("output::output_filepath"                      , output::output_filepath);
    dmrg_config.load_parameter("output::save_profiling"                       , output::save_profiling);
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
    dmrg_config.load_parameter("strategy::compress_mpo_squared"               , strategy::compress_mpo_squared);
    dmrg_config.load_parameter("strategy::chi_quench_when_stuck"              , strategy::chi_quench_when_stuck);
    dmrg_config.load_parameter("strategy::perturb_when_stuck"                 , strategy::perturb_when_stuck);
    dmrg_config.load_parameter("strategy::damping_when_stuck"                 , strategy::damping_when_stuck);
    dmrg_config.load_parameter("strategy::discard_schmidt_when_stuck"         , strategy::discard_schmidt_when_stuck);
    dmrg_config.load_parameter("strategy::project_trial_when_stuck"           , strategy::project_when_stuck);
    dmrg_config.load_parameter("strategy::project_on_every_iter"              , strategy::project_on_every_iter);
    dmrg_config.load_parameter("strategy::project_on_chi_update"              , strategy::project_on_chi_update);
    dmrg_config.load_parameter("strategy::randomize_on_chi_update"            , strategy::randomize_on_chi_update);
    dmrg_config.load_parameter("strategy::randomize_early"                    , strategy::randomize_early);
    dmrg_config.load_parameter("strategy::use_eigenspinors"                   , strategy::use_eigenspinors);
    dmrg_config.load_parameter("strategy::max_resets"                         , strategy::max_resets);
    dmrg_config.load_parameter("strategy::multisite_max_sites"                , strategy::multisite_max_sites);
    dmrg_config.load_parameter("strategy::multisite_move"                     , strategy::multisite_move);
    dmrg_config.load_parameter("strategy::target_sector"                      , strategy::target_sector);
    dmrg_config.load_parameter("strategy::initial_type"                       , strategy::initial_type);
    dmrg_config.load_parameter("strategy::initial_state"                      , strategy::initial_state);
    dmrg_config.load_parameter("strategy::secondary_states"                   , strategy::secondary_states);
    dmrg_config.load_parameter("precision::eig_max_iter"                      , precision::eig_max_iter);
    dmrg_config.load_parameter("precision::eig_threshold"                     , precision::eig_threshold);
    dmrg_config.load_parameter("precision::eig_max_ncv"                       , precision::eig_max_ncv);
    dmrg_config.load_parameter("precision::svd_threshold"                     , precision::svd_threshold);
    dmrg_config.load_parameter("precision::svd_switchsize"                    , precision::svd_switchsize);
    dmrg_config.load_parameter("precision::variance_convergence_threshold"    , precision::variance_convergence_threshold);
    dmrg_config.load_parameter("precision::variance_slope_threshold"          , precision::variance_slope_threshold);
    dmrg_config.load_parameter("precision::entropy_slope_threshold"           , precision::entropy_slope_threshold);
    dmrg_config.load_parameter("precision::subspace_error_factor"             , precision::subspace_error_factor);
    dmrg_config.load_parameter("precision::max_subspace_error"                , precision::max_subspace_error);
    dmrg_config.load_parameter("precision::min_subspace_error"                , precision::min_subspace_error);
    dmrg_config.load_parameter("precision::max_size_full_diag"                , precision::max_size_full_diag);
    dmrg_config.load_parameter("precision::max_size_part_diag"                , precision::max_size_part_diag);
    dmrg_config.load_parameter("precision::max_size_direct"                   , precision::max_size_direct);
    dmrg_config.load_parameter("precision::max_norm_error"                    , precision::max_norm_error);
    dmrg_config.load_parameter("precision::use_reduced_energy"                , precision::use_reduced_energy);
    dmrg_config.load_parameter("threading::num_threads"                       , threading::num_threads);

    //Parameters controlling finite-DMRG
    dmrg_config.load_parameter("fdmrg::on"           , fdmrg::on);
    dmrg_config.load_parameter("fdmrg::max_iters"    , fdmrg::max_iters);
    dmrg_config.load_parameter("fdmrg::min_iters"    , fdmrg::min_iters);
    dmrg_config.load_parameter("fdmrg::chi_lim_max"  , fdmrg::chi_lim_max);
    dmrg_config.load_parameter("fdmrg::chi_lim_grow" , fdmrg::chi_lim_grow);
    dmrg_config.load_parameter("fdmrg::chi_lim_init" , fdmrg::chi_lim_init);
    dmrg_config.load_parameter("fdmrg::print_freq "  , fdmrg::print_freq);
    dmrg_config.load_parameter("fdmrg::store_wavefn" , fdmrg::store_wavefn);

    //Parameters controlling finite-LBIT
    dmrg_config.load_parameter("flbit::on"           , flbit::on);
    dmrg_config.load_parameter("flbit::max_iters"    , flbit::max_iters);
    dmrg_config.load_parameter("flbit::min_iters"    , flbit::min_iters);
    dmrg_config.load_parameter("flbit::chi_lim_max"  , flbit::chi_lim_max);
    dmrg_config.load_parameter("flbit::chi_lim_grow" , flbit::chi_lim_grow);
    dmrg_config.load_parameter("flbit::chi_lim_init" , flbit::chi_lim_init);
    dmrg_config.load_parameter("flbit::print_freq "  , flbit::print_freq);
    dmrg_config.load_parameter("flbit::store_wavefn" , flbit::store_wavefn);

    //Parameters controlling excited state DMRG
    dmrg_config.load_parameter("xdmrg::on"                     , xdmrg::on);
    dmrg_config.load_parameter("xdmrg::max_iters"              , xdmrg::max_iters);
    dmrg_config.load_parameter("xdmrg::min_iters"              , xdmrg::min_iters);
    dmrg_config.load_parameter("xdmrg::olap_iters"             , xdmrg::olap_iters);
    dmrg_config.load_parameter("xdmrg::vsub_iters"             , xdmrg::vsub_iters);
    dmrg_config.load_parameter("xdmrg::chi_lim_max"            , xdmrg::chi_lim_max);
    dmrg_config.load_parameter("xdmrg::chi_lim_grow"           , xdmrg::chi_lim_grow);
    dmrg_config.load_parameter("xdmrg::chi_lim_init"           , xdmrg::chi_lim_init);
    dmrg_config.load_parameter("xdmrg::chi_lim_olap"           , xdmrg::chi_lim_olap);
    dmrg_config.load_parameter("xdmrg::chi_lim_vsub"           , xdmrg::chi_lim_vsub);
    dmrg_config.load_parameter("xdmrg::print_freq "            , xdmrg::print_freq);
    dmrg_config.load_parameter("xdmrg::store_wavefn"           , xdmrg::store_wavefn);
    dmrg_config.load_parameter("xdmrg::energy_density_target"  , xdmrg::energy_density_target);
    dmrg_config.load_parameter("xdmrg::energy_density_window"  , xdmrg::energy_density_window);
    dmrg_config.load_parameter("xdmrg::max_states"             , xdmrg::max_states);
    dmrg_config.load_parameter("xdmrg::finish_if_entanglm_saturated" , xdmrg::finish_if_entanglm_saturated);
    dmrg_config.load_parameter("xdmrg::finish_if_variance_saturated" , xdmrg::finish_if_variance_saturated);

    //Parameters controlling infinite-DMRG
    dmrg_config.load_parameter("idmrg::on"           , idmrg::on);
    dmrg_config.load_parameter("idmrg::max_iters"    , idmrg::max_iters);
    dmrg_config.load_parameter("idmrg::chi_lim_max"  , idmrg::chi_lim_max);
    dmrg_config.load_parameter("idmrg::chi_lim_grow" , idmrg::chi_lim_grow);
    dmrg_config.load_parameter("idmrg::chi_lim_init" , idmrg::chi_lim_init);
    dmrg_config.load_parameter("idmrg::print_freq"   , idmrg::print_freq);


    //Parameters controlling imaginary TEBD (Zero temperature)
    dmrg_config.load_parameter("itebd::on"            , itebd::on       );
    dmrg_config.load_parameter("itebd::max_iters"     , itebd::max_iters);
    dmrg_config.load_parameter("itebd::delta_t0"      , itebd::delta_t0  );
    dmrg_config.load_parameter("itebd::delta_tmin"    , itebd::delta_tmin);
    dmrg_config.load_parameter("itebd::suzuki_order"  , itebd::suzuki_order);
    dmrg_config.load_parameter("itebd::chi_lim_max"   , itebd::chi_lim_max  );
    dmrg_config.load_parameter("itebd::chi_lim_grow"  , itebd::chi_lim_grow);
    dmrg_config.load_parameter("itebd::chi_lim_init"  , itebd::chi_lim_init);
    dmrg_config.load_parameter("itebd::print_freq"    , itebd::print_freq);

    //Profiling
    dmrg_config.load_parameter("profiling::on"        , profiling::on        );
    dmrg_config.load_parameter("profiling::precision" , profiling::precision );
    dmrg_config.load_parameter("profiling::extra"     , profiling::extra     );
    //Console settings
    dmrg_config.load_parameter("console::verbosity"   , console::verbosity);
    dmrg_config.load_parameter("console::timestamp"   , console::timestamp);
    /* clang-format off */
}



void settings::load_config(const std::string & config_filename){
    class_dmrg_config indata(config_filename);
    settings::load_config(indata);
}
