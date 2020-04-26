//
// Created by david on 2018-01-14.
//

#include "nmspc_settings.h"
#include <h5pp/h5pp.h>
#include <io/class_config_reader.h>
#include <simulation/class_dmrg_config.h>
using namespace std;

/*
 * The section below attempts to find and read the parameters from the given inputfile.
 *
 */

void settings::load_config(class_dmrg_config &dmrg_config) {
    if(not dmrg_config.file_exists) throw std::runtime_error(fmt::format("Could not load config [{}]: File does not exist", dmrg_config.file_path));
    dmrg_config.load();
    input::config_filename      = dmrg_config.file_path.string();
    input::config_file_contents = dmrg_config.get_config_file_as_string();
    /* clang-format off */
    dmrg_config.load_parameter<FileCollisionPolicy>     ("output::existing_file_policy"                 , output::file_collision_policy);
    dmrg_config.load_parameter<long>                    ("input::seed"                                  , input::seed);
    dmrg_config.load_parameter<long>                    ("input::state_number"                          , input::bitfield);
    dmrg_config.load_parameter<string>                  ("output::output_filepath"                      , output::output_filepath);
    dmrg_config.load_parameter<bool>                    ("output::save_profiling"                       , output::save_profiling);
    dmrg_config.load_parameter<bool>                    ("output::journal_keep_only_last_iter"          , output::journal_keep_only_last_iter);
    dmrg_config.load_parameter<bool>                    ("output::use_temp_dir"                         , output::use_temp_dir);
    dmrg_config.load_parameter<size_t>                  ("output::copy_from_temp_freq"                  , output::copy_from_temp_freq);
    dmrg_config.load_parameter<string>                  ("output::temp_dir"                             , output::temp_dir);
    dmrg_config.load_parameter<size_t>                  ("output::compression_level"                    , output::compression_level);
    dmrg_config.load_parameter<StorageLevel>            ("output::storage_level_journal"                , output::storage_level_journal);
    dmrg_config.load_parameter<StorageLevel>            ("output::storage_level_results"                , output::storage_level_results);
    dmrg_config.load_parameter<StorageLevel>            ("output::storage_level_chi_update"             , output::storage_level_chi_update);
    dmrg_config.load_parameter<StorageLevel>            ("output::storage_level_projection"             , output::storage_level_projection);
    dmrg_config.load_parameter<std::string>             ("model::model_type"                            , model::model_type);
    dmrg_config.load_parameter<size_t>                  ("model::sites"                                 , model::sites);
    dmrg_config.load_parameter<double>                  ("model::tf_ising::J"                           , model::tf_ising::J);
    dmrg_config.load_parameter<double>                  ("model::tf_ising::g"                           , model::tf_ising::g);
    dmrg_config.load_parameter<double>                  ("model::tf_ising::w"                           , model::tf_ising::w);
    dmrg_config.load_parameter<size_t>                  ("model::tf_ising::d"                           , model::tf_ising::d);
    dmrg_config.load_parameter<double>                  ("model::tf_nn_ising::J1"                       , model::tf_nn_ising::J1);
    dmrg_config.load_parameter<double>                  ("model::tf_nn_ising::J2"                       , model::tf_nn_ising::J2);
    dmrg_config.load_parameter<double>                  ("model::tf_nn_ising::g"                        , model::tf_nn_ising::g);
    dmrg_config.load_parameter<size_t>                  ("model::tf_nn_ising::d"                        , model::tf_nn_ising::d);
    dmrg_config.load_parameter<double>                  ("model::tf_nn_ising::w"                        , model::tf_nn_ising::w);
    dmrg_config.load_parameter<double>                  ("model::selfdual_tf_rf_ising::J_mean"          , model::selfdual_tf_rf_ising::J_mean);
    dmrg_config.load_parameter<double>                  ("model::selfdual_tf_rf_ising::h_mean"          , model::selfdual_tf_rf_ising::h_mean);
    dmrg_config.load_parameter<double>                  ("model::selfdual_tf_rf_ising::J_sigma"         , model::selfdual_tf_rf_ising::J_sigma);
    dmrg_config.load_parameter<double>                  ("model::selfdual_tf_rf_ising::h_sigma"         , model::selfdual_tf_rf_ising::h_sigma);
    dmrg_config.load_parameter<double>                  ("model::selfdual_tf_rf_ising::lambda"          , model::selfdual_tf_rf_ising::lambda);
    dmrg_config.load_parameter<bool>                    ("model::selfdual_tf_rf_ising::parity_sep"      , model::selfdual_tf_rf_ising::parity_sep);
    dmrg_config.load_parameter<std::string>             ("model::selfdual_tf_rf_ising::distribution"    , model::selfdual_tf_rf_ising::distribution);
    dmrg_config.load_parameter<size_t>                  ("model::selfdual_tf_rf_ising::d"               , model::selfdual_tf_rf_ising::d);
    dmrg_config.load_parameter<bool>                    ("strategy::chi_quench_when_stuck"              , strategy::chi_quench_when_stuck);
    dmrg_config.load_parameter<bool>                    ("strategy::perturb_when_stuck"                 , strategy::perturb_when_stuck);
    dmrg_config.load_parameter<bool>                    ("strategy::damping_when_stuck"                 , strategy::damping_when_stuck);
    dmrg_config.load_parameter<bool>                    ("strategy::project_trial_when_stuck"           , strategy::project_when_stuck);
    dmrg_config.load_parameter<bool>                    ("strategy::project_on_every_sweep"             , strategy::project_on_every_sweep);
    dmrg_config.load_parameter<bool>                    ("strategy::project_on_chi_update"              , strategy::project_on_chi_update);
    dmrg_config.load_parameter<bool>                    ("strategy::randomize_on_chi_update"            , strategy::randomize_on_chi_update);
    dmrg_config.load_parameter<bool>                    ("strategy::randomize_early"                    , strategy::randomize_early);
    dmrg_config.load_parameter<bool>                    ("strategy::use_pauli_eigvecs"                  , strategy::use_pauli_eigvecs);
    dmrg_config.load_parameter<std::string>             ("strategy::initial_parity_sector"              , strategy::initial_parity_sector);
    dmrg_config.load_parameter<std::string>             ("strategy::target_parity_sector"               , strategy::target_parity_sector);
    dmrg_config.load_parameter<size_t>                  ("precision::eig_max_iter"                      , precision::eig_max_iter);
    dmrg_config.load_parameter<double>                  ("precision::eig_threshold"                     , precision::eig_threshold);
    dmrg_config.load_parameter<size_t>                  ("precision::eig_max_ncv"                       , precision::eig_max_ncv);
    dmrg_config.load_parameter<double>                  ("precision::svd_threshold"                     , precision::svd_threshold);
    dmrg_config.load_parameter<double>                  ("precision::variance_convergence_threshold"    , precision::variance_convergence_threshold);
    dmrg_config.load_parameter<double>                  ("precision::variance_slope_threshold"          , precision::variance_slope_threshold);
    dmrg_config.load_parameter<double>                  ("precision::entropy_slope_threshold"           , precision::entropy_slope_threshold);
    dmrg_config.load_parameter<double>                  ("precision::subspace_error_factor"             , precision::subspace_error_factor);
    dmrg_config.load_parameter<double>                  ("precision::max_subspace_error"                , precision::max_subspace_error);
    dmrg_config.load_parameter<double>                  ("precision::min_subspace_error"                , precision::min_subspace_error);
    dmrg_config.load_parameter<size_t>                  ("precision::max_size_full_diag"                , precision::max_size_full_diag);
    dmrg_config.load_parameter<size_t>                  ("precision::max_size_part_diag"                , precision::max_size_part_diag);
    dmrg_config.load_parameter<size_t>                  ("precision::max_size_direct"                   , precision::max_size_direct);
    dmrg_config.load_parameter<double>                  ("precision::max_norm_error"                    , precision::max_norm_error);
    dmrg_config.load_parameter<size_t>                  ("precision::max_resets"                        , precision::max_resets);
    dmrg_config.load_parameter<bool>                    ("precision::use_reduced_energy"                , precision::use_reduced_energy);
    dmrg_config.load_parameter<size_t>                  ("precision::max_sites_multidmrg"               , precision::max_sites_multidmrg);
    dmrg_config.load_parameter<std::string>             ("precision::move_sites_multidmrg"              , precision::move_sites_multidmrg);
    dmrg_config.load_parameter<int>                     ("threading::num_threads"                       , threading::num_threads);


    //Parameters controlling infinite-DMRG
    dmrg_config.load_parameter<bool>   ("idmrg::on"         , idmrg::on);
    if(idmrg::on){
        dmrg_config.load_parameter<size_t> ("idmrg::max_steps"  , idmrg::max_steps);
        dmrg_config.load_parameter<long>   ("idmrg::chi_max"    , idmrg::chi_max);
        dmrg_config.load_parameter<bool>   ("idmrg::chi_grow"   , idmrg::chi_grow);
        dmrg_config.load_parameter<long>   ("idmrg::chi_init"   , idmrg::chi_init);
        dmrg_config.load_parameter<size_t> ("idmrg::print_freq" , idmrg::print_freq);
        dmrg_config.load_parameter<size_t> ("idmrg::write_freq" , idmrg::write_freq);
    }

    //Parameters controlling finite-DMRG
    dmrg_config.load_parameter<bool>   ("fdmrg::on"           , fdmrg::on);
    if(fdmrg::on){
        dmrg_config.load_parameter<size_t> ("fdmrg::max_sweeps "  , fdmrg::max_sweeps);
        dmrg_config.load_parameter<size_t> ("fdmrg::min_sweeps "  , fdmrg::min_sweeps);
        dmrg_config.load_parameter<long>   ("fdmrg::chi_max"      , fdmrg::chi_max);
        dmrg_config.load_parameter<bool>   ("fdmrg::chi_grow"     , fdmrg::chi_grow);
        dmrg_config.load_parameter<long>   ("fdmrg::chi_init"     , fdmrg::chi_init);
        dmrg_config.load_parameter<size_t> ("fdmrg::print_freq "  , fdmrg::print_freq);
        dmrg_config.load_parameter<size_t> ("fdmrg::write_freq "  , fdmrg::write_freq);
        dmrg_config.load_parameter<bool>   ("fdmrg::store_wavefn" , fdmrg::store_wavefn);

    }

    //Parameters controlling excited state DMRG
    dmrg_config.load_parameter<bool>   ("xdmrg::on"         , xdmrg::on);
    if(xdmrg::on){
        dmrg_config.load_parameter<size_t> ("xdmrg::max_sweeps "            , xdmrg::max_sweeps);
        dmrg_config.load_parameter<size_t> ("xdmrg::min_sweeps "            , xdmrg::min_sweeps);
        dmrg_config.load_parameter<long>   ("xdmrg::chi_max"                , xdmrg::chi_max);
        dmrg_config.load_parameter<bool>   ("xdmrg::chi_grow"               , xdmrg::chi_grow);
        dmrg_config.load_parameter<long>   ("xdmrg::chi_init"               , xdmrg::chi_init);
        dmrg_config.load_parameter<size_t> ("xdmrg::print_freq "            , xdmrg::print_freq);
        dmrg_config.load_parameter<size_t> ("xdmrg::write_freq "            , xdmrg::write_freq);
        dmrg_config.load_parameter<bool>   ("xdmrg::store_wavefn"           , xdmrg::store_wavefn);
        dmrg_config.load_parameter<double> ("xdmrg::energy_density_target"  , xdmrg::energy_density_target);
        dmrg_config.load_parameter<double> ("xdmrg::energy_density_window"  , xdmrg::energy_density_window);
        dmrg_config.load_parameter<size_t> ("xdmrg::max_states"             , xdmrg::max_states);
    }


    //Parameters controlling imaginary TEBD (Zero temperature)
    dmrg_config.load_parameter<bool>   ("itebd::on"          , itebd::on       );
    if(itebd::on){
        dmrg_config.load_parameter<size_t> ("itebd::max_steps "  , itebd::max_steps);
        dmrg_config.load_parameter<double> ("itebd::delta_t0"    , itebd::delta_t0  );
        dmrg_config.load_parameter<double> ("itebd::delta_tmin"  , itebd::delta_tmin);
        dmrg_config.load_parameter<size_t> ("itebd::suzuki_order", itebd::suzuki_order);
        dmrg_config.load_parameter<long>   ("itebd::chi_max"     , itebd::chi_max  );
        dmrg_config.load_parameter<bool>   ("itebd::chi_grow"    , itebd::chi_grow);
        dmrg_config.load_parameter<long>   ("fdmrg::chi_init"    , itebd::chi_init);
        dmrg_config.load_parameter<size_t> ("itebd::print_freq"  , itebd::print_freq);
        dmrg_config.load_parameter<size_t> ("itebd::write_freq"  , itebd::write_freq);
    }


    //Profiling
    dmrg_config.load_parameter<bool>   ("profiling::on"        , profiling::on        );
    dmrg_config.load_parameter<size_t> ("profiling::precision" , profiling::precision );
    //Console settings
    dmrg_config.load_parameter<size_t> ("console::verbosity"   , console::verbosity);
    dmrg_config.load_parameter<bool>   ("console::timestamp"   , console::timestamp);
    /* clang-format off */
}



void settings::load_config(const std::string & config_filename){
    class_dmrg_config indata(config_filename);
    settings::load_config(indata);
}
