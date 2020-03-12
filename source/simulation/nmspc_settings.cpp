//
// Created by david on 2018-01-14.
//

#include "nmspc_settings.h"
#include <io/class_settings_reader.h>
#include <h5pp/h5pp.h>
using namespace std;

/*
 * The section below attempts to find and read the parameters from the given inputfile.
 *
*/


void settings::load_from_file(class_settings_reader &indata){
    input::input_filename             = indata.get_input_filename();
    input::input_file_raw             = indata.get_input_file_as_string();
    indata.find_parameter<std::string>("model::model_type"                          , model::model_type);
    indata.find_parameter<long>       ("model::seed"                                , model::seed);
    indata.find_parameter<long>       ("model::state_number"                        , model::state_number);
    indata.find_parameter<double>     ("model::tf_ising::J"                         , model::tf_ising::J);
    indata.find_parameter<double>     ("model::tf_ising::g"                         , model::tf_ising::g);
    indata.find_parameter<double>     ("model::tf_ising::w"                         , model::tf_ising::w);
    indata.find_parameter<size_t>     ("model::tf_ising::d"                         , model::tf_ising::d);
    indata.find_parameter<double>     ("model::tf_nn_ising::J1"                     , model::tf_nn_ising::J1);
    indata.find_parameter<double>     ("model::tf_nn_ising::J2"                     , model::tf_nn_ising::J2);
    indata.find_parameter<double>     ("model::tf_nn_ising::g"                      , model::tf_nn_ising::g);
    indata.find_parameter<size_t>     ("model::tf_nn_ising::d"                      , model::tf_nn_ising::d);
    indata.find_parameter<double>     ("model::tf_nn_ising::w"                      , model::tf_nn_ising::w);
    indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::J_mean"        , model::selfdual_tf_rf_ising::J_mean);
    indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::h_mean"        , model::selfdual_tf_rf_ising::h_mean);
    indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::J_sigma"       , model::selfdual_tf_rf_ising::J_sigma);
    indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::h_sigma"       , model::selfdual_tf_rf_ising::h_sigma);
    indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::lambda"        , model::selfdual_tf_rf_ising::lambda);
    indata.find_parameter<bool>       ("model::selfdual_tf_rf_ising::parity_sep"    , model::selfdual_tf_rf_ising::parity_sep);
    indata.find_parameter<std::string>("model::selfdual_tf_rf_ising::distribution"  , model::selfdual_tf_rf_ising::distribution);
    indata.find_parameter<size_t>     ("model::selfdual_tf_rf_ising::d"             , model::selfdual_tf_rf_ising::d);
    indata.find_parameter<bool>       ("strategy::chi_quench_when_stuck"            , strategy::chi_quench_when_stuck);
    indata.find_parameter<bool>       ("strategy::perturb_when_stuck"               , strategy::perturb_when_stuck);
    indata.find_parameter<bool>       ("strategy::damping_when_stuck"               , strategy::damping_when_stuck);
    indata.find_parameter<bool>       ("strategy::project_trial_when_stuck"         , strategy::project_when_stuck);
    indata.find_parameter<bool>       ("strategy::project_on_every_sweep"           , strategy::project_on_every_sweep);
    indata.find_parameter<bool>       ("strategy::project_on_chi_update"            , strategy::project_on_chi_update);
    indata.find_parameter<bool>       ("strategy::randomize_on_chi_update"          , strategy::randomize_on_chi_update);
    indata.find_parameter<bool>       ("strategy::randomize_early"                  , strategy::randomize_early);
    indata.find_parameter<bool>       ("strategy::use_pauli_eigvecs"                , strategy::use_pauli_eigvecs);
    indata.find_parameter<std::string>("strategy::initial_parity_sector"            , strategy::initial_parity_sector);
    indata.find_parameter<std::string>("strategy::target_parity_sector"             , strategy::target_parity_sector);
    indata.find_parameter<size_t>     ("precision::eig_max_iter"                    , precision::eig_max_iter);
    indata.find_parameter<double>     ("precision::eig_threshold"                   , precision::eig_threshold);
    indata.find_parameter<size_t>     ("precision::eig_max_ncv"                     , precision::eig_max_ncv);
    indata.find_parameter<double>     ("precision::svd_threshold"                   , precision::svd_threshold);
    indata.find_parameter<double>     ("precision::variance_convergence_threshold"  , precision::variance_convergence_threshold);
    indata.find_parameter<double>     ("precision::variance_slope_threshold"        , precision::variance_slope_threshold);
    indata.find_parameter<double>     ("precision::entropy_slope_threshold"         , precision::entropy_slope_threshold);
    indata.find_parameter<double>     ("precision::subspace_error_factor"           , precision::subspace_error_factor);
    indata.find_parameter<double>     ("precision::max_subspace_error"              , precision::max_subspace_error);
    indata.find_parameter<double>     ("precision::min_subspace_error"              , precision::min_subspace_error);
    indata.find_parameter<size_t>     ("precision::max_size_full_diag"              , precision::max_size_full_diag);
    indata.find_parameter<size_t>     ("precision::max_size_part_diag"              , precision::max_size_part_diag);
    indata.find_parameter<size_t>     ("precision::max_size_direct"                 , precision::max_size_direct);
    indata.find_parameter<double>     ("precision::max_norm_error"                  , precision::max_norm_error);
    indata.find_parameter<size_t>     ("precision::max_resets"                      , precision::max_resets);
    indata.find_parameter<bool>       ("precision::use_reduced_energy"              , precision::use_reduced_energy);
    indata.find_parameter<size_t>     ("precision::max_sites_multidmrg"             , precision::max_sites_multidmrg);
    indata.find_parameter<std::string>("precision::move_sites_multidmrg"            , precision::move_sites_multidmrg);

    indata.find_parameter<int>        ("threading::num_threads"               , threading::num_threads);


    //Parameters controlling infinite-DMRG
    indata.find_parameter<bool>   ("idmrg::on"         , idmrg::on);
    if(idmrg::on){
        indata.find_parameter<size_t> ("idmrg::max_steps"  , idmrg::max_steps);
        indata.find_parameter<long>   ("idmrg::chi_max"    , idmrg::chi_max);
        indata.find_parameter<bool>   ("idmrg::chi_grow"   , idmrg::chi_grow);
        indata.find_parameter<long>   ("idmrg::chi_init"   , idmrg::chi_init);
        indata.find_parameter<size_t> ("idmrg::print_freq" , idmrg::print_freq);
        indata.find_parameter<size_t> ("idmrg::write_freq" , idmrg::write_freq);
    }

    //Parameters controlling finite-DMRG
    indata.find_parameter<bool>   ("fdmrg::on"           , fdmrg::on);
    if(fdmrg::on){
        indata.find_parameter<size_t> ("fdmrg::num_sites "   , fdmrg::num_sites);
        indata.find_parameter<size_t> ("fdmrg::max_sweeps "  , fdmrg::max_sweeps);
        indata.find_parameter<size_t> ("fdmrg::min_sweeps "  , fdmrg::min_sweeps);
        indata.find_parameter<long>   ("fdmrg::chi_max"      , fdmrg::chi_max);
        indata.find_parameter<bool>   ("fdmrg::chi_grow"     , fdmrg::chi_grow);
        indata.find_parameter<long>   ("fdmrg::chi_init"     , fdmrg::chi_init);
        indata.find_parameter<size_t> ("fdmrg::print_freq "  , fdmrg::print_freq);
        indata.find_parameter<size_t> ("fdmrg::write_freq "  , fdmrg::write_freq);
        indata.find_parameter<bool>   ("fdmrg::store_wavefn" , fdmrg::store_wavefn);

    }

    //Parameters controlling excited state DMRG
    indata.find_parameter<bool>   ("xdmrg::on"         , xdmrg::on);
    if(xdmrg::on){
        indata.find_parameter<size_t> ("xdmrg::num_sites "             , xdmrg::num_sites);
        indata.find_parameter<size_t> ("xdmrg::max_sweeps "            , xdmrg::max_sweeps);
        indata.find_parameter<size_t> ("xdmrg::min_sweeps "            , xdmrg::min_sweeps);
        indata.find_parameter<long>   ("xdmrg::chi_max"                , xdmrg::chi_max);
        indata.find_parameter<bool>   ("xdmrg::chi_grow"               , xdmrg::chi_grow);
        indata.find_parameter<long>   ("xdmrg::chi_init"               , xdmrg::chi_init);
        indata.find_parameter<size_t> ("xdmrg::print_freq "            , xdmrg::print_freq);
        indata.find_parameter<size_t> ("xdmrg::write_freq "            , xdmrg::write_freq);
        indata.find_parameter<bool>   ("xdmrg::store_wavefn"           , xdmrg::store_wavefn);
        indata.find_parameter<double> ("xdmrg::energy_density_target"  , xdmrg::energy_density_target);
        indata.find_parameter<double> ("xdmrg::energy_density_window"  , xdmrg::energy_density_window);
        indata.find_parameter<size_t> ("xdmrg::max_states"             , xdmrg::max_states);
    }


    //Parameters controlling imaginary TEBD (Zero temperature)
    indata.find_parameter<bool>   ("itebd::on"          , itebd::on       );
    if(itebd::on){
        indata.find_parameter<size_t> ("itebd::max_steps "  , itebd::max_steps);
        indata.find_parameter<double> ("itebd::delta_t0"    , itebd::delta_t0  );
        indata.find_parameter<double> ("itebd::delta_tmin"  , itebd::delta_tmin);
        indata.find_parameter<size_t> ("itebd::suzuki_order", itebd::suzuki_order);
        indata.find_parameter<long>   ("itebd::chi_max"     , itebd::chi_max  );
        indata.find_parameter<bool>   ("itebd::chi_grow"    , itebd::chi_grow);
        indata.find_parameter<long>   ("fdmrg::chi_init"    , itebd::chi_init);
        indata.find_parameter<size_t> ("itebd::print_freq"  , itebd::print_freq);
        indata.find_parameter<size_t> ("itebd::write_freq"  , itebd::write_freq);
    }

    //Save data_struct to output
    indata.find_parameter<bool>   ("output::save_logs"               , output::save_logs );
    indata.find_parameter<bool>   ("output::save_profiling"          , output::save_profiling);
    indata.find_parameter<string> ("output::output_filename"         , output::output_filename);
    indata.find_parameter<string> ("output::access_mode"             , output::access_mode);
    indata.find_parameter<string> ("output::create_mode"             , output::create_mode);
    indata.find_parameter<bool>   ("output::use_temp_dir"            , output::use_temp_dir);
    indata.find_parameter<size_t> ("output::copy_from_temp_freq"     , output::copy_from_temp_freq);
    indata.find_parameter<string> ("output::temp_dir"                , output::temp_dir);

//    inline bool         use_temp_dir         = true;                         /*!< If true uses a temporary directory for writes in the local drive (usually /tmp) and copies the results afterwards */
//    inline size_t       copy_from_temp_freq  = 4;                            /*!< How often, in units of iterations, to copy the hdf5 file in tmp dir to target destination */
//    inline std::string  temp_dir             = "/scratch/local";             /*!< Local temp directory on the "local" system. If it doesn't exist we default to /tmp instead (or whatever is the default */
    int storageLevelRead = 2;
    indata.find_parameter<int>    ("output::storage_level"           , storageLevelRead );
    output::storage_level            = static_cast<StorageLevel>     (storageLevelRead);

    //Profiling
    indata.find_parameter<bool>   ("profiling::on"        , profiling::on        );
    indata.find_parameter<size_t> ("profiling::precision" , profiling::precision );
    //Console settings
    indata.find_parameter<size_t> ("console::verbosity"   , console::verbosity);
    indata.find_parameter<bool>   ("console::timestamp"   , console::timestamp);
}

void settings::load_from_hdf5(h5pp::File & h5ppFile){

    std::string settings_from_hdf5;
    std::string temp_filename = "indata_temp.cfg";
    h5ppFile.readDataset(settings_from_hdf5, "/common/input_filepath");

    std::ofstream temp_settings_file(temp_filename);
    temp_settings_file << settings_from_hdf5;
    temp_settings_file.close();
    class_settings_reader indata(temp_filename);
    settings::load_from_file(indata);

}
