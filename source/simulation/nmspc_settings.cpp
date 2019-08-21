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
    input::input_filename                    = indata.get_input_filename();
    input::input_file                        = indata.get_input_file();
    model::model_type                        = indata.find_parameter<std::string>("model::model_type"                          , model::model_type);
    model::seed_init                         = indata.find_parameter<int>        ("model::seed_init"                           , model::seed_init);
    model::seed_state                        = indata.find_parameter<int>        ("model::seed_state"                          , model::seed_state);
    model::use_seed_state_as_enumeration     = indata.find_parameter<bool>       ("model::use_seed_state_as_enumeration"       , model::use_seed_state_as_enumeration);
    model::project_when_saturated            = indata.find_parameter<bool>       ("model::project_when_saturated"              , model::project_when_saturated);
    model::initial_parity_sector             = indata.find_parameter<std::string>("model::initial_parity_sector"               , model::initial_parity_sector);
    model::tf_ising::J                       = indata.find_parameter<double>     ("model::tf_ising::J"                         , model::tf_ising::J);
    model::tf_ising::g                       = indata.find_parameter<double>     ("model::tf_ising::g"                         , model::tf_ising::g);
    model::tf_ising::w                       = indata.find_parameter<double>     ("model::tf_ising::w"                         , model::tf_ising::w);
    model::tf_ising::d                       = indata.find_parameter<int>        ("model::tf_ising::d"                         , model::tf_ising::d);
    model::tf_nn_ising::J1                   = indata.find_parameter<double>     ("model::tf_nn_ising::J1"                     , model::tf_nn_ising::J1);
    model::tf_nn_ising::J2                   = indata.find_parameter<double>     ("model::tf_nn_ising::J2"                     , model::tf_nn_ising::J2);
    model::tf_nn_ising::g                    = indata.find_parameter<double>     ("model::tf_nn_ising::g"                      , model::tf_nn_ising::g);
    model::tf_nn_ising::d                    = indata.find_parameter<int>        ("model::tf_nn_ising::d"                      , model::tf_nn_ising::d);
    model::tf_nn_ising::w                    = indata.find_parameter<double>     ("model::tf_nn_ising::w"                      , model::tf_nn_ising::w);
    model::selfdual_tf_rf_ising::J_log_mean  = indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::J_log_mean"    , model::selfdual_tf_rf_ising::J_log_mean);
    model::selfdual_tf_rf_ising::h_log_mean  = indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::h_log_mean"    , model::selfdual_tf_rf_ising::h_log_mean);
    model::selfdual_tf_rf_ising::J_sigma     = indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::J_sigma"       , model::selfdual_tf_rf_ising::J_sigma);
    model::selfdual_tf_rf_ising::h_sigma     = indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::h_sigma"       , model::selfdual_tf_rf_ising::h_sigma);
    model::selfdual_tf_rf_ising::lambda      = indata.find_parameter<double>     ("model::selfdual_tf_rf_ising::lambda"        , model::selfdual_tf_rf_ising::lambda);
    model::selfdual_tf_rf_ising::d           = indata.find_parameter<int>        ("model::selfdual_tf_rf_ising::d"             , model::selfdual_tf_rf_ising::d);
    precision::eigMaxIter                    = indata.find_parameter<int>        ("precision::eigMaxIter"                      , precision::eigMaxIter);
    precision::eigThreshold                  = indata.find_parameter<double>     ("precision::eigThreshold"                    , precision::eigThreshold);
    precision::eigMaxNcv                     = indata.find_parameter<int>        ("precision::eigMaxNcv"                       , precision::eigMaxNcv);
    precision::SVDThreshold                  = indata.find_parameter<double>     ("precision::SVDThreshold"                    , precision::SVDThreshold);
    precision::VarConvergenceThreshold       = indata.find_parameter<double>     ("precision::VarConvergenceThreshold"         , precision::VarConvergenceThreshold);
    precision::VarSaturationThreshold        = indata.find_parameter<double>     ("precision::VarSaturationThreshold"          , precision::VarSaturationThreshold);
    precision::EntEntrSaturationThreshold    = indata.find_parameter<double>     ("precision::EntEntrSaturationThreshold"      , precision::EntEntrSaturationThreshold);
    precision::SubspaceQualityFactor         = indata.find_parameter<double>     ("precision::SubspaceQualityFactor"           , precision::SubspaceQualityFactor);
    precision::MaxSubspaceQuality            = indata.find_parameter<double>     ("precision::MaxSubspaceQuality"              , precision::MaxSubspaceQuality);
    precision::MaxSitesMultiDmrg             = indata.find_parameter<int>        ("precision::MaxSitesMultiDmrg"               , precision::MaxSitesMultiDmrg);
    precision::MaxSizeFullDiag               = indata.find_parameter<int>        ("precision::MaxSizeFullDiag"                 , precision::MaxSizeFullDiag);
    precision::MaxSizePartDiag               = indata.find_parameter<int>        ("precision::MaxSizePartDiag"                 , precision::MaxSizePartDiag);
    precision::MaxSizeDirect                 = indata.find_parameter<int>        ("precision::MaxSizeDirect"                   , precision::MaxSizeDirect);

    //Parameters controlling infinite-DMRG
    idmrg::on                           = indata.find_parameter<bool>   ("idmrg::on"         , idmrg::on);
    if(idmrg::on){
        idmrg::max_steps                = indata.find_parameter<int>    ("idmrg::max_steps"  , idmrg::max_steps);
        idmrg::chi_max                  = indata.find_parameter<long>   ("idmrg::chi_max"    , idmrg::chi_max);
        idmrg::chi_grow                 = indata.find_parameter<bool>   ("idmrg::chi_grow"   , idmrg::chi_grow);
        idmrg::print_freq               = indata.find_parameter<int>    ("idmrg::print_freq" , idmrg::print_freq);
        idmrg::write_freq               = indata.find_parameter<int>    ("idmrg::write_freq" , idmrg::write_freq);
    }

    //Parameters controlling finite-DMRG
    fdmrg::on                     = indata.find_parameter<bool>   ("fdmrg::on"           , fdmrg::on);
    if(fdmrg::on){
        fdmrg::num_sites          = indata.find_parameter<int>    ("fdmrg::num_sites "   , fdmrg::num_sites);
        fdmrg::max_sweeps         = indata.find_parameter<int>    ("fdmrg::max_sweeps "  , fdmrg::max_sweeps);
        fdmrg::min_sweeps         = indata.find_parameter<int>    ("fdmrg::min_sweeps "  , fdmrg::min_sweeps);
        fdmrg::chi_max            = indata.find_parameter<int>    ("fdmrg::chi_max"      , 8);
        fdmrg::chi_grow           = indata.find_parameter<bool>   ("fdmrg::chi_grow"     , fdmrg::chi_grow);
        fdmrg::print_freq         = indata.find_parameter<int>    ("fdmrg::print_freq "  , fdmrg::print_freq);
        fdmrg::write_freq         = indata.find_parameter<int>    ("fdmrg::write_freq "  , fdmrg::write_freq);
        fdmrg::store_wavefn       = indata.find_parameter<bool>   ("fdmrg::store_wavefn" , fdmrg::store_wavefn);

    }

    //Parameters controlling excited state DMRG
    xdmrg::on                           = indata.find_parameter<bool>   ("xdmrg::on"         , xdmrg::on);
    if(xdmrg::on){
        xdmrg::num_sites                = indata.find_parameter<int>    ("xdmrg::num_sites "             , xdmrg::num_sites);
        xdmrg::max_sweeps               = indata.find_parameter<int>    ("xdmrg::max_sweeps "            , xdmrg::max_sweeps);
        xdmrg::min_sweeps               = indata.find_parameter<int>    ("xdmrg::min_sweeps "            , xdmrg::min_sweeps);
        xdmrg::chi_max                  = indata.find_parameter<int>    ("xdmrg::chi_max"                , xdmrg::chi_max);
        xdmrg::chi_grow                 = indata.find_parameter<bool>   ("xdmrg::chi_grow"               , xdmrg::chi_grow);
        xdmrg::print_freq               = indata.find_parameter<int>    ("xdmrg::print_freq "            , xdmrg::print_freq);
        xdmrg::write_freq               = indata.find_parameter<int>    ("xdmrg::write_freq "            , xdmrg::write_freq);
        xdmrg::store_wavefn             = indata.find_parameter<bool>   ("xdmrg::store_wavefn"           , xdmrg::store_wavefn);
        xdmrg::energy_density_target    = indata.find_parameter<double> ("xdmrg::energy_density_target"  , xdmrg::energy_density_target);
        xdmrg::energy_density_window    = indata.find_parameter<double> ("xdmrg::energy_density_window"  , xdmrg::energy_density_window);
    }


    //Parameters controlling imaginary TEBD (Zero temperature)
    itebd::on                     = indata.find_parameter<bool>   ("itebd::on"          , itebd::on       );
    if(itebd::on){
        itebd::max_steps          = indata.find_parameter<int>    ("itebd::max_steps "  , itebd::max_steps);
        itebd::delta_t0           = indata.find_parameter<double> ("itebd::delta_t0"    , itebd::delta_t0  );
        itebd::delta_tmin         = indata.find_parameter<double> ("itebd::delta_tmin"  , itebd::delta_tmin);
        itebd::suzuki_order       = indata.find_parameter<int>    ("itebd::suzuki_order", itebd::suzuki_order);
        itebd::chi_max            = indata.find_parameter<long>   ("itebd::chi_max"     , itebd::chi_max  );
        itebd::chi_grow           = indata.find_parameter<bool>   ("itebd::chi_grow"    , itebd::chi_grow);
        itebd::print_freq         = indata.find_parameter<int>    ("itebd::print_freq"  , itebd::print_freq);
        itebd::write_freq         = indata.find_parameter<int>    ("itebd::write_freq"  , itebd::write_freq);
    }

    //Save data_struct to output
    output::save_logs                = indata.find_parameter<bool>   ("output::save_logs"               , output::save_logs          );
    output::save_profiling           = indata.find_parameter<bool>   ("output::save_profiling"          , output::save_profiling        );
    output::output_filename          = indata.find_parameter<string> ("output::output_filename"         , output::output_filename);
    output::access_mode              = indata.find_parameter<string> ("output::access_mode"             , output::access_mode);
    output::create_mode              = indata.find_parameter<string> ("output::create_mode"             , output::create_mode);
    int storageLevelRead           = indata.find_parameter<int>    ("output::storage_level"           , 2       );
    output::storage_level            = static_cast<StorageLevel>     (storageLevelRead );

    //Profiling
    profiling::on                  = indata.find_parameter<bool>   ("profiling::on"        , profiling::on        );
    profiling::precision           = indata.find_parameter<int>    ("profiling::precision" , profiling::precision );
    //Console settings
    console::verbosity             = indata.find_parameter<int>    ("console::verbosity"   , console::verbosity);
    console::timestamp             = indata.find_parameter<bool>   ("console::timestamp"   , console::timestamp);
}

void settings::load_from_hdf5(h5pp::File & h5ppFile){

    std::string settings_from_hdf5;
    std::string temp_filename = "indata_temp.cfg";
    h5ppFile.readDataset(settings_from_hdf5, "/common/input_file");

    std::ofstream temp_settings_file(temp_filename);
    temp_settings_file << settings_from_hdf5;
    temp_settings_file.close();
    class_settings_reader indata(temp_filename);
    settings::load_from_file(indata);

}
