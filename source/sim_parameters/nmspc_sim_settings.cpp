//
// Created by david on 2018-01-14.
//

#include "nmspc_sim_settings.h"
#include <io/class_settings_reader.h>
#include <h5pp/h5pp.h>
using namespace std;


/*
 * The section below defines default values for all the parameters,
 * in case an input file is not found, or if values are missing in that file.
 *
*/

namespace settings{
    std::string input::input_file;
    std::string input::input_filename;
    //Common parameters for the model Hamiltonian
    std::string model::model_type        = "tf_ising";                    /*!< The default choice of model type from the enum */
    std::string model::initial_state     = "rps";                         /*!< Choose initial state of the MPS: {upup, updown, GHZ(upup+downdown), W(updown+downup), rps (random product state), random_chi (random state with bond dimension chi)} "cat" or "random". Default "rps". Note that "random_chi" works poorly for finite algorithms */
    int         model::seed_init_mpo     = 1;                             /*!< Seed for the random number generator if you use random fields in the Hamiltonian. */
    int         model::seed_init_mps     = 1;                             /*!< Seed for the random number generator when selecting the initial random product state. */


    //Parameters for the transverse-field Ising model
    double model::tf_ising::J            = 1;                             /*!< Ferromagnetic coupling. J < 0  Gives a ferromagnet. J > 0 an antiferromagnet. */
    double model::tf_ising::g            = 1;                             /*!< Transverse field strength */
    double model::tf_ising::w            = 0;                             /*!< Randomness strength for the random field */
    int    model::tf_ising::d            = 2;                             /*!< Local spin dimension */

    //Parameters for the transvese-field next-nearest neighbor Ising model
    double model::tf_nn_ising::J1        = 1;                             /*!< Ferromagnetic coupling for nearest neighbors.*/
    double model::tf_nn_ising::J2        = 1;                             /*!< Ferromagnetic coupling for next-nearest neighbors.*/
    double model::tf_nn_ising::g         = 1;                             /*!< Transverse field strength */
    double model::tf_nn_ising::w         = 0;                             /*!< Randomness strength for the random field */
    int    model::tf_nn_ising::d         = 2;                             /*!< Local spin dimension */

    //Parameters for the selfdual transvese-field random-field next-neighbor Ising model
    double model::selfdual_tf_rf_ising::J_log_mean = 0;                  /*!< Average ferromagnetic coupling strength.*/
    double model::selfdual_tf_rf_ising::h_log_mean = 0;                  /*!< Average transverse magnetic field strength */
    double model::selfdual_tf_rf_ising::J_sigma    = 1;                  /*!< Standard deviation for the lognormal distribution, i.e. = std(log(J)) , for the ferromagnetic coupling */
    double model::selfdual_tf_rf_ising::h_sigma    = 0;                  /*!< Standard deviation for the lognormal distribution, i.e. = std(log(h))   for the transverse magnetic field */
    double model::selfdual_tf_rf_ising::lambda     = 0;                  /*!< Lambda parameter */
    int    model::selfdual_tf_rf_ising::d          = 2;                  /*!< Local spin dimension */



    int    precision::eigMaxIter                   = 1000  ;
    double precision::eigThreshold                 = 1e-12 ;
    int    precision::eigMaxNcv                    = 16    ;
    double precision::SVDThreshold                 = 1e-12 ;
    double precision::VarConvergenceThreshold      = 1e-8  ;            /*!< Variance convergence threshold. The MPS state is considered good enough when its variance reaches below this value */
    double precision::VarSaturationThreshold       = 1e-4  ;            /*!< Variance saturation  threshold. The variance has saturated when its (absolute) slope reaches below this value */
    double precision::EntEntrSaturationThreshold   = 1e-4  ;            /*!< Entanglement Entropy saturation threshold. The entanglement entropy has saturated when its (absolute) slope reaches below this value*/
    int    precision::MaxSizeFullDiag              = 2048  ;            /*!< Maximum linear size allowed for full diagonalization of the local hamiltonian matrix. */


    //Parameters controlling infinite-DMRG
    bool   idmrg::on                     = true;
    int    idmrg::max_steps              = 5000;
    long   idmrg::chi_max                = 8;
    bool   idmrg::chi_grow               = true;
    int    idmrg::print_freq             = 1000;
    int    idmrg::store_freq             = 100;

    //Parameters controlling Finite-DMRG
    bool   fdmrg::on                     = true;
    int    fdmrg::num_sites              = 30;
    int    fdmrg::max_sweeps             = 10;
    int    fdmrg::min_sweeps             = 4;
    long   fdmrg::chi_max                = 8;
    bool   fdmrg::chi_grow               = true;
    int    fdmrg::print_freq             = 100;
    int    fdmrg::store_freq             = 100;
    bool   fdmrg::store_wavefn           = false;                   /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */


    //Parameters controlling excited state DMRG
    bool   xdmrg::on                     = true;
    int    xdmrg::num_sites              = 200;
    int    xdmrg::max_sweeps             = 10;
    int    xdmrg::min_sweeps             = 4;
    long   xdmrg::chi_max                = 8;
    bool   xdmrg::chi_grow               = true;
    int    xdmrg::print_freq             = 100;
    int    xdmrg::store_freq             = 100;
    bool   xdmrg::store_wavefn           = false;                   /*!< Whether to store the wavefunction. Runs out of memory quick, recommended is false for max_length > 14 */
    double xdmrg::energy_density         = 0.5;                     /*!< Target energy in [0-1], where 0.5 means middle of spectrum. */
    double xdmrg::energy_window          = 0.01;                    /*!< Energy window [0-0.5] Accept states inside of energy_target +- energy_window. */

    //Parameters controlling imaginary TEBD (Zero temperature)
    bool   itebd::on                     = true;
    int    itebd::max_steps              = 100000;
    double itebd::delta_t0               = 0.1;
    double itebd::delta_tmin             = 0.00001;
    int    itebd::suzuki_order           = 1;
    long   itebd::chi_max                = 8;
    bool   itebd::chi_grow               = true;
    int    itebd::print_freq             = 5000;
    int    itebd::store_freq             = 100;

    //Save data_struct to hdf5
//    bool   hdf5::save_to_file             = true;
    bool   hdf5::save_progress            = true;
    string hdf5::output_filename          = "output/default.h5";
    string hdf5::access_mode              = "READWRITE";            /*!< Choose access mode to the file. Choose between READWRITE, READONLY */
    string hdf5::create_mode              = "RENAME" ;              /*!< Choose access mode to the file. Choose between TRUNCATE, OPEN, RENAME */
//    bool   hdf5::full_storage             = true;
    StorageLevel hdf5::storage_level            = StorageLevel::NORMAL;
    bool   hdf5::store_profiling          = true;

    //Profiling
    bool profiling::on                   = false;
    int  profiling::precision            = 5;
    //Console settings
    int  console::verbosity              = 2;
    bool console::timestamp              = false;

}







/*
 * The section below attempts to find and read the parameters from the given inputfile.
 *
*/


void settings::load_from_file(class_settings_reader &indata){
    input::input_filename   = indata.get_input_filename();
    input::input_file       = indata.get_input_file();
    //Parameters for the model Hamiltonian
    model::model_type                   = indata.find_parameter<std::string>("model::model_type"   , model::model_type);
    model::initial_state                = indata.find_parameter<std::string>("model::initial_state", model::initial_state);
    model::seed_init_mpo                = indata.find_parameter<int>        ("model::seed_init_mpo", model::seed_init_mpo);
    model::seed_init_mps                = indata.find_parameter<int>        ("model::seed_init_mps", model::seed_init_mps);

    model::tf_ising::J                  = indata.find_parameter<double> ("model::tf_ising::J"    , model::tf_ising::J);
    model::tf_ising::g                  = indata.find_parameter<double> ("model::tf_ising::g"    , model::tf_ising::g);
    model::tf_ising::w                  = indata.find_parameter<double> ("model::tf_ising::w"    , model::tf_ising::w);
    model::tf_ising::d                  = indata.find_parameter<int>    ("model::tf_ising::d"    , model::tf_ising::d);

    model::tf_nn_ising::J1              = indata.find_parameter<double> ("model::tf_nn_ising::J1", model::tf_nn_ising::J1);
    model::tf_nn_ising::J2              = indata.find_parameter<double> ("model::tf_nn_ising::J2", model::tf_nn_ising::J2);
    model::tf_nn_ising::g               = indata.find_parameter<double> ("model::tf_nn_ising::g" , model::tf_nn_ising::g);
    model::tf_nn_ising::d               = indata.find_parameter<int>    ("model::tf_nn_ising::d" , model::tf_nn_ising::d);
    model::tf_nn_ising::w               = indata.find_parameter<double> ("model::tf_nn_ising::w" , model::tf_nn_ising::w);

    model::selfdual_tf_rf_ising::J_log_mean  = indata.find_parameter<double> ("model::selfdual_tf_rf_ising::J_log_mean"    , model::selfdual_tf_rf_ising::J_log_mean);
    model::selfdual_tf_rf_ising::h_log_mean  = indata.find_parameter<double> ("model::selfdual_tf_rf_ising::h_log_mean"    , model::selfdual_tf_rf_ising::h_log_mean);
    model::selfdual_tf_rf_ising::J_sigma     = indata.find_parameter<double> ("model::selfdual_tf_rf_ising::J_sigma" , model::selfdual_tf_rf_ising::J_sigma);
    model::selfdual_tf_rf_ising::h_sigma     = indata.find_parameter<double> ("model::selfdual_tf_rf_ising::h_sigma" , model::selfdual_tf_rf_ising::h_sigma);
    model::selfdual_tf_rf_ising::lambda      = indata.find_parameter<double> ("model::selfdual_tf_rf_ising::lambda"  , model::selfdual_tf_rf_ising::lambda);
    model::selfdual_tf_rf_ising::d           = indata.find_parameter<int>    ("model::selfdual_tf_rf_ising::d"       , model::selfdual_tf_rf_ising::d);

    //Parmaters that control eigensolver and SVD precision
    precision::eigMaxIter                  = indata.find_parameter<int>    ("precision::eigMaxIter"  , precision::eigMaxIter);
    precision::eigThreshold                = indata.find_parameter<double> ("precision::eigThreshold", precision::eigThreshold);
    precision::eigMaxNcv                   = indata.find_parameter<int>    ("precision::eigMaxNcv"   , precision::eigMaxNcv);
    precision::SVDThreshold                = indata.find_parameter<double> ("precision::SVDThreshold", precision::SVDThreshold);
    precision::VarConvergenceThreshold     = indata.find_parameter<double> ("precision::VarConvergenceThreshold"   , precision::VarConvergenceThreshold);
    precision::VarSaturationThreshold      = indata.find_parameter<double> ("precision::VarSaturationThreshold"    , precision::VarSaturationThreshold);
    precision::EntEntrSaturationThreshold  = indata.find_parameter<double> ("precision::EntEntrSaturationThreshold", precision::EntEntrSaturationThreshold);
    precision::MaxSizeFullDiag             = indata.find_parameter<int>    ("precision::MaxSizeFullDiag"           , precision::MaxSizeFullDiag);

    //Parameters controlling infinite-DMRG
    idmrg::on                           = indata.find_parameter<bool>   ("idmrg::on"         , idmrg::on);
    if(idmrg::on){
        idmrg::max_steps                = indata.find_parameter<int>    ("idmrg::max_steps"  , idmrg::max_steps);
        idmrg::chi_max                  = indata.find_parameter<long>   ("idmrg::chi_max"    , idmrg::chi_max);
        idmrg::chi_grow                 = indata.find_parameter<bool>   ("idmrg::chi_grow"   , idmrg::chi_grow);
        idmrg::print_freq               = indata.find_parameter<int>    ("idmrg::print_freq" , idmrg::print_freq);
        idmrg::store_freq               = indata.find_parameter<int>    ("idmrg::store_freq" , idmrg::store_freq);
    }

    //Parameters controlling Finite-DMRG
    fdmrg::on                     = indata.find_parameter<bool>   ("fdmrg::on"         , fdmrg::on);
    if(fdmrg::on){
        fdmrg::num_sites          = indata.find_parameter<int>    ("fdmrg::num_sites ", fdmrg::num_sites);
        fdmrg::max_sweeps         = indata.find_parameter<int>    ("fdmrg::max_sweeps ", fdmrg::max_sweeps);
        fdmrg::min_sweeps         = indata.find_parameter<int>    ("fdmrg::min_sweeps ", fdmrg::min_sweeps);
        fdmrg::chi_max            = indata.find_parameter<int>    ("fdmrg::chi_max"    , 8);
        fdmrg::chi_grow           = indata.find_parameter<bool>   ("fdmrg::chi_grow"   , fdmrg::chi_grow);
        fdmrg::print_freq         = indata.find_parameter<int>    ("fdmrg::print_freq ", fdmrg::print_freq);
        fdmrg::store_freq         = indata.find_parameter<int>    ("fdmrg::store_freq ", fdmrg::store_freq);
        fdmrg::store_wavefn       = indata.find_parameter<bool>   ("fdmrg::store_wavefn" , fdmrg::store_wavefn);

    }

    //Parameters controlling excited state DMRG
    xdmrg::on                     = indata.find_parameter<bool>   ("xdmrg::on"         , xdmrg::on);
    if(xdmrg::on){
        xdmrg::num_sites          = indata.find_parameter<int>    ("xdmrg::num_sites "     , xdmrg::num_sites);
        xdmrg::max_sweeps         = indata.find_parameter<int>    ("xdmrg::max_sweeps "    , xdmrg::max_sweeps);
        xdmrg::min_sweeps         = indata.find_parameter<int>    ("xdmrg::min_sweeps "    , xdmrg::min_sweeps);
        xdmrg::chi_max            = indata.find_parameter<int>    ("xdmrg::chi_max"        , xdmrg::chi_max);
        xdmrg::chi_grow           = indata.find_parameter<bool>   ("xdmrg::chi_grow"       , xdmrg::chi_grow);
        xdmrg::print_freq         = indata.find_parameter<int>    ("xdmrg::print_freq "    , xdmrg::print_freq);
        xdmrg::store_freq         = indata.find_parameter<int>    ("xdmrg::store_freq "    , xdmrg::store_freq);
        xdmrg::store_wavefn       = indata.find_parameter<bool>   ("xdmrg::store_wavefn"   , xdmrg::store_wavefn);
        xdmrg::energy_density     = indata.find_parameter<double> ("xdmrg::energy_density" , xdmrg::energy_density);
        xdmrg::energy_window      = indata.find_parameter<double> ("xdmrg::energy_window"  , xdmrg::energy_window);
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
        itebd::store_freq         = indata.find_parameter<int>    ("itebd::store_freq"  , itebd::store_freq);
    }

    //Save data_struct to hdf5
//    hdf5::save_to_file             = indata.find_parameter<bool>   ("hdf5::save_to_file"            , hdf5::save_to_file           );
    hdf5::save_progress            = indata.find_parameter<bool>   ("hdf5::save_progress"           , hdf5::save_progress          );
    hdf5::output_filename          = indata.find_parameter<string> ("hdf5::output_filename"         , hdf5::output_filename);
    hdf5::access_mode              = indata.find_parameter<string> ("hdf5::access_mode"             , hdf5::access_mode);
    hdf5::create_mode              = indata.find_parameter<string> ("hdf5::create_mode"             , hdf5::create_mode);
//    hdf5::full_storage             = indata.find_parameter<bool>   ("hdf5::full_storage"            , hdf5::full_storage           );
    int storageLevelRead        = indata.find_parameter<int>  ("hdf5::storage_level"           , 2       );
    hdf5::storage_level            = static_cast<StorageLevel>(storageLevelRead );
    hdf5::store_profiling          = indata.find_parameter<bool>   ("hdf5::store_profiling"         , hdf5::store_profiling        );

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
