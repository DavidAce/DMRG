//
// Created by david on 2018-01-14.
//

#include "n_sim_settings.h"
//Parmaters that control eigensolver and SVD precision
using namespace std;

namespace settings{

int    precision::eigSteps           = 5000;
double precision::eigThreshold       = 1e-12;
int    precision::eig_max_ncv        = 20;
double precision::SVDThreshold       = 1e-12;

//Parameters controlling infinite-DMRG
bool   idmrg::on                     = false;
int    idmrg::max_length             = 200;
long   idmrg::chi_max                = 25;
//Parameters controlling Finite-DMRG
bool   fdmrg::on                     = false;
int    fdmrg::max_length             = 200;
int    fdmrg::max_sweeps             = 2;
long   fdmrg::chi_max                = 8;
//Parameters controlling imaginary TEBD (Zero temperature)
bool   itebd::on                     = false;
int    itebd::max_steps              = 100000;
double itebd::delta_t                = 0.01;
long   itebd::chi_max                = 25;

//Parameters controlling Finite-entanglement scaling (FES) in iTEBD-mode.
bool   fes_itebd::on                 = false;
int    fes_itebd::max_steps          = 100000;
double fes_itebd::delta_t            = 0.01;
long   fes_itebd::chi_min            = 4;
long   fes_itebd::chi_max            = 12;
long   fes_itebd::chi_num            = 3;

//Parameters controlling Finite-entanglement scaling (FES) in iDMRG-mode.
bool   fes_idmrg::on                 = false;
int    fes_idmrg::max_steps          = 4000;
long   fes_idmrg::chi_min            = 4;
long   fes_idmrg::chi_max            = 12;
long   fes_idmrg::chi_num            = 3;

//Save data to hdf5
bool   hdf5::save_to_file             = true;
bool   hdf5::create_dir_if_not_found  = true;
string hdf5::filename                = "data.h5";
string hdf5::path                    = "../output";
bool   hdf5::full_storage             = true;

//Profiling
bool profiling::on                   = false;
int  profiling::precision            = 5;
//Console settings
int  console::verbosity              = 0;

}



void settings::initialize(class_file_reader &indata){
    //Parmaters that control eigensolver and SVD precision
    precision::eigSteps           = indata.find_parameter<int>    ("precision::eigSteps"    , precision::eigSteps);
    precision::eigThreshold       = indata.find_parameter<double> ("precision::eigThreshold", precision::eigThreshold);
    precision::eig_max_ncv        = indata.find_parameter<int>    ("precision::eig_max_ncv" , precision::eig_max_ncv);
    precision::SVDThreshold       = indata.find_parameter<double> ("precision::eigThreshold", precision::SVDThreshold);
    //Parameters controlling infinite-DMRG
    idmrg::on                     = indata.find_parameter<bool>   ("idmrg::on"         , idmrg::on);
    idmrg::max_length             = indata.find_parameter<int>    ("idmrg::max_length ", idmrg::max_length);
    idmrg::chi_max                = indata.find_parameter<long>   ("idmrg::chi_max"    , idmrg::chi_max);
    //Parameters controlling Finite-DMRG
    fdmrg::on                     = indata.find_parameter<bool>   ("fdmrg::on"         , fdmrg::on);
    fdmrg::max_length             = indata.find_parameter<int>    ("fdmrg::max_length ", fdmrg::max_length);
    fdmrg::max_sweeps             = indata.find_parameter<int>    ("fdmrg::max_sweeps ", fdmrg::max_sweeps);
    fdmrg::chi_max                = indata.find_parameter<int>    ("fdmrg::chi_max"    , 8);
    //Parameters controlling imaginary TEBD (Zero temperature)
    itebd::on                     = indata.find_parameter<bool>   ("itebd::on"         , itebd::on       );
    itebd::max_steps              = indata.find_parameter<int>    ("itebd::max_steps " , itebd::max_steps);
    itebd::delta_t                = indata.find_parameter<double> ("itebd::delta_t "   , itebd::delta_t  );
    itebd::chi_max                = indata.find_parameter<long>   ("itebd::chi_max"    , itebd::chi_max  );

    //Parameters controlling Finite-entanglement scaling (FES) in iTEBD-mode.
    fes_itebd::on                 = indata.find_parameter<bool>   ("fes_itebd::on"         , fes_itebd::on        );
    fes_itebd::max_steps          = indata.find_parameter<int>    ("fes_itebd::max_steps " , fes_itebd::max_steps );
    fes_itebd::delta_t            = indata.find_parameter<double> ("fes_itebd::delta_t "   , fes_itebd::delta_t   );
    fes_itebd::chi_min            = indata.find_parameter<long>   ("fes_itebd::chi_min"    , fes_itebd::chi_min   );
    fes_itebd::chi_max            = indata.find_parameter<long>   ("fes_itebd::chi_max"    , fes_itebd::chi_max   );
    fes_itebd::chi_num            = indata.find_parameter<long>   ("fes_itebd::chi_num"    , fes_itebd::chi_num   );

    //Parameters controlling Finite-entanglement scaling (FES) in iDMRG-mode.
    fes_idmrg::on                 = indata.find_parameter<bool>   ("fes_idmrg::on"         , fes_idmrg::on       );
    fes_idmrg::max_steps          = indata.find_parameter<int>    ("fes_idmrg::max_steps " , fes_idmrg::max_steps);
    fes_idmrg::chi_min            = indata.find_parameter<long>   ("fes_idmrg::chi_min"    , fes_idmrg::chi_min  );
    fes_idmrg::chi_max            = indata.find_parameter<long>   ("fes_idmrg::chi_max"    , fes_idmrg::chi_max  );
    fes_idmrg::chi_num            = indata.find_parameter<long>   ("fes_idmrg::chi_num"    , fes_idmrg::chi_num  );

    //Save data to hdf5
    hdf5::save_to_file             = indata.find_parameter<bool>   ("hdf5::save_to_file"            , hdf5::save_to_file           );
    hdf5::create_dir_if_not_found  = indata.find_parameter<bool>   ("hdf5::create_dir_if_not_found" , hdf5::create_dir_if_not_found);
    hdf5::filename                 = indata.find_parameter<string> ("hdf5::filename"                , hdf5::filename               );
    hdf5::path                     = indata.find_parameter<string> ("hdf5::path"                    , hdf5::path                   );
    hdf5::full_storage             = indata.find_parameter<bool>   ("hdf5::full_storage"            , hdf5::full_storage           );

    //Profiling
    profiling::on                  = indata.find_parameter<bool>   ("profiling::on"        , profiling::on        );
    profiling::precision           = indata.find_parameter<int>    ("profiling::precision" , profiling::precision );
    //Console settings
    console::verbosity             = indata.find_parameter<int>    ("console::verbosity"          , console::verbosity);
}
