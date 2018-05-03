//
// Created by david on 7/30/17.
//
#include <sim_parameters/nmspc_sim_settings.h>
#include <algorithms/class_algorithm_launcher.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_finite_chain_storage.h>
#include <mps_routines/class_measurement.h>
#include <IO/class_hdf5_file.h>
#include <IO/class_hdf5_table_buffer.h>

#include <algorithms/class_iDMRG.h>
#include <algorithms/class_fDMRG.h>
#include <algorithms/class_iTEBD.h>
#include <algorithms/class_FES_iDMRG.h>
#include <algorithms/class_FES_iTEBD.h>
#include "class_xDMRG.h"


namespace s = settings;
using namespace std;
using namespace Textra;

class_algorithm_launcher::class_algorithm_launcher(std::shared_ptr<class_hdf5_file> hdf5_):hdf5(std::move(hdf5_))
{

};
class_algorithm_launcher::class_algorithm_launcher()
{
    hdf5 = std::make_shared<class_hdf5_file>(settings::hdf5::output_filename, settings::hdf5::output_folder,true, false);
};


void class_algorithm_launcher::run_infinite_DMRG(){
    if(settings::idmrg::on){
        class_iDMRG iDMRG(hdf5);
        iDMRG.run();
    }
}


void class_algorithm_launcher::run_finite_DMRG(){
    if(settings::fdmrg::on){
        class_fDMRG fDMRG(hdf5);
        fDMRG.run();
    }
}

void class_algorithm_launcher::run_excited_state_DMRG(){
    if(settings::xdmrg::on){
        class_xDMRG xDMRG(hdf5);
        xDMRG.run();
    }
}

void class_algorithm_launcher::run_imaginary_TEBD(){
    if(settings::itebd::on){
        class_iTEBD iTEBD(hdf5);
        iTEBD.run();
    }
}


void class_algorithm_launcher::run_FES_iDMRG(){
    if(settings::fes_idmrg::on){
        class_FES_iDMRG FES_iDMRG(hdf5);
        FES_iDMRG.run();
    }
}

void class_algorithm_launcher::run_FES_iTEBD(){
    if(settings::fes_itebd::on){
        class_FES_iTEBD FES_iTEBD(hdf5);
        FES_iTEBD.run();
    }
}
