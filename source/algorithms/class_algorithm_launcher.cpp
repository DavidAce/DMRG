//
// Created by david on 7/30/17.
//
#include <sim_parameters/nmspc_sim_settings.h>
#include <algorithms/class_algorithm_launcher.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_finite_chain_state.h>
#include <io/class_hdf5_table_buffer2.h>
#include <algorithms/class_iDMRG.h>
#include <algorithms/class_fDMRG.h>
#include <algorithms/class_xDMRG.h>
#include <algorithms/class_iTEBD.h>
#include <spdlog/spdlog.h>
#include <h5pp/h5pp.h>

namespace s = settings;
using namespace std;

class_algorithm_launcher::class_algorithm_launcher(std::shared_ptr<h5pp::File> h5ppFile_):h5ppFile(std::move(h5ppFile_))
{
    hdf5_path = h5ppFile->get_file_path();

}
class_algorithm_launcher::class_algorithm_launcher()
{

    h5pp::CreateMode createMode;
    if(settings::hdf5::create_mode == "TRUNCATE") createMode        = h5pp::CreateMode::TRUNCATE;
    else if (settings::hdf5::create_mode == "RENAME") createMode    = h5pp::CreateMode::RENAME;
    else if (settings::hdf5::create_mode == "OPEN") createMode      = h5pp::CreateMode::OPEN;
    else {throw std::runtime_error("Wrong create mode: " + settings::hdf5::create_mode);}


    h5ppFile = std::make_shared<h5pp::File>(
            settings::hdf5::output_filename,
            h5pp::AccessMode::READWRITE,
            createMode);
    hdf5_path = h5ppFile->get_file_path();
}


void class_algorithm_launcher::run_algorithms(){
    run_iDMRG();
    run_fDMRG();
    run_xDMRG();
    run_iTEBD();
    spdlog::info("All simulations finished");
    bool OK = true;
    h5ppFile->writeDataset(OK, "/common/fileOK");
    spdlog::info("Simulation data written to file: {}", h5ppFile->get_file_path());
}


void class_algorithm_launcher::run_iDMRG(){
    if(settings::idmrg::on){
        class_iDMRG iDMRG(h5ppFile);
        iDMRG.run();
    }
}


void class_algorithm_launcher::run_fDMRG(){
    if(settings::fdmrg::on){
        class_fDMRG fDMRG(h5ppFile);
        fDMRG.run();
    }
}

void class_algorithm_launcher::run_xDMRG(){
    if(settings::xdmrg::on){
        class_xDMRG xDMRG(h5ppFile);
        xDMRG.run();
    }
}

void class_algorithm_launcher::run_iTEBD(){
    if(settings::itebd::on){
        class_iTEBD iTEBD(h5ppFile);
        iTEBD.run();
    }
}

