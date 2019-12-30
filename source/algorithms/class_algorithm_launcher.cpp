//
// Created by david on 7/30/17.
//
#include <simulation/nmspc_settings.h>
#include <algorithms/class_algorithm_launcher.h>
#include <state/class_state_infinite.h>
#include <state/class_state_finite.h>
#include <io/class_h5table_buffer.h>
#include <algorithms/class_iDMRG.h>
#include <algorithms/class_fDMRG.h>
#include <algorithms/class_xDMRG.h>
#include <algorithms/class_iTEBD.h>
#include <spdlog/spdlog.h>
#include <h5pp/h5pp.h>
#include <gitversion.h>
#include <tools/nmspc_tools.h>

//#include <stdlib.h>
#include <cstdlib>

namespace s = settings;
using namespace std;

static std::string hdf5_temp_path;
static std::string hdf5_final_path;

void class_algorithm_launcher::remove_temp_file() {
    tools::common::io::h5tmp::remove_from_temp(hdf5_temp_path);
}



void class_algorithm_launcher::setLogger(std::string name){
    if(spdlog::get(name) == nullptr){
        log = spdlog::stdout_color_mt(name);
        log->set_pattern("[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v");
        log->set_level(spdlog::level::trace);
    }else{
        log = spdlog::get(name);
    }

}



class_algorithm_launcher::class_algorithm_launcher(std::shared_ptr<h5pp::File> h5ppFile_):h5ppFile(std::move(h5ppFile_))
{

    setLogger("DMRG");
    hdf5_temp_path  = h5ppFile->getFilePath();
    hdf5_final_path = h5ppFile->getFilePath();

}
class_algorithm_launcher::class_algorithm_launcher()
{
    setLogger("DMRG");
    std::at_quick_exit(class_algorithm_launcher::remove_temp_file);

    if (settings::output::storage_level == StorageLevel::NONE){return;}

    h5pp::CreateMode createMode;
    if(settings::output::create_mode == "TRUNCATE") createMode        = h5pp::CreateMode::TRUNCATE;
    else if (settings::output::create_mode == "RENAME") createMode    = h5pp::CreateMode::RENAME;
    else if (settings::output::create_mode == "OPEN") createMode      = h5pp::CreateMode::OPEN;
    else {throw std::runtime_error("Wrong create mode: " + settings::output::create_mode);}



    if(settings::output::use_temp_dir){
        hdf5_temp_path = tools::common::io::h5tmp::set_tmp_prefix(settings::output::output_filename);
        tools::common::io::h5tmp::create_directory(hdf5_temp_path);
        h5ppFile = std::make_shared<h5pp::File>(hdf5_temp_path,h5pp::AccessMode::READWRITE,createMode);
    }else{
        tools::common::io::h5tmp::create_directory(settings::output::output_filename);
        h5ppFile = std::make_shared<h5pp::File>(settings::output::output_filename,h5pp::AccessMode::READWRITE,createMode);
    }

    hdf5_final_path = tools::common::io::h5tmp::unset_tmp_prefix(h5ppFile->getFilePath());


    if (createMode == h5pp::CreateMode::TRUNCATE or createMode == h5pp::CreateMode::RENAME){
        //Put git revision in file attribute
        h5ppFile->writeAttributeToFile(GIT::BRANCH      , "GIT BRANCH");
        h5ppFile->writeAttributeToFile(GIT::COMMIT_HASH , "GIT COMMIT");
        h5ppFile->writeAttributeToFile(GIT::REVISION    , "GIT REVISION");
    }


}


void class_algorithm_launcher::run_algorithms(){
    if(h5ppFile) h5ppFile->writeDataset(false, "/common/finOK");
    run_iDMRG();
    run_fDMRG();
    run_xDMRG();
    run_iTEBD();

    if(h5ppFile) {
        h5ppFile->writeDataset(true, "/common/finOK");
        tools::common::io::h5tmp::copy_from_tmp(hdf5_temp_path);
        h5ppFile.reset();// Kill the file
        tools::common::io::h5tmp::remove_from_temp(hdf5_temp_path);
        log->info("Simulation data written to file: {}", hdf5_final_path);
    }
    log->info("All simulations finished");
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

