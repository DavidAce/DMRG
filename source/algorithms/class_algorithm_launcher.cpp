//
// Created by david on 7/30/17.
//
#include <simulation/nmspc_settings.h>
#include <algorithms/class_algorithm_launcher.h>
#include <state/class_state_infinite.h>
#include <state/class_state_finite.h>
#include <algorithms/class_iDMRG.h>
#include <algorithms/class_fDMRG.h>
#include <algorithms/class_xDMRG.h>
#include <algorithms/class_iTEBD.h>
#include <tools/common/log.h>
#include <tools/common/io.h>
#include <tools/common/prof.h>
#include <h5pp/h5pp.h>
#include <memory>
#include <gitversion.h>

//#include <stdlib.h>
#include <cstdlib>

namespace s = settings;
using namespace std;

static std::string hdf5_temp_path;
static std::string hdf5_final_path;

class_algorithm_launcher::class_algorithm_launcher(std::shared_ptr<h5pp::File> h5ppFile_):h5ppFile(std::move(h5ppFile_)){
    setLogger("DMRG");
    setup_temp_path();
    //Called in reverse order
    std::atexit(tools::common::profile::print_mem_usage);
    std::atexit(tools::common::profile::print_profiling);
}


class_algorithm_launcher::class_algorithm_launcher(){
    setLogger("DMRG");
    start_h5pp_file();
    setup_temp_path();
    //Called in reverse order
    std::atexit(tools::common::profile::print_mem_usage);
    std::atexit(tools::common::profile::print_profiling);
    std::at_quick_exit(tools::common::profile::print_mem_usage);
    std::at_quick_exit(tools::common::profile::print_profiling);
}



void class_algorithm_launcher::clean_up() {
    H5garbage_collect();
    H5Eprint(H5E_DEFAULT, stderr);
    H5close();
    std::this_thread::sleep_until(std::chrono::system_clock::now() + std::chrono::seconds(3));
    tools::common::io::h5tmp::copy_from_tmp(hdf5_temp_path);
    tools::common::io::h5tmp::remove_from_tmp(hdf5_temp_path);
}


void class_algorithm_launcher::setLogger(const std::string& name){
    if(spdlog::get(name) == nullptr){
        log = spdlog::stdout_color_mt(name);
        log->set_pattern("[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v");
        log->set_level(spdlog::level::trace);
    }else
        log = spdlog::get(name);
}


void class_algorithm_launcher::start_h5pp_file(){
    if( settings::output::storage_level_journal == StorageLevel::NONE and
        settings::output::storage_level_results == StorageLevel::NONE and
        settings::output::storage_level_chi_update == StorageLevel::NONE and
        settings::output::storage_level_proj_state == StorageLevel::NONE)
        return;

    // There are two possibilities depending on settings::output::output_filename
    // 1) The .h5 file exists
    //    It is important to honor settings::output::file_policy
    //      A) Load simulation state from existing .h5, resume the simulation and continue appending data to the same file.
    //      B) Load simulation state from existing .h5, resume the simulation but add new data to a new file.
    //      C) Remove the .h5 file and start from scratch (i.e. go to 2))

    // 2) The .h5 file does not exist yet -> start new simulation


    if(h5pp::fs::exists(settings::output::output_filepath)){
        switch(settings::output::file_collision_policy){
            case FileCollisionPolicy::RESUME: {
                h5ppFile = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::READWRITE);
                break;
            }
            case FileCollisionPolicy::RENAME: {
                h5ppFile = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::RENAME);
                std::string new_filepath = h5ppFile->getFilePath();
                log->info("Renamed output file: [{}] -> [{}]", settings::output::output_filepath,new_filepath );
                settings::output::output_filepath = new_filepath;
                break;
            }
            case FileCollisionPolicy::BACKUP: {
                h5ppFile = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::BACKUP);
                log->info("Renamed existing file: [{}] -> [{}]", settings::output::output_filepath,settings::output::output_filepath +".bak");
                break;
            }
            case FileCollisionPolicy::REPLACE: {
                h5ppFile = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::REPLACE);
                break;
            }
        }
    }else{
        h5ppFile = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::COLLISION_FAIL);
    }
    h5ppFile->setCompressionLevel(settings::output::compression_level);
    if (h5ppFile->getFilePermission() != h5pp::FilePermission::READWRITE){
        //Put git revision in file attribute
        h5ppFile->writeAttribute(GIT::BRANCH      , "GIT BRANCH", "/");
        h5ppFile->writeAttribute(GIT::COMMIT_HASH , "GIT COMMIT", "/");
        h5ppFile->writeAttribute(GIT::REVISION    , "GIT REVISION", "/");
    }
}



void class_algorithm_launcher::setup_temp_path(){
    if(not h5ppFile)
        return log->warn("Can't set temporary path to a nullptr h5pp file");

    hdf5_temp_path  = h5ppFile->getFilePath();
    hdf5_final_path = h5ppFile->getFilePath();
    if(not settings::output::use_temp_dir) return;
    if(not h5pp::fs::exists(settings::output::output_filepath))
        throw std::runtime_error("Can't set temporary path to non-existent file path: " + settings::output::output_filepath);
    h5ppFile->flush();
    tools::common::io::h5tmp::register_new_file(settings::output::output_filepath);
    tools::common::io::h5tmp::copy_into_tmp(settings::output::output_filepath);
    auto & temp_filepath = tools::common::io::h5tmp::get_temporary_filepath(settings::output::output_filepath);
    h5ppFile = std::make_shared<h5pp::File>(temp_filepath, h5pp::FilePermission::READWRITE);
//    h5ppFile->setLogLevel(0);
    std::at_quick_exit(class_algorithm_launcher::clean_up);
    std::atexit(class_algorithm_launcher::clean_up);
}



void class_algorithm_launcher::run_algorithms(){
    if(h5ppFile) h5ppFile->writeDataset(false, "common/finished_all");
    run_iDMRG();
    run_fDMRG();
    run_xDMRG();
    run_iTEBD();

    if(h5ppFile) {
        h5ppFile->writeDataset(true, "common/finished_all");
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

