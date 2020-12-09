//
// Created by david on 7/30/17.
//
#include <algorithms/class_algorithm_launcher.h>
#include <algorithms/class_fdmrg.h>
#include <algorithms/class_flbit.h>
#include <algorithms/class_idmrg.h>
#include <algorithms/class_itebd.h>
#include <algorithms/class_xdmrg.h>
#include <gitversion.h>
#include <h5pp/h5pp.h>
#include <memory>
#include <config/nmspc_settings.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

//#include <stdlib.h>
#include <cstdlib>


class_algorithm_launcher::class_algorithm_launcher(std::shared_ptr<h5pp::File> h5ppFile_): h5pp_file(std::move(h5ppFile_)){
    tools::log = tools::Logger::setLogger("DMRG++ launch",  settings::console::verbosity, settings::console::timestamp);
    //Called in reverse order
    std::atexit(tools::common::profile::print_mem_usage);
    std::atexit(tools::common::profile::print_profiling_all);

    setup_temp_path();
}


class_algorithm_launcher::class_algorithm_launcher(){
    tools::log = tools::Logger::setLogger("DMRG++ launch", settings::console::verbosity, settings::console::timestamp);
    //Called in reverse order
    std::atexit(tools::common::profile::print_mem_usage);
    std::atexit(tools::common::profile::print_profiling_all);
    std::at_quick_exit(tools::common::profile::print_mem_usage);
    std::at_quick_exit(tools::common::profile::print_profiling_all);

    start_h5pp_file();
    setup_temp_path();
}


void class_algorithm_launcher::start_h5pp_file(){
    if( settings::output::storage_level_model      == StorageLevel::NONE and
        settings::output::storage_level_checkpoint == StorageLevel::NONE and
        settings::output::storage_level_good_state == StorageLevel::NONE and
        settings::output::storage_level_fail_state == StorageLevel::NONE and
        settings::output::storage_level_proj_state == StorageLevel::NONE and
        settings::output::storage_level_init_state == StorageLevel::NONE and
        settings::output::storage_level_emin_state == StorageLevel::NONE and
        settings::output::storage_level_emax_state == StorageLevel::NONE)
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
                    // Inspecting the file in READWRITE mode can update the file modification timestamp.
                    // Therefore we must first inspect the file in READONLY mode, then reopen in READWRITE.
                    // If the file is somehow damaged or invalid, we simply truncate it and start from scratch.
                    // If the file has a fully processed set of simulations, then [common/finished_all] = true,
                    // in which case we consult settings::output::file_resume_policy, to either exit or keep going
                    // and consult .cfg if there is anything more to be done.
                try{
                    h5pp_file = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::READONLY);
                    if (not h5pp_file->fileIsValid()) throw std::runtime_error(fmt::format("HDF5 file is not valid: {}",settings::output::output_filepath));
                    if (not h5pp_file->linkExists("common/finished_all"))
                        throw std::runtime_error(fmt::format("Could not find link common/finished_all in file: {}",settings::output::output_filepath));
                    if (not h5pp_file->linkExists("git/DMRG++"))
                        throw std::runtime_error(fmt::format("Could not find link git/DMRG++ in file: {}",settings::output::output_filepath));
                    if (settings::output::file_resume_policy == FileResumePolicy::FAST and h5pp_file->readDataset<bool>("common/finished_all")){
                        tools::log->info("Detected file_resume_policy == FileResumePolicy::FAST");
                        tools::log->info("Detected [common/finished_all] = true");
                        tools::log->info("All simulations have finished. Nothing more to do.");
                        exit(0);
                    }
                    h5pp_file = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::READWRITE);
                }catch (const std::exception & ex){
                    tools::log->error("Failed to resume simulation: {}", ex.what());
                    tools::log->info("Truncating file [{}]", settings::output::output_filepath);
                    h5pp_file = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::REPLACE);
                }
                break;
            }
            case FileCollisionPolicy::RENAME: {
                h5pp_file                = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::RENAME);
                std::string new_filepath = h5pp_file->getFilePath();
                tools::log->info("Renamed output file: [{}] -> [{}]", settings::output::output_filepath,new_filepath );
                settings::output::output_filepath = new_filepath;
                break;
            }
            case FileCollisionPolicy::BACKUP: {
                h5pp_file = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::BACKUP);
                tools::log->info("Renamed existing file: [{}] -> [{}]", settings::output::output_filepath,settings::output::output_filepath +".bak");
                break;
            }
            case FileCollisionPolicy::REPLACE: {
                h5pp_file = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::REPLACE);
                break;
            }
        }
    }else{
        h5pp_file = std::make_shared<h5pp::File>(settings::output::output_filepath,h5pp::FilePermission::COLLISION_FAIL);
    }
    h5pp_file->setCompressionLevel(settings::output::compression_level);
    if (not h5pp_file->linkExists("git/DMRG++")){
        //Put git metadata in file
        h5pp_file->writeDataset(GIT::BRANCH      , "git/DMRG++/branch");
        h5pp_file->writeDataset(GIT::COMMIT_HASH , "git/DMRG++/commit");
        h5pp_file->writeDataset(GIT::REVISION    , "git/DMRG++/revision");
    }

    if(not h5pp_file->linkExists("common")) {
        tools::log->trace("Copying config to h5pp file: {} --> {}", settings::input::config_filename, h5pp_file->getFileName());
        h5pp_file->writeDataset(settings::input::config_filename, "common/config_filename");
        h5pp_file->writeDataset(settings::input::config_file_contents, "common/config_file_contents");
    }
}



void class_algorithm_launcher::setup_temp_path(){
    if(not h5pp_file)
        return tools::log->warn("Can't set temporary path to a nullptr h5pp file");

    settings::output::tmp::hdf5_final_path = h5pp_file->getFilePath();
    if(not settings::output::use_temp_dir) return;
    if(not h5pp::fs::exists(settings::output::output_filepath))
        throw std::runtime_error("Can't set temporary path to non-existent file path: " + settings::output::output_filepath);

    tools::common::io::h5tmp::register_new_file(settings::output::output_filepath);
    settings::output::tmp::hdf5_final_path  = tools::common::io::h5tmp::get_original_filepath(h5pp_file->getFilePath());
    settings::output::tmp::hdf5_temp_path   = tools::common::io::h5tmp::get_temporary_filepath(h5pp_file->getFilePath());
    tools::log->info("Moving to temporary path [{}] --> [{}]",settings::output::tmp::hdf5_final_path , settings::output::tmp::hdf5_temp_path);
    h5pp_file->moveFileTo(settings::output::tmp::hdf5_temp_path, h5pp::FilePermission::REPLACE);
}



void class_algorithm_launcher::run_algorithms(){
    if(h5pp_file) h5pp_file->writeDataset(false, "common/finished_all");
    run_idmrg();
    run_fdmrg();
    run_flbit();
    run_xdmrg();
    run_itebd();

    if(h5pp_file) {
        h5pp_file->writeDataset(true, "common/finished_all");
        tools::log->info("Simulation data written to file: {}", settings::output::tmp::hdf5_final_path);
    }
    tools::log->info("All simulations finished");
//    h5pp_file->moveFile(tools::common::io::h5tmp::get_original_filepath(h5pp_file->getFilePath()), h5pp::FilePermission::REPLACE);
}


void class_algorithm_launcher::run_idmrg() const{
    if(settings::idmrg::on){
        class_idmrg idmrg(h5pp_file);
        idmrg.run();
    }
}


void class_algorithm_launcher::run_fdmrg() const{
    if(settings::fdmrg::on){
        class_fdmrg fdmrg(h5pp_file);
        fdmrg.run();
    }
}

void class_algorithm_launcher::run_flbit() const{
    if(settings::flbit::on){
        class_flbit flbit(h5pp_file);
        flbit.run();
    }
}


void class_algorithm_launcher::run_xdmrg() const{
    if(settings::xdmrg::on){
        class_xdmrg xdmrg(h5pp_file);
        xdmrg.run();
    }
}

void class_algorithm_launcher::run_itebd() const{
    if(settings::itebd::on){
        class_itebd itebd(h5pp_file);
        itebd.run();
    }
}

