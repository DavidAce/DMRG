//
// Created by david on 2019-11-07.
//
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <simulation/nmspc_settings.h>
#include <h5pp/h5pp.h>
#include <io/nmspc_filesystem.h>
#include <fstream>
#include <stdlib.h>


std::string get_dirname(){
    return "DMRG." + std::string(getenv("USER")) + "/";
}


std::string tools::common::io::h5tmp::set_tmp_prefix(const std::string &output_filename) {
    fs::path temp_path;
    if(fs::exists(settings::output::temp_dir))
         temp_path = settings::output::temp_dir/ fs::path(get_dirname());
    else temp_path = fs::temp_directory_path() / fs::path(get_dirname());

    std::string::size_type pos = output_filename.find(temp_path.string());
    if (pos != std::string::npos)
        return output_filename;
    else
        return fs::absolute(temp_path / output_filename);
}

std::string tools::common::io::h5tmp::unset_tmp_prefix(const std::string &output_filename) {
    fs::path temp_path;
    if(fs::exists(settings::output::temp_dir))
        temp_path = settings::output::temp_dir / fs::path(get_dirname());
    else if (fs::exists("/dev/shm"))
        temp_path = "/dev/shm" / fs::path(get_dirname());
    else if (fs::exists("/scratch/local"))
        temp_path = "/scratch/local" / fs::path(get_dirname());
    else temp_path = fs::temp_directory_path() / fs::path(get_dirname());

    std::string::size_type pos = output_filename.find(temp_path.string());
    if (pos != std::string::npos){
        std::string new_filename = output_filename;
        new_filename.erase(pos, temp_path.string().length());
        return fs::current_path() / fs::path(new_filename);
    }else{
        return output_filename;
    }

}

void tools::common::io::h5tmp::create_directory(const std::string & path){
    if(path.empty()) return;
    if(tools::log == nullptr){
        if(spdlog::get("DMRG") == nullptr)
             tools::log = spdlog::default_logger();
        else tools::log = spdlog::get("DMRG");
    }
    fs::path dir = fs::absolute(path);
    if(dir.has_filename() and dir.has_extension()) dir = dir.parent_path();

    try{
        if (fs::create_directories(dir)){
            tools::log->info("Created directory: {}",dir.string());
        }else{
            tools::log->info("Directory already exists: {}",dir.string());
        }
    }
    catch(std::exception & ex){
        throw std::runtime_error("Failed to create directory: " + std::string(ex.what()));
    }
}


void tools::common::io::h5tmp::copy_from_tmp(const std::string &output_filename){
    if(output_filename.empty()) return;
    copy_file(output_filename,unset_tmp_prefix(output_filename));
}


void tools::common::io::h5tmp::copy_file(const std::string & src, const std::string & tgt ){
    if(src == tgt) return;
    fs::path target_path = src;
    fs::path source_path = tgt;

    if(not fs::exists(target_path.parent_path())){
        tools::common::io::h5tmp::create_directory(target_path);
    }
    if(fs::exists(target_path)){
        std::ifstream target_stream(target_path.string(), std::ios_base::binary);
        std::ifstream source_stream(source_path.string(), std::ios_base::binary);
        typedef std::istreambuf_iterator<char> isbuf_it;
        if (std::equal(isbuf_it(target_stream.rdbuf()), isbuf_it(),
                       isbuf_it(source_stream.rdbuf()), isbuf_it()))
        {
            tools::log->debug("Source [{}] identical to target [{}] ... Skipping copy", src,tgt);
            return;
        }
    }
    tools::log->debug("Copying hdf5 file: {} -> {}",src,tgt);
    fs::copy(source_path, target_path, fs::copy_options::update_existing);
}



void tools::common::io::h5tmp::remove_from_temp(const std::string output_filename){
    if(output_filename.empty()) {std::cout << "Nothing to delete" << std::endl << std::flush; return;}
    fs::path temp_path;
    if(fs::exists(settings::output::temp_dir))
        temp_path = settings::output::temp_dir / fs::path(get_dirname());
    else if (fs::exists("/dev/shm"))
        temp_path = "/dev/shm" / fs::path(get_dirname());
    else if (fs::exists("/scratch/local"))
        temp_path = "/scratch/local" / fs::path(get_dirname());
    else temp_path = fs::temp_directory_path() / fs::path(get_dirname());


    std::string::size_type pos = output_filename.find(temp_path.string());
    if (pos != std::string::npos){
        // Path points to the temp directory!
        if(fs::exists(output_filename)){
            tools::log->debug("Deleting temporary file: {}",output_filename);
            fs::remove(output_filename);
        }else{
            tools::log->debug("Nothing to delete");
        }
    }else{
        tools::log->debug("Temp file is disabled - nothing to delete");
    }
}
