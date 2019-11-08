//
// Created by david on 2019-11-07.
//
#include <simulation/nmspc_settings.h>
#include <tools/nmspc_tools.h>
#include <h5pp/h5pp.h>
#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;

std::string tools::common::io::h5tmp::set_tmp_prefix(const std::string &output_filename) {
    fs::path temp_path;
    if(fs::exists(settings::output::temp_dir))
         temp_path = settings::output::temp_dir/ fs::path("DMRG/");
    else temp_path = fs::temp_directory_path() / fs::path("DMRG/");


    std::string::size_type pos = output_filename.find(temp_path.string());
    if (pos != std::string::npos){
        return output_filename;
    }else{
        fs::path h5pp_path = fs::absolute(output_filename);
        return fs::absolute(temp_path / fs::relative(h5pp_path,fs::current_path()));
    }


}

std::string tools::common::io::h5tmp::unset_tmp_prefix(const std::string &output_filename) {
    fs::path temp_path;
    if(fs::exists(settings::output::temp_dir))
         temp_path = settings::output::temp_dir/ fs::path("DMRG/");
    else temp_path = fs::temp_directory_path() / fs::path("DMRG/");

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
    if(tools::log == nullptr){
        if(spdlog::get("DMRG") == nullptr)
             tools::log = spdlog::default_logger();
        else tools::log = spdlog::get("DMRG");
    }
    fs::path dir = fs::absolute(path);
    if(dir.has_filename()) dir = dir.parent_path();

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
    fs::path target_path = unset_tmp_prefix(output_filename);
    fs::path source_path = output_filename;

    if(target_path == source_path) return;

    if(not fs::exists(target_path.parent_path())){
        tools::common::io::h5tmp::create_directory(target_path.parent_path());
    }
    if(fs::exists(target_path)){
        std::ifstream target_stream(target_path.string(), std::ios_base::binary);
        std::ifstream source_stream(source_path.string(), std::ios_base::binary);
        typedef std::istreambuf_iterator<char> isbuf_it;
        if (std::equal(isbuf_it(target_stream.rdbuf()), isbuf_it(),
                       isbuf_it(source_stream.rdbuf()), isbuf_it()))
        {
            tools::log->debug("Source and target files are equal... Skipping copy");
            return;
        }
    }
    tools::log->debug("Copying hdf5 file to target path: {} -> {}",source_path.string(), target_path.string());
    fs::copy(source_path, target_path, fs::copy_options::update_existing);



}

void tools::common::io::h5tmp::remove_from_temp(const std::string & output_filename){
    fs::path temp_path;
    if(fs::exists(settings::output::temp_dir))
        temp_path = settings::output::temp_dir/ fs::path("DMRG/");
    else temp_path = fs::temp_directory_path() / fs::path("DMRG/");

    std::string::size_type pos = output_filename.find(temp_path.string());
    if (pos != std::string::npos){
        // Its in the temp directory!
        copy_from_tmp(output_filename);
        tools::log->debug("Deleting temporary file: {}",output_filename);
        fs::remove(output_filename);
    }
}
