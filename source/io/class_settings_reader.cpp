//
// Created by david on 2018-01-12.
//

#include <spdlog/spdlog.h>
#include "class_settings_reader.h"


class_settings_reader::class_settings_reader(const fs::path &file_path_): file_path(file_path_) {
//    auto logger_exists = spdlog::get("DMRG");
//    if (logger_exists) {
//        spdlog::set_default_logger(logger_exists);
//        spdlog::debug("Previous logger exists");
//    }else{
//        spdlog::warn("Previous logger does not exist");
//    }

    fs::path found_file =  find_input_file(file_path);
    if (not found_file.empty()){
        try {
            file.open(found_file.c_str());
        }
        catch(std::exception &ex){
            spdlog::critical("Exiting: {}", ex.what() );
            exit(1);
        }
    }
}


bool class_settings_reader::check_if_input_file_exists(const fs::path &path_to_file){
    if (path_to_file.has_filename()){
        if(fs::exists(path_to_file)){
            std::ifstream in(path_to_file.c_str());
            if(in.is_open()){
                in.close();
                return true;
            }
        }
        spdlog::debug("File does not exist: {}",path_to_file.string());
        return false;
    }
    spdlog::debug("Given output_folder does not point to a file: {}", path_to_file.string());
    return false;
}

fs::path class_settings_reader::find_input_file(const fs::path &given_path) {

    //Check if file exists in the given path.
    fs::path complete_path = fs::system_complete(given_path);

    spdlog::debug("Checking for input file: [ {} ] in path: [ {} ]", given_path.string() ,complete_path.string());
    if (check_if_input_file_exists(complete_path)){
        spdlog::info("Found input file: [ {} ] in path: [ {} ]", given_path.string() ,fs::canonical(complete_path).string());
        found_file = true;
        return fs::canonical(complete_path);
    }

    //Check if file exists in the given path (if it is a relative path!), relative to the executable.
    if (given_path.is_relative()) {
        complete_path = fs::system_complete(fs::current_path() / given_path);
        spdlog::debug("Checking for input file: [ {} ] in path: [ {} ]", given_path.string() ,complete_path.string());
        if (check_if_input_file_exists(complete_path)) {
            spdlog::info("Found input file: [ {} ] in path: [ {} ]", given_path.string() ,fs::canonical(complete_path).string());
            found_file = true;
            return fs::canonical(complete_path);
        }
    }

    //Check if file exists in current directory
    complete_path = fs::system_complete(fs::current_path()/given_path.filename());
    spdlog::debug("Checking for input file: [ {} ] in path: [ {} ]", given_path.string() ,complete_path.string());
    if(check_if_input_file_exists(complete_path)){
        spdlog::info("Found input file: [ {} ] in path: [ {} ]", given_path.string() ,fs::canonical(complete_path).string());
        found_file = true;
        return fs::canonical(complete_path);
    }
    spdlog::warn("Input file could not be found: [ {} ]", given_path.string());
    found_file = false;
    return fs::path();
}

void class_settings_reader::remove_spaces(std::string &str){
    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
}

bool class_settings_reader::has_only_digits(const std::string s){
    return s.find_first_not_of( "+-0123456789" ) == std::string::npos;
}

bool class_settings_reader::is_parameterline(const std::string s){
    return s.find("=") != std::string::npos;
}

std::string::size_type class_settings_reader::find_comment_character(const std::string s){
    std::vector<std::string> comment_symbols = {"//", "/*", "#"};
    for(auto &sym : comment_symbols){
        if(s.find(sym) != std::string::npos){
            return s.find(sym);
        }
    }
    return s.npos;
}

