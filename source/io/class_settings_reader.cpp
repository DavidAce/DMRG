//
// Created by david on 2018-01-12.
//

#include "class_settings_reader.h"
#include <io/nmspc_logger.h>


class_settings_reader::class_settings_reader(const fs::path &file_path_,std::string logName): file_path(file_path_) {
    log = Logger::setLogger(logName,2,true);
    fs::path found_file =  find_input_file(file_path);
    std::ifstream          file;
    if (not found_file.empty()){
        try {
            file.open(found_file.c_str());
        }
        catch(std::exception &ex){
            log->critical("Could not open file [ {} ]: {}", found_file.string(), std::string(ex.what()) );
            throw std::runtime_error("Could not open file [ " + found_file.string() + "]: " + std::string(ex.what()));
        }
    }

    // Load the file into memory
    if (file.is_open()) {
        file.clear();
        file.seekg(0, std::ios::beg);
        file_string = std::string ( (std::istreambuf_iterator<char>(file) ),
                                    (std::istreambuf_iterator<char>()     ));

    }
    file.close();

    // Extract the parameters and generate a key-value map
    std::istringstream file_stream (file_string);
    std::string line;
    while (std::getline(file_stream, line)){
        if(is_parameterline(line)){
            log->debug("Loading line: {}",line);
            std::string param_key;
            std::string param_val;
            std::istringstream linestream(line);
            std::getline(linestream,param_key, '=');
            std::getline(linestream,param_val,  '\n');
            param_val = param_val.substr(0, find_comment_character(param_val));
            remove_spaces(param_key);
            remove_spaces(param_val);
            std::transform(param_key.begin(), param_key.end(), param_key.begin(), ::tolower);
            param_map[param_key] = param_val;
        }
    }


    //Generate a parameter  key-value map
//    std::vector<std::istringstream> file_cache;
//
//    for (auto & linestream: file_cache){
//        std::string param_key;
//        std::string param_val;
////        std::istringstream is(line);
//        std::getline(linestream,param_key, '=');
//        std::getline(linestream,param_val,  '\n');
//        param_val = param_val.substr(0, find_comment_character(param_val));
//        remove_spaces(param_key);
//        remove_spaces(param_val);
//        std::transform(param_key.begin(), param_key.end(), param_key.begin(), ::tolower);
//        param_map[param_key] = param_val;
//    }

}

std::string class_settings_reader::get_input_file_as_string(){
    return file_string;
}

std::string class_settings_reader::get_input_filename(){
    return file_path.string();
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
        log->debug("File does not exist: {}",path_to_file.string());
        return false;
    }
    log->debug("Given output_folder does not point to a file: {}", path_to_file.string());
    return false;
}

fs::path class_settings_reader::find_input_file(const fs::path &given_path) {

    //Check if file exists in the given path.
    fs::path complete_path = fs::system_complete(given_path);

    log->debug("Checking for input file: [ {} ] in path: [ {} ]", given_path.string() ,complete_path.string());
    if (check_if_input_file_exists(complete_path)){
        log->info("Found input file: [ {} ] in path: [ {} ]", given_path.string() ,fs::canonical(complete_path).string());
        found_file = true;
        return fs::canonical(complete_path);
    }

    //Check if file exists in the given path (if it is a relative path!), relative to the executable.
    if (given_path.is_relative()) {
        complete_path = fs::system_complete(fs::current_path() / given_path);
        log->debug("Checking for input file: [ {} ] in path: [ {} ]", given_path.string() ,complete_path.string());
        if (check_if_input_file_exists(complete_path)) {
            log->info("Found input file: [ {} ] in path: [ {} ]", given_path.string() ,fs::canonical(complete_path).string());
            found_file = true;
            return fs::canonical(complete_path);
        }
    }

    //Check if file exists in current directory
    complete_path = fs::system_complete(fs::current_path()/given_path.filename());
    log->debug("Checking for input file: [ {} ] in path: [ {} ]", given_path.string() ,complete_path.string());
    if(check_if_input_file_exists(complete_path)){
        log->info("Found input file: [ {} ] in path: [ {} ]", given_path.string() ,fs::canonical(complete_path).string());
        found_file = true;
        return fs::canonical(complete_path);
    }
    log->warn("Input file could not be found: [ {} ]", given_path.string());
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

