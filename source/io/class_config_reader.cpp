//
// Created by david on 2018-01-12.
//

#include "class_config_reader.h"
#include <fstream>
#include <h5pp/h5pp.h>
#include <sstream>
#include <tools/common/log.h>

class_config_reader::class_config_reader(const std::string &file_path_) : file_path(file_path_) {
    auto extension = file_path.extension();
    if(extension == ".config") {
        fs::path      config_file = find_config_file(file_path);
        std::ifstream file;
        if(not config_file.empty()) {
            try {
                file.open(config_file.c_str());
            } catch(std::exception &ex) { throw std::runtime_error(fmt::format("Could not open file [{}]: {} ", config_file.string(), ex.what())); }
        }
        // Load the file into memory
        if(file.is_open()) {
            file.clear();
            file.seekg(0, std::ios::beg);
            file_string = std::string((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));
        }
        file.close();
    } else if(extension == ".h5") {
        h5pp::File file(file_path.string(), h5pp::FilePermission::READONLY);
        file_string = file.readDataset<std::string>("common/config_file_contents");
    }

    // Extract the parameters and generate a key-value map
    std::istringstream file_stream(file_string);
    std::string        line;
    while(std::getline(file_stream, line)) {
        if(is_parameterline(line)) {
            line                  = remove_leading_spaces(line);
            auto comment_position = find_comment_symbols(line);
            if(comment_position == 0) continue;
            tools::log->trace("Reading line: {}", line);
            std::string        param_key;
            std::string        param_val;
            std::istringstream linestream(line);
            std::getline(linestream, param_key, '=');
            std::getline(linestream, param_val, '\n');
            param_val = param_val.substr(0, find_comment_symbols(param_val));
            param_key = remove_spaces(param_key);
            param_val = remove_spaces(param_val);
            std::transform(param_key.begin(), param_key.end(), param_key.begin(), ::tolower);
            param_map[param_key] = param_val;
        }
    }
}

std::string class_config_reader::get_config_file_as_string() { return file_string; }

std::string class_config_reader::get_config_filename() { return file_path.string(); }

bool class_config_reader::check_if_config_file_exists(const fs::path &path_to_file) {
    if(path_to_file.has_filename()) {
        if(fs::exists(path_to_file)) {
            std::ifstream in(path_to_file.c_str());
            if(in.is_open()) {
                in.close();
                return true;
            }
        }
        tools::log->debug("File does not exist: {}", path_to_file.string());
        return false;
    }
    tools::log->debug("Given output_folder does not point to a file: {}", path_to_file.string());
    return false;
}

fs::path class_config_reader::find_config_file(const fs::path &given_path) {
    // Check if file exists in the given path.
    fs::path complete_path = fs::absolute(given_path);
    if(fs::exists(given_path)) { complete_path = fs::canonical(given_path); }

    tools::log->debug("Checking for input file: [ {} ] in path: [ {} ]", given_path.string(), complete_path.string());
    if(check_if_config_file_exists(complete_path)) {
        tools::log->info("Found input file: [ {} ] in path: [ {} ]", given_path.string(), fs::canonical(complete_path).string());
        found_file = true;
        return fs::canonical(complete_path);
    }

    // Check if file exists in the given path (if it is a relative path!), relative to the executable.
    if(given_path.is_relative()) {
        complete_path = fs::absolute(fs::current_path() / given_path);
        tools::log->debug("Checking for input file: [ {} ] in path: [ {} ]", given_path.string(), complete_path.string());
        if(check_if_config_file_exists(complete_path)) {
            tools::log->info("Found input file: [ {} ] in path: [ {} ]", given_path.string(), fs::canonical(complete_path).string());
            found_file = true;
            return fs::canonical(complete_path);
        }
    }

    // Check if file exists in current directory
    complete_path = fs::absolute(fs::current_path() / given_path.filename());
    tools::log->debug("Checking for input file: [ {} ] in path: [ {} ]", given_path.string(), complete_path.string());
    if(check_if_config_file_exists(complete_path)) {
        tools::log->info("Found input file: [ {} ] in path: [ {} ]", given_path.string(), fs::canonical(complete_path).string());
        found_file = true;
        return fs::canonical(complete_path);
    }
    tools::log->warn("Input file could not be found: [ {} ]", given_path.string());
    found_file = false;
    return fs::path();
}

std::string class_config_reader::remove_spaces(std::string str) {
    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
    return str;
}

std::string class_config_reader::remove_leading_spaces(std::string str, const std::string &whitespace) {
    const auto strBegin = str.find_first_not_of(whitespace);
    if(strBegin == std::string::npos) return ""; // no content
    const auto strEnd   = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}

bool class_config_reader::has_only_digits(const std::string &str) { return str.find_first_not_of("+-0123456789") == std::string::npos; }

bool class_config_reader::is_parameterline(const std::string &str) { return str.find('=') != std::string::npos; }

std::string::size_type class_config_reader::find_comment_symbols(const std::string &str) {
    std::vector<std::string> comment_symbols = {"//", "/*", "#"};
    for(auto &sym : comment_symbols)
        if(str.find(sym) != std::string::npos) return str.find(sym);
    return std::string::npos;
}
