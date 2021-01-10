//
// Created by david on 2020-04-24.
//

#include "class_dmrg_config.h"
#include <fstream>
#include <h5pp/h5pp.h>
class_dmrg_config::class_dmrg_config(const fs::path &cfg_or_h5_file) {
    file_exists = fs::exists(cfg_or_h5_file);
    file_path   = fs::absolute(cfg_or_h5_file);
    if(not file_exists) return;
    file_path = fs::canonical(cfg_or_h5_file);
    file_date = fs::last_write_time(cfg_or_h5_file);
}

void class_dmrg_config::load() {
    // Extract the parameters and generate a key-value map
    if(not file_exists) return tools::log->error("Can't load config file [{}]: it does not exist", file_path.string());
    if(not param_map.empty()) return;

    std::istringstream file_stream(get_config_file_as_string());
    std::string        line;
    while(std::getline(file_stream, line)) {
        if(is_parameterline(line)) {
            line                  = remove_leading_spaces(line);
            auto comment_position = find_comment_character(line);
            if(comment_position == 0) continue;
            //            tools::log->trace("Reading line: {}", line);
            std::string        param_key;
            std::string        param_val;
            std::istringstream linestream(line);
            std::getline(linestream, param_key, '=');
            std::getline(linestream, param_val, '\n');
            param_val = param_val.substr(0, find_comment_character(param_val));
            param_key = remove_spaces(param_key);
            param_val = remove_spaces(param_val);
            std::transform(param_key.begin(), param_key.end(), param_key.begin(), ::tolower);
            param_map[param_key] = param_val;
        }
    }
}

void class_dmrg_config::unload() { param_map.clear(); }

std::string class_dmrg_config::get_config_file_as_string() {
    // Load the file into memory
    if(file_path.extension() == ".cfg") {
        std::ifstream file;
        try {
            file.open(file_path);
        } catch(std::exception &ex) { throw std::runtime_error(fmt::format("Could not open file [ {} ]: {}", file_path.string(), ex.what())); }
        if(file.is_open()) {
            file.clear();
            file.seekg(0, std::ios::beg);
            return std::string((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));
        }
        file.close();
    } else if(file_path.extension() == ".h5") {
        h5pp::File file(file_path.string(), h5pp::FilePermission::READONLY);
        return file.readDataset<std::string>("common/config_file_contents");
    } else {
        throw std::runtime_error(
            fmt::format("Tried to read config from file with unknown extension [{}], expected .h5 or .config", file_path.filename().string()));
    }
    return std::string();
}

std::string class_dmrg_config::remove_spaces(std::string str) {
    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
    return str;
}

std::string class_dmrg_config::remove_leading_spaces(std::string str, const std::string &whitespace) {
    const auto strBegin = str.find_first_not_of(whitespace);
    if(strBegin == std::string::npos) return ""; // no content
    const auto strEnd   = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}

bool class_dmrg_config::is_parameterline(const std::string &str) { return str.find('=') != std::string::npos; }

std::string::size_type class_dmrg_config::find_comment_character(const std::string &str) {
    std::vector<std::string> comment_symbols = {"//", "/*", "#"};
    for(auto &sym : comment_symbols)
        if(str.find(sym) != std::string::npos) return str.find(sym);
    return str.npos;
}
