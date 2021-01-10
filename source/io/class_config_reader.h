//
// Created by david on 2018-01-12.
//

#pragma once

#include <algorithm>
#include <config/enums.h>
#include <io/nmspc_filesystem.h>
#include <string>
#include <tools/common/log.h>

class class_config_reader {
    private:
    fs::path                                     file_path;
    std::string                                  file_string;
    std::unordered_map<std::string, std::string> param_map;

    [[nodiscard]] bool                          check_if_config_file_exists(const fs::path &path_to_file);
    [[nodiscard]] fs::path                      find_config_file(const fs::path &given_path);
    [[nodiscard]] std::string                   remove_spaces(std::string str);
    [[nodiscard]] std::string                   remove_leading_spaces(std::string str, const std::string &whitespace = " \t");
    [[nodiscard]] bool                          has_only_digits(const std::string &s);
    [[nodiscard]] static bool                   is_parameterline(const std::string &s);
    [[nodiscard]] static std::string::size_type find_comment_symbols(const std::string &str);

    public:
    bool found_file       = false;
    class_config_reader() = default;
    explicit class_config_reader(const std::string &file_path_);

    [[nodiscard]] std::string get_config_file_as_string();
    [[nodiscard]] std::string get_config_filename();

    template<typename T>
    void find_parameter(T &param_value, const std::string &param_name) {
        try {
            T new_value = find_parameter<T>(param_name);
            param_value = new_value;
            if constexpr(std::is_enum_v<T>)
                tools::log->debug("Loaded parameter: {:<48} = {:<20}", param_name, enum2str<T>(param_value));
            else
                tools::log->debug("Loaded parameter: {:<48} = {:<20}", param_name, param_value);

        } catch(std::exception &ex) {
            tools::log->info("Failed to read parameter [{}]: {}", param_name, ex.what());
            tools::log->info("Using default parameter: [{}] = [{}]", param_name, param_value);
        }
    }

    private:
    template<typename T>
    [[nodiscard]] T parse_param(const std::string &param_val, const std::string &param_name) {
        try {
            if(param_val.empty()) throw std::range_error("Parameter [" + param_name + "] has no value");
            if constexpr(std::is_same_v<T, int>) return static_cast<T>(std::stoi(param_val));
            if constexpr(std::is_same_v<T, long>) return static_cast<T>(std::stol(param_val));
            if constexpr(std::is_same_v<T, size_t>) return static_cast<T>(std::stol(param_val));
            if constexpr(std::is_same_v<T, double>) return static_cast<T>(std::stod(param_val));
            if constexpr(std::is_enum_v<T>) return str2enum<T>(param_val);
            if constexpr(std::is_same<T, std::string>::value) return param_val;
            if constexpr(std::is_same<T, bool>::value) {
                if(param_val == "true") return true;
                if(param_val == "false") return false;
                throw std::runtime_error(fmt::format("Expected true or false, got {}", param_val));
            }
            throw std::runtime_error("Type mismatch on parameter: " + param_val);
        } catch(std::exception &ex) { throw std::runtime_error("Error parsing param: " + std::string(ex.what())); } catch(...) {
            throw std::runtime_error("Error parsing param: Unknown error");
        }
    }

    template<typename T>
    [[nodiscard]] T find_parameter(std::string param_requested) {
        param_requested = remove_spaces(param_requested);
        std::transform(param_requested.begin(), param_requested.end(), param_requested.begin(), ::tolower);
        try {
            if(param_map.find(param_requested) != param_map.end()) {
                return parse_param<T>(param_map[param_requested], param_requested);
            } else
                throw std::range_error(fmt::format("Config file does not specify the requested parameter: [{}]", param_requested));
        } catch(std::range_error &ex) { throw std::runtime_error(ex.what()); } catch(std::exception &ex) {
            throw std::runtime_error(fmt::format("Error parsing the requested parameter: [{}]: {}", param_requested, ex.what()));
        }
    }
};
