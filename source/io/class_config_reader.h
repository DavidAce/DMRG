//
// Created by david on 2018-01-12.
//

#pragma once

#include <algorithm>
#include <cctype>
#include <fstream>
#include <io/nmspc_filesystem.h>
#include <io/nmspc_logger.h>
#include <iomanip>
#include <iostream>
#include <simulation/enums.h>
#include <string>

class class_config_reader {
private:
    tools::fs::path    file_path;
    std::string file_string = "";
    bool check_if_config_file_exists(const tools::fs::path &path_to_file);
    tools::fs::path find_config_file(const tools::fs::path &given_path);
    void remove_spaces(std::string &str);
    bool has_only_digits(const std::string s);
    bool is_parameterline(const std::string s);
    std::unordered_map<std::string,std::string> param_map;
    std::string::size_type find_comment_character(const std::string s);
    std::shared_ptr<spdlog::logger> log;
public:
    bool found_file = false;
    class_config_reader() = default;
    explicit class_config_reader(const std::string &file_path_, std::string logName="DMRG");

    std::string get_config_file_as_string();
    std::string get_config_filename();



    template <typename T>
    void find_parameter(std::string param_name, T &param_value){
        try{
            T new_value = find_parameter<T>(param_name);
            param_value = new_value;
            if constexpr (std::is_enum_v<T>)
                log->trace("Successfully parsed parameter: {:<40} = {:<20}", param_name, enum2str<T>(param_value));
            else
                log->trace("Successfully parsed parameter: {:<40} = {:<20}", param_name, param_value);

        }catch (std::exception & ex){
            log->error("Failed to read parameter [{}]: {}",param_name, ex.what() );
            log->error("Using default parameter: [{}] = {}",param_name,param_value );

        }
    }


    template <typename T>
    T parse_param(const std::string &param_val){
        try{
            if constexpr (std::is_same_v<T,int>)    return (T) std::stoi(param_val);
            if constexpr (std::is_same_v<T,long>)   return (T) std::stol(param_val);
            if constexpr (std::is_same_v<T,size_t>) return (T) std::stol(param_val);
            if constexpr (std::is_same_v<T,double>) return (T) std::stod(param_val);
            if constexpr (std::is_enum_v<T>)        return str2enum<T>(param_val);
            if constexpr (std::is_same<T,std::string>::value) return param_val;
            if constexpr (std::is_same<T,bool>::value){
                if (param_val == "true") return true;
                if (param_val == "false")return false;
                throw std::runtime_error(fmt::format("Expected true or false, got {}", param_val));
            }
            throw std::runtime_error("Type mismatch on parameter: " + param_val);
        }
        catch (std::exception &ex){
            throw std::runtime_error("Error parsing param: " + std::string(ex.what()));
        }
        catch (...){
            throw std::runtime_error("Error parsing param: Unknown error");
        }
    }


    template<typename T>
    T find_parameter(std::string param_requested){

        remove_spaces(param_requested);
        std::transform(param_requested.begin(), param_requested.end(), param_requested.begin(), ::tolower);
        try{
            std::string param_val = param_map[param_requested];
            return parse_param<T>(param_val);
        }
        catch(std::exception &ex){
            throw std::runtime_error(fmt::format("Error parsing parameter. Requested [{}]. Reason {}", param_requested, ex.what()));
        }


    }
};

