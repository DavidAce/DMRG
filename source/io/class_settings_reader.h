//
// Created by david on 2018-01-12.
//

#ifndef DMRG_CLASS_FILE_READER_H
#define DMRG_CLASS_FILE_READER_H



#include<iostream>
#include<fstream>
#include<string>
#include <iomanip>
#include <cctype>

#include <experimental/filesystem>
#include <algorithm>
#include <io/nmspc_logger.h>
namespace fs = std::experimental::filesystem;


class class_settings_reader {
private:
    fs::path    file_path;
    std::ifstream    file;
    bool check_if_input_file_exists(const fs::path &path_to_file);
    fs::path find_input_file(const fs::path &given_path);
    void remove_spaces(std::string &str);
    bool has_only_digits(const std::string s);
    bool is_parameterline(const std::string s);
    std::string::size_type find_comment_character(const std::string s);
    std::shared_ptr<spdlog::logger> log;
public:
    bool found_file = false;
    class_settings_reader() = default;
    explicit class_settings_reader(const fs::path &file_path_, std::string logName="DMRG");

    std::string get_input_file();
    std::string get_input_filename();

    template <typename T>
    T find_parameter(std::string param_requested, T default_value){
        if (file.is_open()){
            return find_parameter<T>(param_requested);
        }else{
            log->warn("Missing input file: Using default value: {}", default_value);
            return default_value;
        }
    }

    template<typename T>
    T find_parameter(std::string param_requested){
        if (file.is_open()) {
            file.clear();
            file.seekg(0, std::ios::beg);


            std::string param_key;
            std::string param_val;
            std::string line;

            while (!file.eof()) {
                getline(file, line);
                if(!is_parameterline(line)){continue;}
                std::istringstream is(line);
                getline(is,param_key, '=');
                getline(is,param_val,  '\n');
                param_val = param_val.substr(0, find_comment_character(param_val));
                remove_spaces(param_requested);
                remove_spaces(param_key);
                remove_spaces(param_val);

                std::transform(param_key.begin(), param_key.end(), param_key.begin(), ::tolower);
                std::transform(param_requested.begin(), param_requested.end(), param_requested.begin(), ::tolower);
                if (param_requested == param_key && !param_key.empty()) {
                    log->debug("Loading line: {}",line);
                    if constexpr (std::is_same<T,int>::value){
                        if (has_only_digits(param_val)) {
                            try {
                                T parsed_val = std::stoi(param_val);
                                return parsed_val;
                            }
                            catch (...) {log->error("Error reading parameter from file: Unknown error."); }
                        }else{
                            log->error("Error reading parameter from file. Wrong format: [{}]. Expected an integer", param_val);
                        }
                    }
                    if constexpr (std::is_same<T,long>::value){
                        if (has_only_digits(param_val)) {
                            try {
                                T parsed_val = std::stol(param_val);
                                return parsed_val;
                            }
                            catch (...) {log->error("Error reading parameter from file: Unknown error."); }

                        }else{
                            log->error("Error reading parameter from file. Wrong format: [{}]. Expected a long integer", param_val);

                        }
                    }

                    if constexpr (std::is_same<T,double>::value){
                        try {
                            return std::stod(param_val);
                        }
                        catch (...) {
                            log->error("Error reading parameter from file. Wrong format: [{}]. Expected double", param_val);
                        }
                    }
                    if constexpr (std::is_same<T,bool>::value) {
                        if (param_val == "true") { return true; }
                        if (param_val == "false") { return false;}
                    }

                    if constexpr (std::is_same<T,std::string>::value){
                        return param_val;
                    }
                    log->critical("Critical error when reading parameter from file. Possible type mismatch.");
                    log->critical("Requested : [{}] with type [{}]", param_requested ,typeid(T).name());
                    log->critical("Found key : [{}]", param_key);
                    log->critical("Found val : [{}]", param_val);
                    log->critical("Exiting...");
                    exit(1);
                }
            }

            log->critical("Input file does not contain a parameter matching your query: [{}]", param_requested);
            log->critical("Exiting...");
            exit(1);
        }
        else {
            log->critical("Error: Input file has not been found yet and no default value was given. Exiting...");
            exit(1);
        }
    }
};



#endif //DMRG_CLASS_FILE_READER_H
