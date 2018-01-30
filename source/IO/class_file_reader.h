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
#include <IO/class_custom_cout.h>

namespace fs = std::experimental::filesystem::v1;


class class_file_reader {
private:
    fs::path    file_path;
    std::ifstream    file;
    bool check_if_input_file_exists(const fs::path &path_to_file);
    fs::path find_input_file(const fs::path &given_path);
    void remove_spaces(std::string &str);
    bool has_only_digits(const std::string s);
    bool is_parameterline(const std::string s);
    std::string::size_type find_comment_character(const std::string s);
    class_custom_cout ccout;
public:
    class_file_reader() = default;
    explicit class_file_reader(const fs::path &file_path_): file_path(file_path_) {
        try {
            file.open(find_input_file(file_path).c_str());
        }
        catch(std::exception &ex){
            std::cerr << "Exiting: " << ex.what() << std::endl;
            exit(1);
        }
    };

    template <typename T>
    T find_parameter(std::string param_requested, T default_value){
        if (file.is_open()){
            return find_parameter<T>(param_requested);
        }else{
            ccout(1) << "Missing input file: Using default value." << std::endl;
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
                    ccout(2) << "Loading line:     " << line << std::endl;
                    if constexpr (std::is_same<T,int>::value){
                        if (has_only_digits(param_val)) {
                            try {
                                T parsed_val = std::stoi(param_val);
                                return parsed_val;
                            }
                            catch (...) {std::cerr << "Error reading parameter from file: Unknown error." <<std::endl; }
                        }else{
                            std::cerr << "Error reading parameter from file. Wrong format: [" << param_val << "]. Expected an integer."<< std::endl;
                        }
                    }
                    if constexpr (std::is_same<T,long>::value){
                        if (has_only_digits(param_val)) {
                            try {
                                T parsed_val = std::stol(param_val);
                                return parsed_val;
                            }
                            catch (...) {std::cerr << "Error reading parameter from file: Unknown error." << std::endl; }
                        }else{
                            std::cerr << "Error reading parameter from file. Wrong format: [" << param_val << "]. Expected a long integer."<< std::endl;
                        }
                    }

                    if constexpr (std::is_same<T,double>::value){
                        try {
                            return std::stod(param_val);
                        }
                        catch (...) {
                            std::cerr << "Error reading parameter from file: Unknown error: [" << param_val << "]. Expected a double." << std::endl;
                        }
                    }
                    if constexpr (std::is_same<T,bool>::value) {
                        if (param_val == "true") { return true; }
                        if (param_val == "false") { return false;}
                    }

                    if constexpr (std::is_same<T,std::string>::value){
                        return param_val;
                    }

                    std::cerr << "\nCritical error when reading parameter from file. Possible type mismatch." << std::endl;
                    std::cerr << "Requested : [" << param_requested << "]" << " with type: [" << typeid(T).name() << "]"<< std::endl;
                    std::cerr << "Found key : [" << param_key << "]" << std::endl;
                    std::cerr << "Found val : [" << param_val << "]" << std::endl;
                    std::cerr << "Exiting..." << std::endl << std::flush;
                    exit(1);
                }
            }
            std::cerr << "Input file does not contain a parameter matching your query: ["  << param_requested << "]" << std::endl;
            std::cerr << "Exiting..." << std::endl  << std::flush;
            exit(1);
        }
        else {
            std::cerr << "Error: Input file has not been found yet and no default value was given. Exiting..." << std::endl;
            exit(1);
        }
    }
};



#endif //DMRG_CLASS_FILE_READER_H
