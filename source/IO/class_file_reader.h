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
#include <directory.h>


namespace fs = std::experimental::filesystem::v1;


class class_file_reader {
private:
    fs::path    file_path;
    std::ifstream    file;

    bool check_if_input_file_exists(const fs::path &path_to_file){
        if (path_to_file.has_filename()){
            if(fs::exists(path_to_file)){
                std::ifstream in(path_to_file.c_str());
                if(in.is_open()){
                    in.close();
                    std::cout << "Found input file: " << path_to_file << '\n';
                    return true;
                }
            }
            std::cout << "File does not exist: " << path_to_file << std::endl;
            return false;
        }
        std::cout << "Given path does not point to a file: "  << path_to_file << std::endl;
        return false;
    }

    fs::path find_input_file(const fs::path &given_path) {

        //Check if file exists in path as given, relative to the executable.
        fs::path complete_path = fs::system_complete(given_path);
        if (check_if_input_file_exists(complete_path)){
            return fs::canonical(complete_path);
        }

        //Check if file exists in path as given, relative to the project root folder.
        complete_path = fs::system_complete(fs::path(directory::SOURCE_DIR) / given_path);
        if (check_if_input_file_exists(complete_path)){
            return fs::canonical(complete_path);
        }

        //As a last resort, search in input/ folder
        std::vector<fs::path> matching_files;
        fs::path input_abs_path(directory::INPUT_DIR);
        fs::path given_filename = given_path.filename();
        std::cout << "Searching for file " << given_filename << " in folder PROJECT_ROOT/input" << std::endl;

        for(auto& p: fs::recursive_directory_iterator(input_abs_path)) {
            if (p.path().filename() == given_filename ) {
                if(check_if_input_file_exists(p)) {
                    matching_files.emplace_back(p.path());
                }
            }
        }

        if(matching_files.size() > 1){
            std::cout << std::flush;
            std::cerr << "Found multiple files with the given filename: [" << given_filename << "]." << std::endl;
            std::cerr << "Exiting..." << std::endl;
            exit(1);
        }else{
            std::cout << std::flush;
            return fs::canonical(matching_files[0]);
        }
    }


    void remove_spaces(std::string &str){
        str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
    }
    bool has_only_digits(const std::string s){
        return s.find_first_not_of( "+-0123456789" ) == std::string::npos;
    }

    bool is_parameterline(const std::string s){
        return s.find("=") != std::string::npos;
    }

    auto find_comment_character(const std::string s){
        return s.find_first_of( "/#!*{}()&$@;" );
    }

public:
    class_file_reader(){}

    explicit class_file_reader(const fs::path &file_path_): file_path(file_path_) {
        try {
            file.open(find_input_file(file_path).c_str());
        }
        catch(std::exception &ex){
            std::cerr << "Exiting: " << ex.what() << std::endl;
            exit(1);
        }
    };

    void set_inputfile(const fs::path &file_path_){
        file_path = file_path_;
        try {
            file.open(find_input_file(file_path).c_str());
        }
        catch(std::exception &ex){
            std::cerr << "Exiting: " << ex.what() << std::endl;
            exit(1);
        }    }

    template <typename T>
    T find_parameter(std::string param_requested, T default_value){
        if (file.is_open()){
            return find_parameter<T>(param_requested);
        }else{
            std::cout << "Missing input file: Using default value." << std::endl;
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
                    std::cout << "Loading line:     " << line << std::endl;
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

                    std::cerr << "Critical error when reading parameter from file. Possible type mismatch." << std::endl;
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
