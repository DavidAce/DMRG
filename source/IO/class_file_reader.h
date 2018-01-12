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


using namespace std;
namespace fs = std::experimental::filesystem::v1;

class class_file_reader {
private:
    fs::path    file_path;
    ifstream    file;

    fs::path find_file(const fs::path &given_path) {
        if (given_path.has_filename()){
            if(fs::exists(given_path)){
                ifstream in(given_path.c_str());
                if(in.is_open()){
                    in.close();
                    std::cout << "Found input file: " << given_path << '\n';
                    return given_path;
                }
            }
            cout << "Given path is valid but points to the wrong folder." << endl;
            if(given_path.is_absolute()){
                throw(std::runtime_error("Could not find input file"));
            }else{
                cout << "Searching for file in parent folders..." << endl << flush;
            }
        }


        fs::path current_abs_folder     = fs::current_path();
        fs::path current_abs_file_path  = current_abs_folder / given_path;


        //Try to find the file
        while (true) {
            std::cout << "Trying path: " << current_abs_file_path << '\n';
            if (fs::exists(current_abs_file_path)) {
                current_abs_file_path = fs::canonical(current_abs_file_path);
                std::cout << "Found file: " << current_abs_file_path << '\n';
                break;
            } else if (current_abs_folder.has_parent_path()) {
                current_abs_folder = current_abs_folder.parent_path();
                current_abs_file_path = current_abs_folder / given_path;
            } else {
                throw(std::runtime_error(string("Given path [") + given_path.c_str() + string("] cannot be found.\n Try relative or absolute paths to an existing file.")));
            }
        }
        cout << flush;
        return current_abs_file_path;
    }
    void remove_spaces(std::string &str){
        str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
    }
    bool has_only_digits(const string s){
        return s.find_first_not_of( "+-0123456789" ) == string::npos;
    }

    bool is_parameterline(const string s){
        return s.find("=") != string::npos;
    }

    auto find_comment_character(const string s){
        return s.find_first_of( "/#!*{}()&$@;" );
    }

public:
    explicit class_file_reader(const fs::path &file_path_): file_path(file_path_) {
        try {
            file.open(find_file(file_path).c_str());
        }
        catch(std::exception &ex){
            cerr << "Exiting: " << ex.what() << endl;
            exit(1);
        }
    };

    template<typename T>
    T find_parameter(std::string param_requested){
        if (file.is_open()) {
            file.clear();
            file.seekg(0, ios::beg);


            string param_key;
            string param_val;
            string line;

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
                    std::cout << "Reading Line:     " << line << endl;
                    if constexpr (std::is_same<T,int>::value){
                        if (has_only_digits(param_val)) {
                            try {
                                T parsed_val = std::stoi(param_val);
                                return parsed_val;
                            }
                            catch (...) {cerr << "Error reading parameter from file: Unknown error." << endl; }
                        }else{
                            cerr << "Error reading parameter from file. Wrong format: [" << param_val << "]. Expected an integer."<< endl;
                        }
                    }
                    if constexpr (std::is_same<T,double>::value){
                        try {
                            return std::stod(param_val);
                        }
                        catch (...) {
                            cerr << "Error reading parameter from file: Unknown error: [" << param_val << "]. Expected a double." << endl;
                        }
                    }
                    if constexpr (std::is_same<T,bool>::value) {
                        if (param_val == "true") { return true; }
                        if (param_val == "false") { return false;}
                    }

                    if constexpr (std::is_same<T,std::string>::value){
                        return param_val;
                    }

                    cerr << "Critical error when reading parameter from file. Possible type mismatch." << endl;
                    cerr << "Requested : [" << param_requested << "]" << " with type: [" << typeid(T).name() << "]"<< endl;
                    cerr << "Found key : [" << param_key << "]" << endl;
                    cerr << "Found val : [" << param_val << "]" << endl;
                    cerr << "Exiting..." << endl << flush;
                    exit(1);

                }
            }
            std::cerr << "Input file does not contain a parameter matching your query: ["  << param_requested << "]" << endl;
            std::cerr << "Exiting..." << endl  << flush;
            exit(1);
        }
    }
};



#endif //DMRG_CLASS_FILE_READER_H
