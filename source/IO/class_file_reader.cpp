//
// Created by david on 2018-01-12.
//

#include "class_file_reader.h"
//#include <directory.h>
bool class_file_reader::check_if_input_file_exists(const fs::path &path_to_file){
    if (path_to_file.has_filename()){
        if(fs::exists(path_to_file)){
            std::ifstream in(path_to_file.c_str());
            if(in.is_open()){
                in.close();
                ccout(1) << "Found input file: " << fs::canonical(path_to_file) << '\n';
                return true;
            }
        }
        ccout(1) << "File does not exist: " << path_to_file << std::endl;
        return false;
    }
    ccout(1) << "Given output_folder does not point to a file: "  << path_to_file << std::endl;
    return false;
}

fs::path class_file_reader::find_input_file(const fs::path &given_path) {

    //Check if file exists in the given path.
    fs::path complete_path = fs::system_complete(given_path);
    std::cout << "Checking for input file: ["<< given_path << "] in path: " << complete_path << std::endl;
    if (check_if_input_file_exists(complete_path)){
        return fs::canonical(complete_path);
    }

    //Check if file exists in the given path (if it is a relative path!), relative to the executable.
    if (given_path.is_relative()) {
        complete_path = fs::system_complete(fs::current_path() / given_path);
        std::cout << "Checking for input file: ["<< given_path << "] in path: " << complete_path << std::endl;
        if (check_if_input_file_exists(complete_path)) {
            return fs::canonical(complete_path);
        }
    }

    //Check if file exists in current directory
    complete_path = fs::system_complete(fs::current_path()/given_path.filename());
    std::cout << "Checking for input file: ["<< given_path << "] in path: " << complete_path << std::endl;
    if(check_if_input_file_exists(complete_path)){
        return fs::canonical(complete_path);
    }



    //Search recursively
    std::vector<fs::path> matching_files;
    fs::path recurse_from_path ;// fs::current_path() / given_path.has_parent_path();
    if(given_path.is_relative() and given_path.has_parent_path()){
            recurse_from_path = fs::system_complete(fs::current_path()/given_path.parent_path());
    }else if(given_path.has_parent_path()) {
        recurse_from_path = given_path.parent_path();
    }else{
        recurse_from_path = fs::current_path();
    }
    ccout(1) << "Searching recursively for file [" << given_path.filename() << "] in folder: " << recurse_from_path << std::endl;
    for(auto& p: fs::recursive_directory_iterator(recurse_from_path)) {
//        std::cout << "Trying path: " << p.path() << std::endl;
        if (p.path().filename() == given_path.filename()  ) {
            if(check_if_input_file_exists(p)) {
                return fs::canonical(p.path());
            }
        }
    }

//    //Search recursively inside current_path/
//    recurse_from_path = fs::current_path();
//    ccout(1) << "Searching recursively for file [" << given_path.filename() << "] in folder EXECUTABLE_PATH/..." << std::endl;
//    for(auto& p: fs::recursive_directory_iterator(recurse_from_path)) {
////        std::cout << "Trying path: " << p.path() << std::endl;
//        if (p.path().filename() == given_path.filename() ) {
//            if(check_if_input_file_exists(p)) {
//                return fs::canonical(p.path());
//            }
//        }
//    }

    std::cerr << "Input file could not be found. Exiting"  << std::flush << std::endl;
    exit(1);
//    if(matching_files.size() > 1){
//        std::cout << std::flush;
//        std::cerr << "Found multiple files with the given output_filename: [" << given_filename << "]." << std::endl;
//        std::cerr << "Exiting..." << std::endl;
//        exit(1);
//    }else{
//        std::cout << std::flush;
//        return fs::canonical(matching_files[0]);
//    }

}

void class_file_reader::remove_spaces(std::string &str){
    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
}

bool class_file_reader::has_only_digits(const std::string s){
    return s.find_first_not_of( "+-0123456789" ) == std::string::npos;
}

bool class_file_reader::is_parameterline(const std::string s){
    return s.find("=") != std::string::npos;
}

std::string::size_type class_file_reader::find_comment_character(const std::string s){
    std::vector<std::string> comment_symbols = {"//", "/*", "#"};
    for(auto &sym : comment_symbols){
        if(s.find(sym) != std::string::npos){
            return s.find(sym);
        }
    }
    return s.npos;
}

