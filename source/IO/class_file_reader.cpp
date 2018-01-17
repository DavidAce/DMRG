//
// Created by david on 2018-01-12.
//

#include "class_file_reader.h"
#include <directory.h>
bool class_file_reader::check_if_input_file_exists(const fs::path &path_to_file){
    if (path_to_file.has_filename()){
        if(fs::exists(path_to_file)){
            std::ifstream in(path_to_file.c_str());
            if(in.is_open()){
                in.close();
                ccout(1) << "Found input file: " << path_to_file << '\n';
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
    if (check_if_input_file_exists(complete_path)){
        return fs::canonical(complete_path);
    }

    //Check if file exists in the given path (if it is relative path!), relative to the project root folder.
    if (given_path.is_relative()) {
        complete_path = fs::system_complete(fs::path(directory::PROJECT_DIR) / given_path);
        if (check_if_input_file_exists(complete_path)) {
            return fs::canonical(complete_path);
        }
    }

    //As a last resort, search in input/ folder
    std::vector<fs::path> matching_files;
    fs::path input_abs_path(directory::PROJECT_DIR + "/input");
    fs::path given_filename = given_path.filename();

    ccout(1) << "Searching for file " << given_filename << " in folder PROJECT_ROOT/input" << std::endl;

    for(auto& p: fs::recursive_directory_iterator(input_abs_path)) {
        if (p.path().filename() == given_filename ) {
            if(check_if_input_file_exists(p)) {
                matching_files.emplace_back(p.path());
            }
        }
    }

    if(matching_files.size() > 1){
        std::cout << std::flush;
        std::cerr << "Found multiple files with the given output_filename: [" << given_filename << "]." << std::endl;
        std::cerr << "Exiting..." << std::endl;
        exit(1);
    }else{
        std::cout << std::flush;
        return fs::canonical(matching_files[0]);
    }
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
    return s.find_first_of( "/#!*{}()&$@;" );
}

