
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_io_class_settings_reader.h:

Program Listing for File class_settings_reader.h
================================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_io_class_settings_reader.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/io/class_settings_reader.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-01-12.
   //
   
   #pragma once
   
   
   #include<iostream>
   #include<fstream>
   #include<string>
   #include <iomanip>
   #include <cctype>
   
   #include <filesystem>
   #include <algorithm>
   #include <io/nmspc_logger.h>
   namespace fs = std::filesystem;
   
   
   class class_settings_reader {
   private:
       fs::path    file_path;
       std::string file_string = "";
       bool check_if_input_file_exists(const fs::path &path_to_file);
       fs::path find_input_file(const fs::path &given_path);
       void remove_spaces(std::string &str);
       bool has_only_digits(const std::string s);
       bool is_parameterline(const std::string s);
       std::unordered_map<std::string,std::string> param_map;
       std::string::size_type find_comment_character(const std::string s);
       std::shared_ptr<spdlog::logger> log;
   public:
       bool found_file = false;
       class_settings_reader() = default;
       explicit class_settings_reader(const fs::path &file_path_, std::string logName="DMRG");
   
       std::string get_input_file_as_string();
       std::string get_input_filename();
   
   
   
       template <typename T>
       void find_parameter(std::string param_name, T &param_value){
           try{
               T new_value = find_parameter<T>(param_name);
               param_value = new_value;
               log->trace("Successfully parsed parameter: {:<40} = {:<20}", param_name, param_value);
           }catch (std::exception & ex){
               log->error("Failed to read parameter [{}]: {}",param_name, ex.what() );
               log->error("Using default parameter: [{}] = {}",param_name,param_value );
   
           }
       }
   
   
       template <typename T>
       T parse_param(const std::string &param_val){
           try{
               if constexpr (std::is_same<T,int>::value)    return (T) std::stoi(param_val);
               if constexpr (std::is_same<T,long>::value)   return (T) std::stol(param_val);
               if constexpr (std::is_same<T,size_t>::value) return (T) std::stol(param_val);
               if constexpr (std::is_same<T,double>::value) return (T) std::stod(param_val);
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
   
   
   
   //        if (file.is_open()) {
   //            file.clear();
   //            file.seekg(0, std::ios::beg);
   //
   //
   //            std::string param_key;
   //            std::string param_val;
   //            std::string line;
   //
   //            while (!file.eof()) {
   //                getline(file, line);
   //                if(!is_parameterline(line)){continue;}
   //                std::istringstream is(line);
   //                getline(is,param_key, '=');
   //                getline(is,param_val,  '\n');
   //                param_val = param_val.substr(0, find_comment_character(param_val));
   //                remove_spaces(param_requested);
   //                remove_spaces(param_key);
   //                remove_spaces(param_val);
   //
   //                std::transform(param_key.begin(), param_key.end(), param_key.begin(), ::tolower);
   //                std::transform(param_requested.begin(), param_requested.end(), param_requested.begin(), ::tolower);
   //                if (param_requested == param_key && !param_key.empty()) {
   //                    log->debug("Loading line: {}",line);
   //                    try {
   //                        return parse_param<T>(param_val);
   //                    }catch (std::exception &ex){
   //                        throw std::runtime_error(fmt::format("Error parsing parameter. Requested [{}]. Found key [{}] with value [{}]. Reason {}", param_requested,param_key,param_val, ex.what()));
   //                    }
   //                }
   //            }
   //            throw std::runtime_error(fmt::format("Input file does not contain a parameter matching your query: [{}]", param_requested));
   //        }
   //        else {
   //            throw std::runtime_error(fmt::format("Error: Input file could not be opened"));
   //        }
       }
   };
   
