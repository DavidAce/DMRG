#pragma once
#include <io/nmspc_filesystem.h>
#include <simulation/enums.h>
#include <string>
#include <tools/common/log.h>
#include <unordered_map>

class class_dmrg_config {
    public:
    bool               file_exists;
    fs::path           file_path;
    fs::file_time_type file_date;

    explicit class_dmrg_config(const fs::path &cfg_or_h5_file);


    void                      load();
    void                      unload();
    [[nodiscard]] std::string get_config_file_as_string();

    template<typename T>
    void load_parameter(const std::string &param_name, T &param_value) {
        try {
            T new_value = find_parameter<T>(param_name);
            param_value = new_value;
            if constexpr(std::is_enum_v<T>)
                tools::log->debug("Loaded parameter: {:<48} = {:<20}", param_name, enum2str<T>(param_value));
            else
                tools::log->debug("Loaded parameter: {:<48} = {:<20}", param_name, param_value);

        } catch(std::exception &ex) {
            tools::log->info("Failed to read parameter [{}]: {}", param_name, ex.what());
        }
    }

    private:
    std::unordered_map<std::string, std::string> param_map;
    [[nodiscard]] static std::string             remove_spaces(std::string str);
    [[nodiscard]] static std::string             remove_leading_spaces(std::string str, const std::string &whitespace = " \t");
    [[nodiscard]] static bool                    is_parameterline(const std::string &s);
    [[nodiscard]] static std::string::size_type  find_comment_character(const std::string &s);

    template<typename T>
    [[nodiscard]] T parse_param(const std::string &param_val, const std::string &param_name) {
        try {
            if(param_val.empty()) throw std::range_error("Parameter [" + param_name + "] has no value");
            if constexpr(std::is_same_v<T, unsigned int>) {
                int val = (T) std::stoi(param_val);
                if(val < 0) throw std::runtime_error("Read negative value for unsigned parameter: " + std::to_string(val));
                return (T) val;
            }
            if constexpr(std::is_same_v<T, unsigned long>) return (T) std::stoul(param_val);
            if constexpr(std::is_same_v<T, unsigned long long>) return (T) std::stoull(param_val);
            if constexpr(std::is_same_v<T, int>) return (T) std::stoi(param_val);
            if constexpr(std::is_same_v<T, long>) return (T) std::stol(param_val);
            if constexpr(std::is_same_v<T, long long>) return (T) std::stoll(param_val);
            if constexpr(std::is_same_v<T, size_t>) return (T) std::stol(param_val);
            if constexpr(std::is_same_v<T, double>) return (T) std::stod(param_val);
            if constexpr(std::is_enum_v<T>) return str2enum<T>(param_val);
            if constexpr(std::is_same<T, std::string>::value) return param_val;
            if constexpr(std::is_same<T, bool>::value) {
                if(param_val == "true") return true;
                if(param_val == "false") return false;
                throw std::runtime_error(fmt::format("Expected true or false, got {}", param_val));
            }
            throw std::runtime_error("Type mismatch on parameter: " + param_val);
        } catch(std::exception &ex) {
            throw std::runtime_error("Error parsing param: " + std::string(ex.what()));
        } catch(...) {
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
        } catch(std::range_error &ex) {
            throw std::runtime_error(ex.what());
        } catch(std::exception &ex) {
            throw std::runtime_error(fmt::format("Error parsing the requested parameter: [{}]: {}", param_requested, ex.what()));
        }
    }
};