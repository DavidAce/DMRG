#pragma once
#include <h5pp/details/h5ppFormat.h>
#include <stdexcept>
#include <string>

namespace tools::parse {

    template<typename T>
    [[nodiscard]] T parse_param(std::string_view param_val, const std::string &param_name) {
        try {
            if(param_val.empty()) throw std::range_error("Parameter [" + param_name + "] has no value");
            if constexpr(std::is_same_v<T, int>) return (T) std::stoi(param_val.data());
            if constexpr(std::is_same_v<T, long>) return (T) std::stol(param_val.data());
            if constexpr(std::is_same_v<T, size_t>) return (T) std::stol(param_val.data());
            if constexpr(std::is_same_v<T, double>) return (T) std::stod(param_val.data());
            if constexpr(std::is_same<T, std::string>::value) return param_val;
            if constexpr(std::is_same<T, bool>::value) {
                if(param_val == "true") return true;
                if(param_val == "false") return false;
                throw std::range_error(h5pp::format("Expected true or false, got {}", param_val));
            }
            throw std::runtime_error(h5pp::format("Type mismatch on parameter: {}", param_val));
        } catch(std::exception &ex) { throw std::runtime_error(h5pp::format("Error parsing param: {}", ex.what())); } catch(...) {
            throw std::runtime_error("Error parsing param: Unknown error");
        }
    }

    template<typename T>
    [[nodiscard]] T extract_digits_from_h5_filename(const std::string &input) {
        std::string stem = h5pp::fs::path(input).stem().string();
        std::string seed_str;
        for(const auto &c : stem)
            if(std::isdigit(c)) seed_str.push_back(c);
        try {
            return parse_param<T>(seed_str, input);
        } catch(const std::exception &err) { throw std::runtime_error(h5pp::format("Could not convert {} to a number: {}", input, err.what())); }
    }

    template<typename T>
    [[nodiscard]] T extract_parameter_from_path(const std::string &input, const std::string &param_name) {
        std::string param_str = input.substr(input.find(param_name) + param_name.size(), input.find('/') - 1);
        try {
            return parse_param<T>(param_str, input);
        } catch(const std::exception &err) { throw std::runtime_error(h5pp::format("Could not convert {} to a number: {}", input, err.what())); }
    }
}
