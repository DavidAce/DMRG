#pragma once

#include <optional>
#include <string>
#include <string_view>

namespace text {
    extern bool        endsWith(std::string_view str, std::string_view suffix);
    extern bool        startsWith(std::string_view str, std::string_view prefix);
    extern std::string replace(std::string_view str, std::string_view from, std::string_view to);

    template<typename T>
    [[nodiscard]] std::optional<T> extract_value(std::string_view input) {
        if constexpr(std::is_same_v<T, bool>) {
            if(input.find("true") != std::string_view::npos) return true;
            if(input.find("True") != std::string_view::npos) return true;
            if(input.find("false") != std::string_view::npos) return false;
            if(input.find("False") != std::string_view::npos) return false;
        } else if constexpr(std::is_arithmetic_v<T> or std::is_same_v<T, std::string>) {
            static constexpr char const *digits = "0123456789";
            static constexpr char const *figits = "0123456789.";
            auto                         bgn    = input.find_first_of(digits);
            if(bgn != std::string::npos) {
                auto end = input.find_first_not_of(figits, bgn);
                auto num = input.substr(bgn, end != std::string::npos ? end - bgn : end);
                try {
                    if constexpr(std::is_same_v<T, unsigned>) return (T) std::abs(std::stoi(num.data()));
                    if constexpr(std::is_same_v<T, int>) return (T) std::stoi(num.data());
                    if constexpr(std::is_same_v<T, long>) return (T) std::stol(num.data());
                    if constexpr(std::is_same_v<T, size_t>) return (T) std::abs(std::stol(num.data()));
                    if constexpr(std::is_same_v<T, double>) return (T) std::stod(num.data());
                    if constexpr(std::is_same<T, std::string>::value) return num;
                } catch(const std::exception &err) { std::printf("Failed to extract number from [%s]:%s", input.data(), err.what()); }
            }
        }
        return std::nullopt;
    }

    extern bool natcomp(std::string_view sa, std::string_view sb);

    template<typename T1, typename T2>
    std::optional<std::string> match(const T1 &v1, const T2 &v2) {
        for(const auto &s1 : v1)
            for(const auto &s2 : v2)
                if(s1 == s2) return std::string(s1);
        return std::nullopt;
    }

}
