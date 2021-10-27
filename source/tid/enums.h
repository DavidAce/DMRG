#pragma once
#include <stdexcept>
#include <string_view>

namespace tid {
    enum level : int {
        parent = -1,
        normal = 0,
        detail = 1,
        pedant = 2,
    };

    constexpr std::string_view level2sv(level l) noexcept {
        switch(l) {
            case normal: return "normal";
            case detail: return "detail";
            case pedant: return "pedant";
            case parent: return "parent";
            default: return "unknown";
        }
    }
}

template<typename T>
extern constexpr std::string_view enum2sv(const T &item);

template<typename T>
extern constexpr auto sv2enum(std::string_view item);

template<>
constexpr std::string_view enum2sv(const tid::level &item) {
    switch(item) {
        case(tid::level::normal): return "normal";
        case(tid::level::detail): return "detail";
        case(tid::level::pedant): return "pedant";
        case(tid::level::parent): return "parent";
        default: throw std::runtime_error("Unrecognized tid::level enum");
    }
}

template<>
constexpr auto sv2enum<tid::level>(std::string_view item) {
    if(item == "normal") return tid::level::normal;
    if(item == "detail") return tid::level::detail;
    if(item == "pedant") return tid::level::pedant;
    if(item == "parent") return tid::level::parent;
    throw std::runtime_error("Given item is not a tid::level enum");
}
