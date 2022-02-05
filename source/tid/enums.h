#pragma once
#include <stdexcept>
#include <string_view>

namespace tid {
    enum level : int {
        parent   = -1,
        normal   = 0,
        extra    = 1,
        detailed = 2,
    };

    constexpr std::string_view level2sv(level l) noexcept {
        switch(l) {
            case parent: return "parent  ";
            case normal: return "normal  ";
            case extra: return "extra   ";
            case detailed: return "detailed";
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
        case(tid::level::extra): return "extra";
        case(tid::level::detailed): return "detailed";
        case(tid::level::parent): return "parent";
        default: throw std::runtime_error("Unrecognized tid::level enum");
    }
}

template<>
constexpr auto sv2enum<tid::level>(std::string_view item) {
    if(item == "normal") return tid::level::normal;
    if(item == "extra") return tid::level::extra;
    if(item == "detailed") return tid::level::detailed;
    if(item == "parent") return tid::level::parent;
    throw std::runtime_error("Given item is not a tid::level enum");
}
