#pragma once
#include <stdexcept>
#include <string_view>
enum class FileIdStatus { UPTODATE, STALE, MISSING };
enum class Model { SDUAL, MAJORANA, LBIT };
enum class StrictTableSize { TRUE, FALSE };
enum class SlabSelect { FULL, MIDCOL };
template<typename T>
requires std::is_enum_v<T>
constexpr std::string_view enum2str(const T &item) {
    if constexpr(std::is_same_v<T, FileIdStatus>) {
        if(item == FileIdStatus::UPTODATE) return "UPTODATE";
        if(item == FileIdStatus::STALE) return "STALE";
        if(item == FileIdStatus::MISSING) return "MISSING";
    }
    if constexpr(std::is_same_v<T, Model>) {
        if(item == Model::MAJORANA) return "MAJORANA";
        if(item == Model::SDUAL) return "SDUAL";
        if(item == Model::LBIT) return "LBIT";
    }
    if constexpr(std::is_same_v<T, StrictTableSize>) {
        if(item == StrictTableSize::TRUE) return "TRUE";
        if(item == StrictTableSize::FALSE) return "FALSE";
    }
    if constexpr(std::is_same_v<T, SlabSelect>) {
        if(item == SlabSelect::FULL) return "FULL";
        if(item == SlabSelect::MIDCOL) return "MIDCOL";
    }
    throw std::runtime_error("Given invalid enum item");
}

template<typename T>
constexpr auto str2enum(std::string_view item) {
    if constexpr(std::is_same_v<T, FileIdStatus>) {
        if(item == "UPTODATE") return FileIdStatus::UPTODATE;
        if(item == "STALE") return FileIdStatus::STALE;
        if(item == "MISSING") return FileIdStatus::MISSING;
    }
    if constexpr(std::is_same_v<T, Model>) {
        if(item == "xdmrg") return Model::SDUAL;
        if(item == "sdual") return Model::SDUAL;
        if(item == "SDUAL") return Model::SDUAL;
        if(item == "Sdual") return Model::SDUAL;
        if(item == "s-dual") return Model::SDUAL;
        if(item == "Majorana") return Model::MAJORANA;
        if(item == "majorana") return Model::MAJORANA;
        if(item == "ising-majorana") return Model::MAJORANA;
        if(item == "lbit") return Model::LBIT;
        if(item == "LBIT") return Model::LBIT;
        if(item == "flbit") return Model::LBIT;
        if(item == "f-lbit") return Model::LBIT;
        if(item == "FLBIT") return Model::LBIT;
    }
    if constexpr(std::is_same_v<T, StrictTableSize>) {
        if(item == "TRUE") return StrictTableSize::TRUE;
        if(item == "FALSE") return StrictTableSize::FALSE;
    }
    if constexpr(std::is_same_v<T, SlabSelect>) {
        if(item == "FULL") return SlabSelect::FULL;
        if(item == "MIDCOL") return SlabSelect::MIDCOL;
    }
    throw std::runtime_error("str2enum given invalid string item: " + std::string(item));
}