#pragma once
#include <h5pp/details/h5ppPermissions.h>
#include <type_traits>
#include <string_view>
#include <stdexcept>
enum class SimulationType { iDMRG, fDMRG, xDMRG, iTEBD };
enum class StorageLevel { NONE, LIGHT, NORMAL, FULL };
enum class PerturbMode {
    PERCENTAGE,                // J_ptb = couplingPtb * J_rnd
    ABSOLUTE,                  // J_ptb = couplingPtb
    UNIFORM_RANDOM_PERCENTAGE, // J_ptb = std::random_uniform(-couplingPtb, couplingPtb) * J_rnd
    UNIFORM_RANDOM_ABSOLUTE,   // J_ptb = std::random_uniform(-couplingPtb, couplingPtb)
};

template<typename T >
constexpr auto enum2str(T & item) {
    if constexpr (std::is_same_v<T,SimulationType>){
        if (item == SimulationType::iDMRG) return "iDMRG";
        if (item == SimulationType::fDMRG) return "fDMRG";
        if (item == SimulationType::xDMRG) return "xDMRG";
        if (item == SimulationType::iTEBD) return "iDMRG";
    }
    if constexpr (std::is_same_v<T,StorageLevel>) {
        if (item == StorageLevel::NONE)     return "NONE";
        if (item == StorageLevel::LIGHT)    return "LIGHT";
        if (item == StorageLevel::NORMAL)   return "NORMAL";
        if (item == StorageLevel::FULL)     return "FULL";
    }
    if constexpr (std::is_same_v<T,PerturbMode>){
        if (item == PerturbMode::PERCENTAGE               ) return "PERCENTAGE";
        if (item == PerturbMode::ABSOLUTE                 ) return "ABSOLUTE";
        if (item == PerturbMode::UNIFORM_RANDOM_PERCENTAGE) return "UNIFORM_RANDOM_PERCENTAGE";
        if (item == PerturbMode::UNIFORM_RANDOM_ABSOLUTE  ) return "UNIFORM_RANDOM_ABSOLUTE";
    }
    if constexpr (std::is_same_v<T,h5pp::AccessMode>){
        if (item == h5pp::AccessMode::READWRITE) return "READWRITE";
        if (item == h5pp::AccessMode::READONLY) return "READONLY";
    }
    if constexpr (std::is_same_v<T,h5pp::CreateMode>){
        if (item == h5pp::CreateMode::OPEN) return "OPEN";
        if (item == h5pp::CreateMode::TRUNCATE) return "TRUNCATE";
        if (item == h5pp::CreateMode::RENAME) return "RENAME";
    }

    throw std::runtime_error("Given invalid enum item");
}

template<typename T>
constexpr auto str2enum(std::string_view item) {
    if constexpr (std::is_same_v<T,SimulationType>){
        if (item == "iDMRG") return SimulationType::iDMRG;
        if (item == "fDMRG") return SimulationType::fDMRG;
        if (item == "xDMRG") return SimulationType::xDMRG;
        if (item == "iDMRG") return SimulationType::iTEBD;
    }
    if constexpr (std::is_same_v<T,StorageLevel>) {
        if (item == "NONE"   ) return StorageLevel::NONE;
        if (item == "LIGHT"  ) return StorageLevel::LIGHT;
        if (item == "NORMAL" ) return StorageLevel::NORMAL;
        if (item == "FULL"   ) return StorageLevel::FULL;
    }
    if constexpr (std::is_same_v<T,PerturbMode>){
        if (item == "PERCENTAGE"                ) return  PerturbMode::PERCENTAGE;
        if (item == "ABSOLUTE"                  ) return  PerturbMode::ABSOLUTE;
        if (item == "UNIFORM_RANDOM_PERCENTAGE" ) return  PerturbMode::UNIFORM_RANDOM_PERCENTAGE;
        if (item == "UNIFORM_RANDOM_ABSOLUTE"   ) return  PerturbMode::UNIFORM_RANDOM_ABSOLUTE;
    }
    if constexpr (std::is_same_v<T,h5pp::AccessMode>){
        if (item == "READWRITE") return  h5pp::AccessMode::READWRITE;
        if (item == "READONLY")  return  h5pp::AccessMode::READONLY ;
    }
    if constexpr (std::is_same_v<T,h5pp::CreateMode>){
        if (item == "OPEN")     return h5pp::CreateMode::OPEN;
        if (item == "TRUNCATE") return h5pp::CreateMode::TRUNCATE;
        if (item == "RENAME")   return h5pp::CreateMode::RENAME;
    }
    throw std::runtime_error("Given invalid string item");
}