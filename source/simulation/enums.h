#pragma once
#include <h5pp/details/h5ppPermissions.h>
#include <stdexcept>
#include <string_view>
#include <type_traits>
enum class SimulationType { iDMRG, fDMRG, xDMRG, iTEBD };
enum class StorageLevel { NONE, LIGHT, NORMAL, FULL };
enum class StorageReason { JOURNAL, RESULTS, CHI_UPDATE, PROJ_STATE, INIT_STATE, EMIN_STATE, EMAX_STATE };
enum class StopReason { SUCCEEDED, SATURATED, MAX_ITERS, MAX_RESET, RANDOMIZE, NONE };
enum class ResumeReason {NOT_FINISHED, NOT_CONVERGED, NEW_STORAGE_LEVELS, NEW_SIMULATION_PARAMS, APPEND_STATES, NONE };
enum class FileCollisionPolicy { RESUME, BACKUP, RENAME, REPLACE};
enum class PerturbMode {
    PERCENTAGE,                // J_ptb = couplingPtb * J_rnd
    ABSOLUTE,                  // J_ptb = couplingPtb
    UNIFORM_RANDOM_PERCENTAGE, // J_ptb = std::random_uniform(-couplingPtb, couplingPtb) * J_rnd
    UNIFORM_RANDOM_ABSOLUTE,   // J_ptb = std::random_uniform(-couplingPtb, couplingPtb)
};

/* clang-format off */
template<typename T>
constexpr auto enum2str(const T &item) {
    if constexpr(std::is_same_v<T, SimulationType>) {
        if(item == SimulationType::iDMRG) return "iDMRG";
        if(item == SimulationType::fDMRG) return "fDMRG";
        if(item == SimulationType::xDMRG) return "xDMRG";
        if(item == SimulationType::iTEBD) return "iDMRG";
    }
    if constexpr(std::is_same_v<T, StopReason>) {
        if(item == StopReason::SUCCEEDED) return "SUCCEEDED";
        if(item == StopReason::SATURATED) return "SATURATED";
        if(item == StopReason::MAX_ITERS) return "MAX_ITERS";
        if(item == StopReason::MAX_RESET) return "MAX_RESET";
        if(item == StopReason::RANDOMIZE) return "RANDOMIZE";
        if(item == StopReason::NONE) return "NONE";
    }
    if constexpr(std::is_same_v<T, StorageLevel>) {
        if(item == StorageLevel::NONE)      return "NONE";
        if(item == StorageLevel::LIGHT)     return "LIGHT";
        if(item == StorageLevel::NORMAL)    return "NORMAL";
        if(item == StorageLevel::FULL)      return "FULL";
    }
    if constexpr(std::is_same_v<T, StorageReason>) {
        if(item == StorageReason::JOURNAL)      return "JOURNAL";
        if(item == StorageReason::RESULTS)      return "RESULTS";
        if(item == StorageReason::CHI_UPDATE)   return "CHI_UPDATE";
        if(item == StorageReason::PROJ_STATE)   return "PROJ_STATE";
        if(item == StorageReason::INIT_STATE)   return "INIT_STATE";
        if(item == StorageReason::EMIN_STATE)   return "EMIN_STATE";
        if(item == StorageReason::EMAX_STATE)   return "EMAX_STATE";
    }
    if constexpr(std::is_same_v<T, PerturbMode>) {
        if(item == PerturbMode::PERCENTAGE)                 return "PERCENTAGE";
        if(item == PerturbMode::ABSOLUTE)                   return "ABSOLUTE";
        if(item == PerturbMode::UNIFORM_RANDOM_PERCENTAGE)  return "UNIFORM_RANDOM_PERCENTAGE";
        if(item == PerturbMode::UNIFORM_RANDOM_ABSOLUTE)    return "UNIFORM_RANDOM_ABSOLUTE";
    }
    if constexpr(std::is_same_v<T, FileCollisionPolicy>) {
        if(item == FileCollisionPolicy::RESUME)     return "RESUME";
        if(item == FileCollisionPolicy::RENAME)     return "RENAME";
        if(item == FileCollisionPolicy::BACKUP)     return "BACKUP";
        if(item == FileCollisionPolicy::REPLACE)    return "REPLACE";
    }

    if constexpr(std::is_same_v<T, h5pp::FilePermission>) {
        if(item == h5pp::FilePermission::READONLY)          return "READONLY";
        if(item == h5pp::FilePermission::COLLISION_FAIL)    return "COLLISION_FAIL";
        if(item == h5pp::FilePermission::RENAME)            return "RENAME";
        if(item == h5pp::FilePermission::READWRITE)         return "READWRITE";
        if(item == h5pp::FilePermission::BACKUP)            return "BACKUP";
        if(item == h5pp::FilePermission::REPLACE)           return "REPLACE";
    }

    throw std::runtime_error("Given invalid enum item");
}

template<typename T>
constexpr auto str2enum(std::string_view item) {
    if constexpr(std::is_same_v<T, SimulationType>) {
        if(item == "iDMRG") return SimulationType::iDMRG;
        if(item == "fDMRG") return SimulationType::fDMRG;
        if(item == "xDMRG") return SimulationType::xDMRG;
        if(item == "iDMRG") return SimulationType::iTEBD;
    }
    if constexpr(std::is_same_v<T, StopReason>) {
        if(item == "SUCCEEDED") return StopReason::SUCCEEDED;
        if(item == "SATURATED") return StopReason::SATURATED;
        if(item == "MAX_ITERS") return StopReason::MAX_ITERS;
        if(item == "MAX_RESET") return StopReason::MAX_RESET;
        if(item == "RANDOMIZE") return StopReason::RANDOMIZE;
        if(item == "NONE")      return StopReason::NONE;
    }
    if constexpr(std::is_same_v<T, StorageLevel>) {
        if(item == "NONE")      return StorageLevel::NONE;
        if(item == "LIGHT")     return StorageLevel::LIGHT;
        if(item == "NORMAL")    return StorageLevel::NORMAL;
        if(item == "FULL")      return StorageLevel::FULL;
    }
    if constexpr(std::is_same_v<T, StorageReason>) {
        if(item == "JOURNAL")   return StorageReason::JOURNAL;
        if(item == "RESULTS")   return StorageReason::RESULTS;
        if(item == "CHI_UPDATE")return StorageReason::CHI_UPDATE;
        if(item == "PROJ_STATE")return StorageReason::PROJ_STATE;
        if(item == "INIT_STATE")return StorageReason::INIT_STATE;
        if(item == "EMIN_STATE")return StorageReason::EMIN_STATE;
        if(item == "EMAX_STATE")return StorageReason::EMAX_STATE;
    }
    if constexpr(std::is_same_v<T, PerturbMode>) {
        if(item == "PERCENTAGE")                return PerturbMode::PERCENTAGE;
        if(item == "ABSOLUTE")                  return PerturbMode::ABSOLUTE;
        if(item == "UNIFORM_RANDOM_PERCENTAGE") return PerturbMode::UNIFORM_RANDOM_PERCENTAGE;
        if(item == "UNIFORM_RANDOM_ABSOLUTE")   return PerturbMode::UNIFORM_RANDOM_ABSOLUTE;
    }
    if constexpr(std::is_same_v<T, FileCollisionPolicy>) {
        if(item == "RESUME")    return FileCollisionPolicy::RESUME;
        if(item == "RENAME")    return FileCollisionPolicy::RENAME;
        if(item == "BACKUP")    return FileCollisionPolicy::BACKUP;
        if(item == "REPLACE")   return FileCollisionPolicy::REPLACE;
    }
    if constexpr(std::is_same_v<T, h5pp::FilePermission>) {
        if(item == "READONLY")          return h5pp::FilePermission::READONLY;
        if(item == "COLLISION_FAIL")    return h5pp::FilePermission::COLLISION_FAIL;
        if(item == "RENAME")            return h5pp::FilePermission::RENAME;
        if(item == "READWRITE")         return h5pp::FilePermission::READWRITE;
        if(item == "BACKUP")            return h5pp::FilePermission::BACKUP;
        if(item == "REPLACE")           return h5pp::FilePermission::REPLACE;
    }

    throw std::runtime_error("Given invalid string item: " + std::string(item));
}
/* clang-format on */
