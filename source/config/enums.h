#pragma once
#include <h5pp/details/h5ppEnums.h>
#include <stdexcept>
#include <string_view>
#include <type_traits>
enum class AlgorithmType { iDMRG, fDMRG, xDMRG, iTEBD };
enum class MultisiteMove { ONE, MID, MAX };
enum class StateRitz { LR, SR }; // Smallest Real or Largest Real, i.e. ground state or max state. Relevant for fdmrg.
enum class ModelType { ising_tf_rf, ising_sdual };
enum class StorageLevel { NONE, LIGHT, NORMAL, FULL };
enum class StorageReason { CHECKPOINT, FINISHED, CHI_UPDATE, PROJ_STATE, INIT_STATE, EMIN_STATE, EMAX_STATE, MODEL };
enum class StopReason { SUCCEEDED, SATURATED, MAX_ITERS, MAX_RESET, RANDOMIZE, NONE };
enum class ResetReason { INIT, FIND_WINDOW, SATURATED, NEW_STATE, CHI_UPDATE};
enum class NormPolicy { ALWAYS, IFNEEDED }; // Rules of engagement
enum class FileCollisionPolicy { RESUME, BACKUP, RENAME, REPLACE };
enum class StateType {
    RANDOM_PRODUCT_STATE,
    RANDOM_ENTANGLED_STATE,
    RANDOMIZE_GIVEN_STATE,
    PRODUCT_STATE
};


enum class PerturbMode {
    PERCENTAGE,                // J_ptb = couplingPtb * J_rnd
    ABSOLUTE,                  // J_ptb = couplingPtb
    UNIFORM_RANDOM_PERCENTAGE, // J_ptb = std::random_uniform(-couplingPtb, couplingPtb) * J_rnd
    UNIFORM_RANDOM_ABSOLUTE,   // J_ptb = std::random_uniform(-couplingPtb, couplingPtb)
};



enum class fdmrg_task {
    INIT_RANDOMIZE_MODEL,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE,
    INIT_RANDOMIZE_INTO_ENTANGLED_STATE,
    INIT_BOND_DIM_LIMITS,
    INIT_WRITE_MODEL,
    INIT_CLEAR_STATUS,
    INIT_DEFAULT,
    FIND_GROUND_STATE,
    FIND_HIGHEST_STATE,
    POST_WRITE_RESULT,
    POST_PRINT_RESULT,
    POST_DEFAULT,
};

enum class xdmrg_task {
    INIT_RANDOMIZE_MODEL,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE,
    INIT_RANDOMIZE_INTO_ENTANGLED_STATE,
    INIT_RANDOMIZE_INTO_STATE_IN_WIN,
    INIT_RANDOMIZE_FROM_CURRENT_STATE,
    INIT_BOND_DIM_LIMITS,
    INIT_ENERGY_LIMITS,
    INIT_WRITE_MODEL,
    INIT_CLEAR_STATUS,
    INIT_DEFAULT,
    NEXT_RANDOMIZE_INTO_PRODUCT_STATE,
    NEXT_RANDOMIZE_INTO_ENTANGLED_STATE,
    NEXT_RANDOMIZE_FROM_CURRENT_STATE,
    NEXT_RANDOMIZE_INTO_STATE_IN_WIN,
    FIND_ENERGY_RANGE,
    FIND_EXCITED_STATE,
    POST_WRITE_RESULT,
    POST_PRINT_RESULT,
    POST_DEFAULT,
};

/* clang-format off */
template<typename T>
constexpr std::string_view enum2str(const T &item) {
    if constexpr(std::is_same_v<T, AlgorithmType>) {
        if(item == AlgorithmType::iDMRG) return "iDMRG";
        if(item == AlgorithmType::fDMRG) return "fDMRG";
        if(item == AlgorithmType::xDMRG) return "xDMRG";
        if(item == AlgorithmType::iTEBD) return "iDMRG";
    }
    if constexpr(std::is_same_v<T, MultisiteMove>) {
        if(item == MultisiteMove::ONE) return "ONE";
        if(item == MultisiteMove::MID) return "MID";
        if(item == MultisiteMove::MAX) return "MAX";
    }
    if constexpr(std::is_same_v<T, StateRitz>) {
        if(item == StateRitz::SR) return "SR";
        if(item == StateRitz::LR) return "LR";
    }
    if constexpr(std::is_same_v<T, ModelType>) {
        if(item == ModelType::ising_tf_rf) return "ising_tf_rf";
        if(item == ModelType::ising_sdual) return "ising_sdual";
    }
    if constexpr(std::is_same_v<T, StopReason>) {
        if(item == StopReason::SUCCEEDED) return "SUCCEEDED";
        if(item == StopReason::SATURATED) return "SATURATED";
        if(item == StopReason::MAX_ITERS) return "MAX_ITERS";
        if(item == StopReason::MAX_RESET) return "MAX_RESET";
        if(item == StopReason::RANDOMIZE) return "RANDOMIZE";
        if(item == StopReason::NONE)      return "NONE";
    }
    if constexpr(std::is_same_v<T, ResetReason>) {
        if(item == ResetReason::INIT) return "INIT";
        if(item == ResetReason::FIND_WINDOW)        return "FIND_WINDOW";
        if(item == ResetReason::SATURATED)          return "SATURATED";
        if(item == ResetReason::NEW_STATE)  return "NEW_STATE";
        if(item == ResetReason::CHI_UPDATE)         return "CHI_UPDATE";
    }
    if constexpr(std::is_same_v<T, NormPolicy>) {
        if(item == NormPolicy::ALWAYS)      return "ALWAYS";
        if(item == NormPolicy::IFNEEDED)    return "IFNEEDED";
    }
    if constexpr(std::is_same_v<T, StorageLevel>) {
        if(item == StorageLevel::NONE)      return "NONE";
        if(item == StorageLevel::LIGHT)     return "LIGHT";
        if(item == StorageLevel::NORMAL)    return "NORMAL";
        if(item == StorageLevel::FULL)      return "FULL";
    }
    if constexpr(std::is_same_v<T, StorageReason>) {
        if(item == StorageReason::CHECKPOINT)   return "CHECKPOINT";
        if(item == StorageReason::FINISHED)     return "FINISHED";
        if(item == StorageReason::CHI_UPDATE)   return "CHI_UPDATE";
        if(item == StorageReason::PROJ_STATE)   return "PROJ_STATE";
        if(item == StorageReason::INIT_STATE)   return "INIT_STATE";
        if(item == StorageReason::EMIN_STATE)   return "EMIN_STATE";
        if(item == StorageReason::EMAX_STATE)   return "EMAX_STATE";
        if(item == StorageReason::MODEL)        return "MODEL";
    }
    if constexpr(std::is_same_v<T, StateType>) {
        if(item == StateType::RANDOM_PRODUCT_STATE)     return "RANDOM_PRODUCT_STATE";
        if(item == StateType::RANDOM_ENTANGLED_STATE)   return "RANDOM_ENTANGLED_STATE";
        if(item == StateType::RANDOMIZE_GIVEN_STATE)    return "RANDOMIZE_GIVEN_STATE";
        if(item == StateType::PRODUCT_STATE)            return "PRODUCT_STATE";
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
    if constexpr(std::is_same_v<T, AlgorithmType>) {
        if(item == "iDMRG") return AlgorithmType::iDMRG;
        if(item == "fDMRG") return AlgorithmType::fDMRG;
        if(item == "xDMRG") return AlgorithmType::xDMRG;
        if(item == "iDMRG") return AlgorithmType::iTEBD;
    }
    if constexpr(std::is_same_v<T, MultisiteMove>) {
        if(item == "ONE") return MultisiteMove::ONE;
        if(item == "MID") return MultisiteMove::MID;
        if(item == "MAX") return MultisiteMove::MAX;
    }
    if constexpr(std::is_same_v<T, StateRitz>) {
        if(item == "SR") return StateRitz::SR;
        if(item == "LR") return StateRitz::LR;
    }
    if constexpr(std::is_same_v<T, ModelType>) {
        if(item == "ising_tf_rf") return ModelType::ising_tf_rf;
        if(item == "ising_sdual") return ModelType::ising_sdual;
    }
    if constexpr(std::is_same_v<T, StopReason>) {
        if(item == "SUCCEEDED") return StopReason::SUCCEEDED;
        if(item == "SATURATED") return StopReason::SATURATED;
        if(item == "MAX_ITERS") return StopReason::MAX_ITERS;
        if(item == "MAX_RESET") return StopReason::MAX_RESET;
        if(item == "RANDOMIZE") return StopReason::RANDOMIZE;
        if(item == "NONE")      return StopReason::NONE;
    }
    if constexpr(std::is_same_v<T, ResetReason>) {
        if(item == "INIT")              return ResetReason::INIT;
        if(item == "FIND_WINDOW")       return ResetReason::FIND_WINDOW;
        if(item == "SATURATED")         return ResetReason::SATURATED;
        if(item == "NEW_STATE") return ResetReason::NEW_STATE;
        if(item == "CHI_UPDATE")        return ResetReason::CHI_UPDATE;
    }
    if constexpr(std::is_same_v<T, NormPolicy>) {
        if(item == "ALWAYS") return NormPolicy::ALWAYS;
        if(item == "IFNEEDED") return NormPolicy::IFNEEDED;
    }
    if constexpr(std::is_same_v<T, StorageLevel>) {
        if(item == "NONE")      return StorageLevel::NONE;
        if(item == "LIGHT")     return StorageLevel::LIGHT;
        if(item == "NORMAL")    return StorageLevel::NORMAL;
        if(item == "FULL")      return StorageLevel::FULL;
    }
    if constexpr(std::is_same_v<T, StorageReason>) {
        if(item == "CHECKPOINT")return StorageReason::CHECKPOINT;
        if(item == "FINISHED")  return StorageReason::FINISHED;
        if(item == "CHI_UPDATE")return StorageReason::CHI_UPDATE;
        if(item == "PROJ_STATE")return StorageReason::PROJ_STATE;
        if(item == "INIT_STATE")return StorageReason::INIT_STATE;
        if(item == "EMIN_STATE")return StorageReason::EMIN_STATE;
        if(item == "EMAX_STATE")return StorageReason::EMAX_STATE;
        if(item == "MODEL")     return StorageReason::MODEL;
    }
    if constexpr(std::is_same_v<T, StateType>) {
        if(item == "RANDOM_PRODUCT_STATE")     return StateType::RANDOM_PRODUCT_STATE;
        if(item == "RANDOM_ENTANGLED_STATE")   return StateType::RANDOM_ENTANGLED_STATE;
        if(item == "RANDOMIZE_GIVEN_STATE")    return StateType::RANDOMIZE_GIVEN_STATE;
        if(item == "PRODUCT_STATE")            return StateType::PRODUCT_STATE;
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

    throw std::runtime_error("str2enum given invalid string item: " + std::string(item));
}
/* clang-format on */
