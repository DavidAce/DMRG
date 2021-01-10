#pragma once
#include <stdexcept>
#include <string_view>
#include <type_traits>
enum class AlgorithmType { iDMRG, fDMRG, xDMRG, iTEBD, fLBIT, ANY };
enum class MultisiteMove { ONE, MID, MAX };
enum class StateRitz { LR, SR }; // Smallest Real or Largest Real, i.e. ground state or max state. Relevant for fdmrg.
enum class SVDMode { EIGEN, LAPACKE };
enum class ModelType { ising_tf_rf, ising_sdual, lbit };
enum class EdgeStatus {STALE, FRESH};
enum class StorageLevel { NONE, LIGHT, NORMAL, FULL };
enum class StorageReason { SAVEPOINT, CHECKPOINT, FINISHED, CHI_UPDATE, PROJ_STATE, INIT_STATE, EMIN_STATE, EMAX_STATE, MODEL };
enum class CopyPolicy { FORCE, TRY, OFF };
enum class StopReason { SUCCEEDED, SATURATED, MAX_ITERS, MAX_RESET, RANDOMIZE, NONE };
enum class ResetReason { INIT, FIND_WINDOW, SATURATED, NEW_STATE, CHI_UPDATE };
enum class NormPolicy { ALWAYS, IFNEEDED }; // Rules of engagement
enum class FileCollisionPolicy { RESUME, BACKUP, RENAME, REPLACE };
enum class FileResumePolicy { FULL, FAST };
enum class LogPolicy { NORMAL, QUIET };
enum class RandomizerMode {SHUFFLE, SELECT1, ASIS};
enum class OptType { REAL, CPLX };
enum class OptMode { VARIANCE, OVERLAP };
enum class OptSpace { SUBSPACE_ONLY, SUBSPACE_AND_DIRECT, DIRECT };
enum class OptWhen {
    ALWAYS,
    NEVER,
    PREV_FAIL,
};
enum class OptMark {
    PASS,
    FAIL,
};
enum class OptInit { CURRENT_STATE, LAST_RESULT };
enum class StateInitType { REAL, CPLX };
enum class StateInit {
    RANDOM_PRODUCT_STATE,
    RANDOM_ENTANGLED_STATE,
    RANDOMIZE_PREVIOUS_STATE,
    PRODUCT_STATE_ALIGNED,
    PRODUCT_STATE_NEEL,
    //    ALL_DOWN_ONE_UP,
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
    POST_PRINT_PROFILING,
    POST_DEFAULT,
    PROF_RESET,
};

enum class flbit_task {
    INIT_RANDOMIZE_MODEL,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE,
    INIT_RANDOMIZE_INTO_ENTANGLED_STATE,
    INIT_BOND_DIM_LIMITS,
    INIT_WRITE_MODEL,
    INIT_CLEAR_STATUS,
    INIT_DEFAULT,
    TIME_EVOLVE,
    POST_WRITE_RESULT,
    POST_PRINT_RESULT,
    POST_PRINT_PROFILING,
    POST_DEFAULT,
    PROF_RESET,
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
    NEXT_RANDOMIZE_PREVIOUS_STATE,
    NEXT_RANDOMIZE_INTO_STATE_IN_WIN,
    FIND_ENERGY_RANGE,
    FIND_EXCITED_STATE,
    POST_WRITE_RESULT,
    POST_PRINT_RESULT,
    POST_PRINT_PROFILING,
    POST_DEFAULT,
    PROF_RESET,
};

/* clang-format off */
template<typename T>
constexpr std::string_view enum2str(const T &item) {
    if constexpr(std::is_same_v<T, AlgorithmType>) {
        if(item == AlgorithmType::iDMRG)                                return "iDMRG";
        if(item == AlgorithmType::fDMRG)                                return "fDMRG";
        if(item == AlgorithmType::fLBIT)                                return "fLBIT";
        if(item == AlgorithmType::xDMRG)                                return "xDMRG";
        if(item == AlgorithmType::iTEBD)                                return "iTEBD";
        if(item == AlgorithmType::ANY)                                  return "ANY";
    }
    if constexpr(std::is_same_v<T, MultisiteMove>) {
        if(item == MultisiteMove::ONE)                                  return "ONE";
        if(item == MultisiteMove::MID)                                  return "MID";
        if(item == MultisiteMove::MAX)                                  return "MAX";
    }
    if constexpr(std::is_same_v<T, StateRitz>) {
        if(item == StateRitz::SR)                                       return "SR";
        if(item == StateRitz::LR)                                       return "LR";
    }
    if constexpr(std::is_same_v<T, SVDMode>) {
        if(item == SVDMode::EIGEN)                                      return "EIGEN";
        if(item == SVDMode::LAPACKE)                                    return "LAPACKE";
    }
    if constexpr(std::is_same_v<T, ModelType>) {
        if(item == ModelType::ising_tf_rf)                              return "ising_tf_rf";
        if(item == ModelType::ising_sdual)                              return "ising_sdual";
        if(item == ModelType::lbit)                                     return "lbit";
    }
    if constexpr(std::is_same_v<T, EdgeStatus>) {
        if(item == EdgeStatus::STALE)                                   return "STALE";
        if(item == EdgeStatus::FRESH)                                   return "FRESH";
    }
    if constexpr(std::is_same_v<T, StopReason>) {
        if(item == StopReason::SUCCEEDED)                               return "SUCCEEDED";
        if(item == StopReason::SATURATED)                               return "SATURATED";
        if(item == StopReason::MAX_ITERS)                               return "MAX_ITERS";
        if(item == StopReason::MAX_RESET)                               return "MAX_RESET";
        if(item == StopReason::RANDOMIZE)                               return "RANDOMIZE";
        if(item == StopReason::NONE)                                    return "NONE";
    }
    if constexpr(std::is_same_v<T, ResetReason>) {
        if(item == ResetReason::INIT)                                   return "INIT";
        if(item == ResetReason::FIND_WINDOW)                            return "FIND_WINDOW";
        if(item == ResetReason::SATURATED)                              return "SATURATED";
        if(item == ResetReason::NEW_STATE)                              return "NEW_STATE";
        if(item == ResetReason::CHI_UPDATE)                             return "CHI_UPDATE";
    }
    if constexpr(std::is_same_v<T, NormPolicy>) {
        if(item == NormPolicy::ALWAYS)                                  return "ALWAYS";
        if(item == NormPolicy::IFNEEDED)                                return "IFNEEDED";
    }
    if constexpr(std::is_same_v<T, LogPolicy>) {
        if(item == LogPolicy::NORMAL)                                   return "NORMAL";
        if(item == LogPolicy::QUIET)                                    return "QUIET";
    }
    if constexpr(std::is_same_v<T, StorageLevel>) {
        if(item == StorageLevel::NONE)                                  return "NONE";
        if(item == StorageLevel::LIGHT)                                 return "LIGHT";
        if(item == StorageLevel::NORMAL)                                return "NORMAL";
        if(item == StorageLevel::FULL)                                  return "FULL";
    }
    if constexpr(std::is_same_v<T, StorageReason>) {
        if(item == StorageReason::SAVEPOINT)                            return "SAVEPOINT";
        if(item == StorageReason::CHECKPOINT)                           return "CHECKPOINT";
        if(item == StorageReason::FINISHED)                             return "FINISHED";
        if(item == StorageReason::CHI_UPDATE)                           return "CHI_UPDATE";
        if(item == StorageReason::PROJ_STATE)                           return "PROJ_STATE";
        if(item == StorageReason::INIT_STATE)                           return "INIT_STATE";
        if(item == StorageReason::EMIN_STATE)                           return "EMIN_STATE";
        if(item == StorageReason::EMAX_STATE)                           return "EMAX_STATE";
        if(item == StorageReason::MODEL)                                return "MODEL";
    }
    if constexpr(std::is_same_v<T, StateInit>) {
        if(item == StateInit::RANDOM_PRODUCT_STATE)                     return "RANDOM_PRODUCT_STATE";
        if(item == StateInit::RANDOM_ENTANGLED_STATE)                   return "RANDOM_ENTANGLED_STATE";
        if(item == StateInit::RANDOMIZE_PREVIOUS_STATE)                 return "RANDOMIZE_PREVIOUS_STATE";
        if(item == StateInit::PRODUCT_STATE_ALIGNED)                    return "PRODUCT_STATE_ALIGNED";
        if(item == StateInit::PRODUCT_STATE_NEEL)                       return "PRODUCT_STATE_NEEL";
    }
    if constexpr(std::is_same_v<T, StateInitType>) {
        if(item == StateInitType::REAL)   return "REAL";
        if(item == StateInitType::CPLX)   return "CPLX";
    }
    if constexpr(std::is_same_v<T, PerturbMode>) {
        if(item == PerturbMode::PERCENTAGE)                             return "PERCENTAGE";
        if(item == PerturbMode::ABSOLUTE)                               return "ABSOLUTE";
        if(item == PerturbMode::UNIFORM_RANDOM_PERCENTAGE)              return "UNIFORM_RANDOM_PERCENTAGE";
        if(item == PerturbMode::UNIFORM_RANDOM_ABSOLUTE)                return "UNIFORM_RANDOM_ABSOLUTE";
    }
    if constexpr(std::is_same_v<T, FileCollisionPolicy>) {
        if(item == FileCollisionPolicy::RESUME)                         return "RESUME";
        if(item == FileCollisionPolicy::RENAME)                         return "RENAME";
        if(item == FileCollisionPolicy::BACKUP)                         return "BACKUP";
        if(item == FileCollisionPolicy::REPLACE)                        return "REPLACE";
    }
    if constexpr(std::is_same_v<T, FileResumePolicy>) {
        if(item == FileResumePolicy::FULL)                              return "FULL";
        if(item == FileResumePolicy::FAST)                              return "FAST";
    }
    if constexpr(std::is_same_v<T,fdmrg_task>){
        if(item == fdmrg_task::INIT_RANDOMIZE_MODEL)                    return "INIT_RANDOMIZE_MODEL";
        if(item == fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE)       return "INIT_RANDOMIZE_INTO_PRODUCT_STATE";
        if(item == fdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE)     return "INIT_RANDOMIZE_INTO_ENTANGLED_STATE";
        if(item == fdmrg_task::INIT_BOND_DIM_LIMITS)                    return "INIT_BOND_DIM_LIMITS";
        if(item == fdmrg_task::INIT_WRITE_MODEL)                        return "INIT_WRITE_MODEL";
        if(item == fdmrg_task::INIT_CLEAR_STATUS)                       return "INIT_CLEAR_STATUS";
        if(item == fdmrg_task::INIT_DEFAULT)                            return "INIT_DEFAULT";
        if(item == fdmrg_task::FIND_GROUND_STATE)                       return "FIND_GROUND_STATE";
        if(item == fdmrg_task::FIND_HIGHEST_STATE)                      return "FIND_HIGHEST_STATE";
        if(item == fdmrg_task::POST_WRITE_RESULT)                       return "POST_WRITE_RESULT";
        if(item == fdmrg_task::POST_PRINT_RESULT)                       return "POST_PRINT_RESULT";
        if(item == fdmrg_task::POST_DEFAULT)                            return "POST_DEFAULT";
        if(item == fdmrg_task::PROF_RESET)                              return "PROF_RESET";

    }

    if constexpr(std::is_same_v<T,xdmrg_task>){
        if(item == xdmrg_task::INIT_RANDOMIZE_MODEL)                   return "INIT_RANDOMIZE_MODEL";
        if(item == xdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE)      return "INIT_RANDOMIZE_INTO_PRODUCT_STATE";
        if(item == xdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE)    return "INIT_RANDOMIZE_INTO_ENTANGLED_STATE";
        if(item == xdmrg_task::INIT_RANDOMIZE_INTO_STATE_IN_WIN)       return "INIT_RANDOMIZE_INTO_STATE_IN_WIN";
        if(item == xdmrg_task::INIT_RANDOMIZE_FROM_CURRENT_STATE)      return "INIT_RANDOMIZE_FROM_CURRENT_STATE";
        if(item == xdmrg_task::INIT_BOND_DIM_LIMITS)                   return "INIT_BOND_DIM_LIMITS";
        if(item == xdmrg_task::INIT_ENERGY_LIMITS)                     return "INIT_ENERGY_LIMITS";
        if(item == xdmrg_task::INIT_WRITE_MODEL)                       return "INIT_WRITE_MODEL";
        if(item == xdmrg_task::INIT_CLEAR_STATUS)                      return "INIT_CLEAR_STATUS";
        if(item == xdmrg_task::INIT_DEFAULT)                           return "INIT_DEFAULT";
        if(item == xdmrg_task::NEXT_RANDOMIZE_INTO_PRODUCT_STATE)      return "NEXT_RANDOMIZE_INTO_PRODUCT_STATE";
        if(item == xdmrg_task::NEXT_RANDOMIZE_INTO_ENTANGLED_STATE)    return "NEXT_RANDOMIZE_INTO_ENTANGLED_STATE";
        if(item == xdmrg_task::NEXT_RANDOMIZE_PREVIOUS_STATE)          return "NEXT_RANDOMIZE_PREVIOUS_STATE";
        if(item == xdmrg_task::NEXT_RANDOMIZE_INTO_STATE_IN_WIN)       return "NEXT_RANDOMIZE_INTO_STATE_IN_WIN";
        if(item == xdmrg_task::FIND_ENERGY_RANGE)                      return "FIND_ENERGY_RANGE";
        if(item == xdmrg_task::FIND_EXCITED_STATE)                     return "FIND_EXCITED_STATE";
        if(item == xdmrg_task::POST_WRITE_RESULT)                      return "POST_WRITE_RESULT";
        if(item == xdmrg_task::POST_PRINT_RESULT)                      return "POST_PRINT_RESULT";
        if(item == xdmrg_task::POST_DEFAULT)                           return "POST_DEFAULT";
        if(item == xdmrg_task::PROF_RESET)                             return "PROF_RESET";
    }
    if constexpr(std::is_same_v<T,CopyPolicy>){
        if(item == CopyPolicy::FORCE)                                  return "FORCE";
        if(item == CopyPolicy::TRY)                                    return "TRY";
        if(item == CopyPolicy::OFF)                                    return "OFF";
    }
    if constexpr(std::is_same_v<T,RandomizerMode>){
        if(item == RandomizerMode::SHUFFLE)                            return "SHUFFLE";
        if(item == RandomizerMode::SELECT1)                            return "SELECT1";
        if(item == RandomizerMode::ASIS)                               return "ASIS";
    }
    if constexpr(std::is_same_v<T,OptType>){
        if(item == OptType::REAL)                                      return "REAL";
        if(item == OptType::CPLX)                                      return "CPLX";
    }
    if constexpr(std::is_same_v<T,OptMode>){
        if(item == OptMode::VARIANCE)                                  return "VARIANCE";
        if(item == OptMode::OVERLAP)                                   return "OVERLAP";
    }
    if constexpr(std::is_same_v<T,OptSpace>){
        if(item == OptSpace::SUBSPACE_ONLY)                            return "SUBSPACE_ONLY";
        if(item == OptSpace::SUBSPACE_AND_DIRECT)                      return "SUBSPACE_AND_DIRECT";
        if(item == OptSpace::DIRECT)                                   return "DIRECT";
    }
    if constexpr(std::is_same_v<T,OptWhen>){
        if(item == OptWhen::ALWAYS)                                    return "ALWAYS";
        if(item == OptWhen::NEVER)                                     return "NEVER";
        if(item == OptWhen::PREV_FAIL)                                 return "PREV_FAIL";
    }
    if constexpr(std::is_same_v<T,OptMark>){
        if(item == OptMark::PASS)                                      return "PASS";
        if(item == OptMark::FAIL)                                      return "FAIL";
    }
    if constexpr(std::is_same_v<T,OptInit>){
        if(item == OptInit::CURRENT_STATE)                             return "CURRENT_STATE";
        if(item == OptInit::LAST_RESULT)                               return "LAST_RESULT";
    }
    throw std::runtime_error("Given invalid enum item");
}

template<typename T>
constexpr auto str2enum(std::string_view item) {
    if constexpr(std::is_same_v<T, AlgorithmType>) {
        if(item == "iDMRG")                                 return AlgorithmType::iDMRG;
        if(item == "fDMRG")                                 return AlgorithmType::fDMRG;
        if(item == "fLBIT")                                 return AlgorithmType::fLBIT;
        if(item == "xDMRG")                                 return AlgorithmType::xDMRG;
        if(item == "iTEBD")                                 return AlgorithmType::iTEBD;
        if(item == "ANY")                                   return AlgorithmType::ANY;
    }
    if constexpr(std::is_same_v<T, MultisiteMove>) {
        if(item == "ONE")                                   return MultisiteMove::ONE;
        if(item == "MID")                                   return MultisiteMove::MID;
        if(item == "MAX")                                   return MultisiteMove::MAX;
    }
    if constexpr(std::is_same_v<T, StateRitz>) {
        if(item == "SR")                                    return StateRitz::SR;
        if(item == "LR")                                    return StateRitz::LR;
    }
    if constexpr(std::is_same_v<T, SVDMode>) {
        if(item == "EIGEN")                                 return SVDMode::EIGEN;
        if(item == "LAPACKE")                               return SVDMode::LAPACKE;
    }
    if constexpr(std::is_same_v<T, ModelType>) {
        if(item == "ising_tf_rf")                           return ModelType::ising_tf_rf;
        if(item == "ising_sdual")                           return ModelType::ising_sdual;
        if(item == "lbit")                                  return ModelType::lbit;
    }
    if constexpr(std::is_same_v<T, EdgeStatus>) {
        if(item == "STALE")                                 return EdgeStatus::STALE ;
        if(item == "FRESH")                                 return EdgeStatus::FRESH ;
    }
    if constexpr(std::is_same_v<T, StopReason>) {
        if(item == "SUCCEEDED")                             return StopReason::SUCCEEDED;
        if(item == "SATURATED")                             return StopReason::SATURATED;
        if(item == "MAX_ITERS")                             return StopReason::MAX_ITERS;
        if(item == "MAX_RESET")                             return StopReason::MAX_RESET;
        if(item == "RANDOMIZE")                             return StopReason::RANDOMIZE;
        if(item == "NONE")                                  return StopReason::NONE;
    }
    if constexpr(std::is_same_v<T, ResetReason>) {
        if(item == "INIT")                                  return ResetReason::INIT;
        if(item == "FIND_WINDOW")                           return ResetReason::FIND_WINDOW;
        if(item == "SATURATED")                             return ResetReason::SATURATED;
        if(item == "NEW_STATE")                             return ResetReason::NEW_STATE;
        if(item == "CHI_UPDATE")                            return ResetReason::CHI_UPDATE;
    }
    if constexpr(std::is_same_v<T, NormPolicy>) {
        if(item == "ALWAYS")                                return NormPolicy::ALWAYS;
        if(item == "IFNEEDED")                              return NormPolicy::IFNEEDED;
    }
    if constexpr(std::is_same_v<T, LogPolicy>) {
        if(item == "NORMAL")                                return LogPolicy::NORMAL;
        if(item == "QUIET")                                 return LogPolicy::QUIET;
    }
    if constexpr(std::is_same_v<T, StorageLevel>) {
        if(item == "NONE")                                  return StorageLevel::NONE;
        if(item == "LIGHT")                                 return StorageLevel::LIGHT;
        if(item == "NORMAL")                                return StorageLevel::NORMAL;
        if(item == "FULL")                                  return StorageLevel::FULL;
    }
    if constexpr(std::is_same_v<T, StorageReason>) {
        if(item == "SAVEPOINT")                             return StorageReason::SAVEPOINT;
        if(item == "CHECKPOINT")                            return StorageReason::CHECKPOINT;
        if(item == "FINISHED")                              return StorageReason::FINISHED;
        if(item == "CHI_UPDATE")                            return StorageReason::CHI_UPDATE;
        if(item == "PROJ_STATE")                            return StorageReason::PROJ_STATE;
        if(item == "INIT_STATE")                            return StorageReason::INIT_STATE;
        if(item == "EMIN_STATE")                            return StorageReason::EMIN_STATE;
        if(item == "EMAX_STATE")                            return StorageReason::EMAX_STATE;
        if(item == "MODEL")                                 return StorageReason::MODEL;
    }
    if constexpr(std::is_same_v<T, StateInit>) {
        if(item == "RANDOM_PRODUCT_STATE")                  return StateInit::RANDOM_PRODUCT_STATE;
        if(item == "RANDOM_ENTANGLED_STATE")                return StateInit::RANDOM_ENTANGLED_STATE;
        if(item == "RANDOMIZE_PREVIOUS_STATE")              return StateInit::RANDOMIZE_PREVIOUS_STATE;
        if(item == "PRODUCT_STATE_ALIGNED")                 return StateInit::PRODUCT_STATE_ALIGNED;
        if(item == "PRODUCT_STATE_NEEL")                    return StateInit::PRODUCT_STATE_NEEL;
    }
    if constexpr(std::is_same_v<T, StateInitType>) {
        if(item ==  "REAL")                                 return StateInitType::REAL;
        if(item ==  "CPLX")                                 return StateInitType::CPLX;
    }
    if constexpr(std::is_same_v<T, PerturbMode>) {
        if(item == "PERCENTAGE")                            return PerturbMode::PERCENTAGE;
        if(item == "ABSOLUTE")                              return PerturbMode::ABSOLUTE;
        if(item == "UNIFORM_RANDOM_PERCENTAGE")             return PerturbMode::UNIFORM_RANDOM_PERCENTAGE;
        if(item == "UNIFORM_RANDOM_ABSOLUTE")               return PerturbMode::UNIFORM_RANDOM_ABSOLUTE;
    }
    if constexpr(std::is_same_v<T, FileCollisionPolicy>) {
        if(item == "RESUME")                                return FileCollisionPolicy::RESUME;
        if(item == "RENAME")                                return FileCollisionPolicy::RENAME;
        if(item == "BACKUP")                                return FileCollisionPolicy::BACKUP;
        if(item == "REPLACE")                               return FileCollisionPolicy::REPLACE;
    }
    if constexpr(std::is_same_v<T, FileResumePolicy>) {
        if(item == "FULL")                                  return FileResumePolicy::FULL;
        if(item == "FAST")                                  return FileResumePolicy::FAST;
    }
    if constexpr(std::is_same_v<T,fdmrg_task>){
        if(item == "INIT_RANDOMIZE_MODEL")                  return fdmrg_task::INIT_RANDOMIZE_MODEL;
        if(item == "INIT_RANDOMIZE_INTO_PRODUCT_STATE")     return fdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE;
        if(item == "INIT_RANDOMIZE_INTO_ENTANGLED_STATE")   return fdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE;
        if(item == "INIT_BOND_DIM_LIMITS")                  return fdmrg_task::INIT_BOND_DIM_LIMITS;
        if(item == "INIT_WRITE_MODEL")                      return fdmrg_task::INIT_WRITE_MODEL;
        if(item == "INIT_CLEAR_STATUS")                     return fdmrg_task::INIT_CLEAR_STATUS;
        if(item == "INIT_DEFAULT")                          return fdmrg_task::INIT_DEFAULT;
        if(item == "FIND_GROUND_STATE")                     return fdmrg_task::FIND_GROUND_STATE;
        if(item == "FIND_HIGHEST_STATE")                    return fdmrg_task::FIND_HIGHEST_STATE;
        if(item == "POST_WRITE_RESULT")                     return fdmrg_task::POST_WRITE_RESULT;
        if(item == "POST_PRINT_RESULT")                     return fdmrg_task::POST_PRINT_RESULT;
        if(item == "POST_DEFAULT")                          return fdmrg_task::POST_DEFAULT;
        if(item == "PROF_RESET")                            return fdmrg_task::PROF_RESET;
    }

    if constexpr(std::is_same_v<T,xdmrg_task>){
        if(item == "INIT_RANDOMIZE_MODEL")                  return xdmrg_task::INIT_RANDOMIZE_MODEL;
        if(item == "INIT_RANDOMIZE_INTO_PRODUCT_STATE")     return xdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE;
        if(item == "INIT_RANDOMIZE_INTO_ENTANGLED_STATE")   return xdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE;
        if(item == "INIT_RANDOMIZE_INTO_STATE_IN_WIN")      return xdmrg_task::INIT_RANDOMIZE_INTO_STATE_IN_WIN;
        if(item == "INIT_RANDOMIZE_FROM_CURRENT_STATE")     return xdmrg_task::INIT_RANDOMIZE_FROM_CURRENT_STATE;
        if(item == "INIT_BOND_DIM_LIMITS")                  return xdmrg_task::INIT_BOND_DIM_LIMITS;
        if(item == "INIT_ENERGY_LIMITS")                    return xdmrg_task::INIT_ENERGY_LIMITS;
        if(item == "INIT_WRITE_MODEL")                      return xdmrg_task::INIT_WRITE_MODEL;
        if(item == "INIT_CLEAR_STATUS")                     return xdmrg_task::INIT_CLEAR_STATUS;
        if(item == "INIT_DEFAULT")                          return xdmrg_task::INIT_DEFAULT;
        if(item == "NEXT_RANDOMIZE_INTO_PRODUCT_STATE")     return xdmrg_task::NEXT_RANDOMIZE_INTO_PRODUCT_STATE;
        if(item == "NEXT_RANDOMIZE_INTO_ENTANGLED_STATE")   return xdmrg_task::NEXT_RANDOMIZE_INTO_ENTANGLED_STATE;
        if(item == "NEXT_RANDOMIZE_PREVIOUS_STATE")         return xdmrg_task::NEXT_RANDOMIZE_PREVIOUS_STATE;
        if(item == "NEXT_RANDOMIZE_INTO_STATE_IN_WIN")      return xdmrg_task::NEXT_RANDOMIZE_INTO_STATE_IN_WIN;
        if(item == "FIND_ENERGY_RANGE")                     return xdmrg_task::FIND_ENERGY_RANGE;
        if(item == "FIND_EXCITED_STATE")                    return xdmrg_task::FIND_EXCITED_STATE;
        if(item == "POST_WRITE_RESULT")                     return xdmrg_task::POST_WRITE_RESULT;
        if(item == "POST_PRINT_RESULT")                     return xdmrg_task::POST_PRINT_RESULT;
        if(item == "POST_DEFAULT")                          return xdmrg_task::POST_DEFAULT;
        if(item == "PROF_RESET")                            return xdmrg_task::PROF_RESET;
    }
    if constexpr(std::is_same_v<T,CopyPolicy>){
        if(item == "FORCE")                                 return CopyPolicy::FORCE;
        if(item == "TRY")                                   return CopyPolicy::TRY;
        if(item == "OFF")                                   return CopyPolicy::OFF;
    }
    if constexpr(std::is_same_v<T,RandomizerMode>){
        if(item == "SHUFFLE")                               return RandomizerMode::SHUFFLE;
        if(item == "SELECT1")                               return RandomizerMode::SELECT1;
        if(item == "ASIS")                                  return RandomizerMode::ASIS;
    }
    if constexpr(std::is_same_v<T,OptType>){
        if(item == "REAL")                                  return OptType::REAL;
        if(item == "CPLX")                                  return OptType::CPLX;
    }
    if constexpr(std::is_same_v<T,OptMode>){
        if(item == "VARIANCE")                              return OptMode::VARIANCE;
        if(item == "OVERLAP")                               return OptMode::OVERLAP;
    }
    if constexpr(std::is_same_v<T,OptSpace>){
        if(item == "SUBSPACE_ONLY")                         return OptSpace::SUBSPACE_ONLY;
        if(item == "SUBSPACE_AND_DIRECT")                   return OptSpace::SUBSPACE_AND_DIRECT;
        if(item == "DIRECT")                                return OptSpace::DIRECT;
    }
    if constexpr(std::is_same_v<T,OptWhen>){
        if(item == "ALWAYS")                                return OptWhen::ALWAYS;
        if(item == "NEVER")                                 return OptWhen::NEVER;
        if(item == "PREV_FAIL")                             return OptWhen::PREV_FAIL;
    }
    if constexpr(std::is_same_v<T,OptMark>){
        if(item == "PASS")                                  return OptMark::PASS;
        if(item == "FAIL")                                  return OptMark::FAIL;
    }
    if constexpr(std::is_same_v<T,OptInit>){
        if(item == "CURRENT_STATE")                         return OptInit::CURRENT_STATE;
        if(item == "LAST_RESULT")                           return OptInit::LAST_RESULT;
    }
    throw std::runtime_error("str2enum given invalid string item: " + std::string(item));
}
/* clang-format on */
