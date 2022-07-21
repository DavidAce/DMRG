#pragma once
#include <array>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

enum class AlgorithmType : int { iDMRG, fDMRG, xDMRG, iTEBD, fLBIT, ANY };
enum class AlgorithmStop : int { SUCCESS, SATURATED, MAX_ITERS, MAX_RESET, RANDOMIZE, NONE };
enum class MultisiteMove { ONE, MID, MAX };
enum class MultisiteWhen { OFF, SATURATED, ALWAYS };
enum class SVDLibrary { EIGEN, LAPACKE, RSVD };
enum class UpdateWhen { NEVER, TRUNCATED, SATURATED, ITERATION };
enum class GateMove { OFF, ON, AUTO };
enum class ModelType { ising_tf_rf, ising_sdual, ising_majorana, lbit };
enum class EdgeStatus { STALE, FRESH };
enum class StorageLevel { NONE, LIGHT, NORMAL, FULL };
enum class StorageReason { SAVEPOINT, CHECKPOINT, FINISHED, PROJ_STATE, INIT_STATE, EMIN_STATE, EMAX_STATE, MODEL, BOND_INCREASE, TRNC_DECREASE, FES };
enum class CopyPolicy { FORCE, TRY, OFF };
enum class ResetReason { INIT, FIND_WINDOW, SATURATED, NEW_STATE, BOND_UPDATE };
enum class NormPolicy { ALWAYS, IFNEEDED }; // Rules of engagement
enum class FileCollisionPolicy {
    RESUME, /*!< If finished -> exit, else resume simulation from the latest "FULL" storage state. Throw if none is found. */
    BACKUP, /*!< Backup the existing file by appending .bak, then start with a new file. */
    RENAME, /*!< Rename the current file by appending .# to avoid collision with existing. */
    REVIVE, /*!< Try RESUME, but do REPLACE on error instead of throwing */
    REPLACE /*!< Just erase/truncate the existing file and start from the beginning. */
};
enum class FileResumePolicy { FULL, FAST };
enum class LogPolicy { NORMAL, QUIET };
enum class RandomizerMode { SHUFFLE, SELECT1, ASIS };
enum class OptType { REAL, CPLX };
enum class OptMode { ENERGY, VARIANCE, OVERLAP, SUBSPACE, SIMPS };
enum class OptSolver { EIGS, BFGS };
enum class OptRitz { LR, SR, SM }; // Smallest Real or Largest Real, i.e. ground state or max state. Use SM for middle of spectrum states.
enum class OptWhen : int {
    NEVER              = 0,
    PREV_FAIL_GRADIENT = 1,
    PREV_FAIL_RESIDUAL = 2,
    PREV_FAIL_OVERLAP  = 4,
    PREV_FAIL_NOCHANGE = 8,
    PREV_FAIL_WORSENED = 16,
    PREV_FAIL_ERROR    = 32,
    ALWAYS             = 64
};
enum class OptEigs { ALWAYS, WHEN_SATURATED }; // When to prefer eigs over bfgs (the default)

enum class OptExit : int {
    SUCCESS       = 0,
    FAIL_GRADIENT = 1,
    FAIL_RESIDUAL = 2,
    FAIL_OVERLAP  = 4,
    FAIL_NOCHANGE = 8,
    FAIL_WORSENED = 16,
    FAIL_ERROR    = 32,
    NONE          = 64,
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
};

enum class fdmrg_task {
    INIT_RANDOMIZE_MODEL,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE,
    INIT_RANDOMIZE_INTO_ENTANGLED_STATE,
    INIT_BOND_LIMITS,
    INIT_TRNC_LIMITS,
    INIT_WRITE_MODEL,
    INIT_CLEAR_STATUS,
    INIT_CLEAR_CONVERGENCE,
    INIT_DEFAULT,
    FIND_GROUND_STATE,
    FIND_HIGHEST_STATE,
    POST_WRITE_RESULT,
    POST_PRINT_RESULT,
    POST_PRINT_TIMERS,
    POST_FES_ANALYSIS,
    POST_DEFAULT,
    TIMER_RESET,
};

enum class flbit_task {
    INIT_RANDOMIZE_MODEL,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE,
    INIT_RANDOMIZE_INTO_ENTANGLED_STATE,
    INIT_BOND_LIMITS,
    INIT_TRNC_LIMITS,
    INIT_WRITE_MODEL,
    INIT_CLEAR_STATUS,
    INIT_CLEAR_CONVERGENCE,
    INIT_DEFAULT,
    INIT_GATES,
    INIT_TIME,
    TRANSFORM_TO_LBIT,
    TRANSFORM_TO_REAL,
    TIME_EVOLVE,
    POST_WRITE_RESULT,
    POST_PRINT_RESULT,
    POST_PRINT_TIMERS,
    POST_FES_ANALYSIS,
    POST_DEFAULT,
    TIMER_RESET,
};

enum class xdmrg_task {
    INIT_RANDOMIZE_MODEL,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE,
    INIT_RANDOMIZE_INTO_ENTANGLED_STATE,
    INIT_RANDOMIZE_INTO_STATE_IN_WIN,
    INIT_RANDOMIZE_FROM_CURRENT_STATE,
    INIT_BOND_LIMITS,
    INIT_TRNC_LIMITS,
    INIT_ENERGY_LIMITS,
    INIT_WRITE_MODEL,
    INIT_CLEAR_STATUS,
    INIT_CLEAR_CONVERGENCE,
    INIT_DEFAULT,
    NEXT_RANDOMIZE_INTO_PRODUCT_STATE,
    NEXT_RANDOMIZE_INTO_ENTANGLED_STATE,
    NEXT_RANDOMIZE_PREVIOUS_STATE,
    NEXT_RANDOMIZE_INTO_STATE_IN_WIN,
    FIND_ENERGY_RANGE,
    FIND_EXCITED_STATE,
    POST_WRITE_RESULT,
    POST_PRINT_RESULT,
    POST_PRINT_TIMERS,
    POST_FES_ANALYSIS,
    POST_DEFAULT,
    TIMER_RESET,
};

template<typename T, bool B = std::is_enum<T>::value>
struct is_scoped_enum : std::false_type {};
template<typename T>
struct is_scoped_enum<T, true> : std::integral_constant<bool, !std::is_convertible<T, typename std::underlying_type<T>::type>::value> {};
template<typename T>
constexpr bool is_scoped_enum_v = is_scoped_enum<T>();

template<typename T>
constexpr bool is_scoped_enum_int_v() {
    if constexpr(is_scoped_enum_v<T>)
        return std::is_same_v<typename std::underlying_type<T>::type, int>;
    else
        return false;
}

template<typename E, typename = std::enable_if_t<is_scoped_enum_int_v<E>()>>
constexpr typename std::underlying_type<E>::type enum2int(E e) noexcept {
    return static_cast<typename std::underlying_type<E>::type>(e);
}

template<typename E, typename = std::enable_if_t<is_scoped_enum_int_v<E>()>>
constexpr inline E operator|(E lhs, E rhs) {
    using T = std::underlying_type_t<E>;
    return static_cast<E>(static_cast<T>(lhs) | static_cast<T>(rhs));
}

template<typename E, typename = std::enable_if_t<is_scoped_enum_int_v<E>()>>
constexpr inline E &operator|=(E &lhs, E rhs) {
    lhs = lhs | rhs;
    return lhs;
}

template<typename E, typename = std::enable_if_t<is_scoped_enum_int_v<E>()>>
inline bool has_flag(E target, E check) {
    return (static_cast<int>(target) & static_cast<int>(check)) == static_cast<int>(check);
}
template<typename E1, typename E2, typename = std::enable_if_t<is_scoped_enum_int_v<E1>() and is_scoped_enum_int_v<E2>()>>
inline bool have_common(E1 lhs, E2 rhs) {
    return (static_cast<int>(lhs) & static_cast<int>(rhs)) != 0;
}

template<typename T>
constexpr std::string_view enum2sv(const T &item) {
    static_assert(std::is_enum_v<T> and "enum2sv<T>: T must be an enum");
    /* clang-format off */
    if constexpr(std::is_same_v<T, AlgorithmType>) switch(item) {
        case AlgorithmType::iDMRG:                                       return "iDMRG";
        case AlgorithmType::fDMRG:                                       return "fDMRG";
        case AlgorithmType::fLBIT:                                       return "fLBIT";
        case AlgorithmType::xDMRG:                                       return "xDMRG";
        case AlgorithmType::iTEBD:                                       return "iTEBD";
        case AlgorithmType::ANY:                                         return "ANY";
    }
    if constexpr(std::is_same_v<T, MultisiteMove>) switch(item){
        case MultisiteMove::ONE :                                       return "ONE";
        case MultisiteMove::MID :                                       return "MID";
        case MultisiteMove::MAX :                                       return "MAX";
    }
    if constexpr(std::is_same_v<T, MultisiteWhen>) switch(item){
        case MultisiteWhen::OFF :                                       return "OFF";
        case MultisiteWhen::SATURATED :                                 return "SATURATED";
        case MultisiteWhen::ALWAYS :                                    return "ALWAYS";
    }
    if constexpr(std::is_same_v<T, OptRitz>) {
        if(item == OptRitz::SR)                                         return "SR";
        if(item == OptRitz::LR)                                         return "LR";
        if(item == OptRitz::SM)                                         return "SM";
    }
    if constexpr(std::is_same_v<T, SVDLibrary>) {
        if(item == SVDLibrary::EIGEN)                                      return "EIGEN";
        if(item == SVDLibrary::LAPACKE)                                    return "LAPACKE";
        if(item == SVDLibrary::RSVD)                                       return "RSVD";
    }
    if constexpr(std::is_same_v<T, UpdateWhen>) {
        if(item == UpdateWhen::NEVER)                                  return "NEVER";
        if(item == UpdateWhen::TRUNCATED)                              return "TRUNCATED";
        if(item == UpdateWhen::SATURATED)                              return "SATURATED";
        if(item == UpdateWhen::ITERATION)                              return "ITERATION";
    }
    if constexpr(std::is_same_v<T, GateMove>) {
        if(item == GateMove::OFF)                                       return "OFF";
        if(item == GateMove::ON)                                        return "ON";
        if(item == GateMove::AUTO)                                      return "AUTO";
    }
    if constexpr(std::is_same_v<T, ModelType>) {
        if(item == ModelType::ising_tf_rf)                              return "ising_tf_rf";
        if(item == ModelType::ising_sdual)                              return "ising_sdual";
        if(item == ModelType::ising_majorana)                           return "ising_majorana";
        if(item == ModelType::lbit)                                     return "lbit";
    }
    if constexpr(std::is_same_v<T, EdgeStatus>) {
        if(item == EdgeStatus::STALE)                                   return "STALE";
        if(item == EdgeStatus::FRESH)                                   return "FRESH";
    }
    if constexpr(std::is_same_v<T, AlgorithmStop>) {
        if(item == AlgorithmStop::SUCCESS)                              return "SUCCESS";
        if(item == AlgorithmStop::SATURATED)                            return "SATURATED";
        if(item == AlgorithmStop::MAX_ITERS)                            return "MAX_ITERS";
        if(item == AlgorithmStop::MAX_RESET)                            return "MAX_RESET";
        if(item == AlgorithmStop::RANDOMIZE)                            return "RANDOMIZE";
        if(item == AlgorithmStop::NONE)                                 return "NONE";
    }
    if constexpr(std::is_same_v<T, ResetReason>) {
        if(item == ResetReason::INIT)                                   return "INIT";
        if(item == ResetReason::FIND_WINDOW)                            return "FIND_WINDOW";
        if(item == ResetReason::SATURATED)                              return "SATURATED";
        if(item == ResetReason::NEW_STATE)                              return "NEW_STATE";
        if(item == ResetReason::BOND_UPDATE)                            return "BOND_UPDATE";
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
        if(item == StorageReason::PROJ_STATE)                           return "PROJ_STATE";
        if(item == StorageReason::INIT_STATE)                           return "INIT_STATE";
        if(item == StorageReason::EMIN_STATE)                           return "EMIN_STATE";
        if(item == StorageReason::EMAX_STATE)                           return "EMAX_STATE";
        if(item == StorageReason::MODEL)                                return "MODEL";
        if(item == StorageReason::BOND_INCREASE)                        return "BOND_INCREASE";
        if(item == StorageReason::TRNC_DECREASE)                        return "TRNC_DECREASE";
        if(item == StorageReason::FES)                                  return "FES";
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
    if constexpr(std::is_same_v<T, FileCollisionPolicy>) {
        if(item == FileCollisionPolicy::RESUME)                         return "RESUME";
        if(item == FileCollisionPolicy::RENAME)                         return "RENAME";
        if(item == FileCollisionPolicy::BACKUP)                         return "BACKUP";
        if(item == FileCollisionPolicy::REVIVE)                         return "REVIVE";
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
        if(item == fdmrg_task::INIT_BOND_LIMITS)                        return "INIT_BOND_LIMITS";
        if(item == fdmrg_task::INIT_TRNC_LIMITS)                        return "INIT_TRNC_LIMITS";
        if(item == fdmrg_task::INIT_WRITE_MODEL)                        return "INIT_WRITE_MODEL";
        if(item == fdmrg_task::INIT_CLEAR_STATUS)                       return "INIT_CLEAR_STATUS";
        if(item == fdmrg_task::INIT_CLEAR_CONVERGENCE)                  return "INIT_CLEAR_CONVERGENCE";
        if(item == fdmrg_task::INIT_DEFAULT)                            return "INIT_DEFAULT";
        if(item == fdmrg_task::FIND_GROUND_STATE)                       return "FIND_GROUND_STATE";
        if(item == fdmrg_task::FIND_HIGHEST_STATE)                      return "FIND_HIGHEST_STATE";
        if(item == fdmrg_task::POST_WRITE_RESULT)                       return "POST_WRITE_RESULT";
        if(item == fdmrg_task::POST_PRINT_RESULT)                       return "POST_PRINT_RESULT";
        if(item == fdmrg_task::POST_DEFAULT)                            return "POST_DEFAULT";
        if(item == fdmrg_task::TIMER_RESET)                             return "TIMER_RESET";

    }
    if constexpr(std::is_same_v<T,flbit_task>){
        if(item == flbit_task::INIT_RANDOMIZE_MODEL)                  return "INIT_RANDOMIZE_MODEL";
        if(item == flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE)     return "INIT_RANDOMIZE_INTO_PRODUCT_STATE";
        if(item == flbit_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE)   return "INIT_RANDOMIZE_INTO_ENTANGLED_STATE";
        if(item == flbit_task::INIT_BOND_LIMITS)                      return "INIT_BOND_LIMITS";
        if(item == flbit_task::INIT_TRNC_LIMITS)                      return "INIT_TRNC_LIMITS";
        if(item == flbit_task::INIT_WRITE_MODEL)                      return "INIT_WRITE_MODEL";
        if(item == flbit_task::INIT_CLEAR_STATUS)                     return "INIT_CLEAR_STATUS";
        if(item == flbit_task::INIT_CLEAR_CONVERGENCE)                return "INIT_CLEAR_CONVERGENCE";
        if(item == flbit_task::INIT_DEFAULT)                          return "INIT_DEFAULT";
        if(item == flbit_task::INIT_GATES)                            return "INIT_GATES";
        if(item == flbit_task::INIT_TIME)                             return "INIT_TIME";
        if(item == flbit_task::TRANSFORM_TO_LBIT)                     return "TRANSFORM_TO_LBIT";
        if(item == flbit_task::TRANSFORM_TO_REAL)                     return "TRANSFORM_TO_REAL";
        if(item == flbit_task::TIME_EVOLVE)                           return "TIME_EVOLVE";
        if(item == flbit_task::POST_WRITE_RESULT)                     return "POST_WRITE_RESULT";
        if(item == flbit_task::POST_PRINT_RESULT)                     return "POST_PRINT_RESULT";
        if(item == flbit_task::POST_PRINT_TIMERS)                     return "POST_PRINT_TIMERS";
        if(item == flbit_task::POST_FES_ANALYSIS)                     return "POST_FES_ANALYSIS";
        if(item == flbit_task::POST_DEFAULT)                          return "POST_DEFAULT";
        if(item == flbit_task::TIMER_RESET)                           return "TIMER_RESET";
    }
    if constexpr(std::is_same_v<T,xdmrg_task>){
        if(item == xdmrg_task::INIT_RANDOMIZE_MODEL)                   return "INIT_RANDOMIZE_MODEL";
        if(item == xdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE)      return "INIT_RANDOMIZE_INTO_PRODUCT_STATE";
        if(item == xdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE)    return "INIT_RANDOMIZE_INTO_ENTANGLED_STATE";
        if(item == xdmrg_task::INIT_RANDOMIZE_INTO_STATE_IN_WIN)       return "INIT_RANDOMIZE_INTO_STATE_IN_WIN";
        if(item == xdmrg_task::INIT_RANDOMIZE_FROM_CURRENT_STATE)      return "INIT_RANDOMIZE_FROM_CURRENT_STATE";
        if(item == xdmrg_task::INIT_BOND_LIMITS)                       return "INIT_BOND_LIMITS";
        if(item == xdmrg_task::INIT_TRNC_LIMITS)                       return "INIT_TRNC_LIMITS";
        if(item == xdmrg_task::INIT_ENERGY_LIMITS)                     return "INIT_ENERGY_LIMITS";
        if(item == xdmrg_task::INIT_WRITE_MODEL)                       return "INIT_WRITE_MODEL";
        if(item == xdmrg_task::INIT_CLEAR_STATUS)                      return "INIT_CLEAR_STATUS";
        if(item == xdmrg_task::INIT_CLEAR_CONVERGENCE)                 return "INIT_CLEAR_CONVERGENCE";
        if(item == xdmrg_task::INIT_DEFAULT)                           return "INIT_DEFAULT";
        if(item == xdmrg_task::NEXT_RANDOMIZE_INTO_PRODUCT_STATE)      return "NEXT_RANDOMIZE_INTO_PRODUCT_STATE";
        if(item == xdmrg_task::NEXT_RANDOMIZE_INTO_ENTANGLED_STATE)    return "NEXT_RANDOMIZE_INTO_ENTANGLED_STATE";
        if(item == xdmrg_task::NEXT_RANDOMIZE_PREVIOUS_STATE)          return "NEXT_RANDOMIZE_PREVIOUS_STATE";
        if(item == xdmrg_task::NEXT_RANDOMIZE_INTO_STATE_IN_WIN)       return "NEXT_RANDOMIZE_INTO_STATE_IN_WIN";
        if(item == xdmrg_task::FIND_ENERGY_RANGE)                      return "FIND_ENERGY_RANGE";
        if(item == xdmrg_task::FIND_EXCITED_STATE)                     return "FIND_EXCITED_STATE";
        if(item == xdmrg_task::POST_WRITE_RESULT)                      return "POST_WRITE_RESULT";
        if(item == xdmrg_task::POST_PRINT_RESULT)                      return "POST_PRINT_RESULT";
        if(item == xdmrg_task::POST_PRINT_TIMERS)                      return "POST_PRINT_TIMERS";
        if(item == xdmrg_task::POST_FES_ANALYSIS)                      return "POST_FES_ANALYSIS";
        if(item == xdmrg_task::POST_DEFAULT)                           return "POST_DEFAULT";
        if(item == xdmrg_task::TIMER_RESET)                            return "TIMER_RESET";
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
        if(item == OptMode::ENERGY  )                                  return "ENERGY";
        if(item == OptMode::VARIANCE)                                  return "VARIANCE";
        if(item == OptMode::OVERLAP)                                   return "OVERLAP";
        if(item == OptMode::SUBSPACE)                                  return "SUBSPACE";
    }
    if constexpr(std::is_same_v<T,OptSolver>){
        if(item == OptSolver::EIGS)                                    return "EIGS";
        if(item == OptSolver::BFGS)                                    return "BFGS";
    }
    if constexpr(std::is_same_v<T,OptWhen>){
        if(item == OptWhen::NEVER)                                     return "NEVER";
        if(item == OptWhen::PREV_FAIL_GRADIENT)                        return "PREV_FAIL_GRADIENT";
        if(item == OptWhen::PREV_FAIL_RESIDUAL)                        return "PREV_FAIL_RESIDUAL";
        if(item == OptWhen::PREV_FAIL_OVERLAP)                         return "PREV_FAIL_OVERLAP";
        if(item == OptWhen::PREV_FAIL_NOCHANGE)                        return "PREV_FAIL_NOCHANGE";
        if(item == OptWhen::PREV_FAIL_WORSENED)                        return "PREV_FAIL_WORSENED";
        if(item == OptWhen::PREV_FAIL_ERROR)                           return "PREV_FAIL_ERROR";
        if(item == OptWhen::ALWAYS)                                    return "ALWAYS";
    }
    if constexpr(std::is_same_v<T,OptEigs>){
        if(item == OptEigs::ALWAYS)                                    return "ALWAYS";
        if(item == OptEigs::WHEN_SATURATED)                            return "WHEN_SATURATED";
    }
    if constexpr(std::is_same_v<T,OptMark>){
        if(item == OptMark::PASS)                                      return "PASS";
        if(item == OptMark::FAIL)                                      return "FAIL";
    }
    if constexpr(std::is_same_v<T,OptInit>){
        if(item == OptInit::CURRENT_STATE)                             return "CURRENT_STATE";
        if(item == OptInit::LAST_RESULT)                               return "LAST_RESULT";
    }
    if constexpr(std::is_same_v<T,OptExit>){
        if(item == OptExit::SUCCESS)                                   return "SUCCESS";
        if(item == OptExit::FAIL_GRADIENT)                             return "FAIL_GRADIENT";
        if(item == OptExit::FAIL_RESIDUAL)                             return "FAIL_RESIDUAL";
        if(item == OptExit::FAIL_OVERLAP)                              return "FAIL_OVERLAP";
        if(item == OptExit::FAIL_NOCHANGE)                             return "FAIL_NOCHANGE";
        if(item == OptExit::FAIL_WORSENED)                             return "FAIL_WORSENED";
        if(item == OptExit::FAIL_ERROR)                                return "FAIL_ERROR";
        if(item == OptExit::NONE)                                      return "NONE";
    }
    throw std::runtime_error("Given invalid enum item");
}

template<typename T>
std::string flag2str(const T &item) {
    static_assert((std::is_same_v<T, OptExit> or std::is_same_v<T,OptWhen>) and
        "flag2str only works with OptWhen or OptExit");

    std::vector<std::string> v;
    if constexpr(std::is_same_v<T,OptExit>){
        if(has_flag(item,OptExit::FAIL_GRADIENT)) v.emplace_back("FAIL_GRADIENT");
        if(has_flag(item,OptExit::FAIL_RESIDUAL)) v.emplace_back("FAIL_RESIDUAL");
        if(has_flag(item,OptExit::FAIL_OVERLAP))  v.emplace_back("FAIL_OVERLAP");
        if(has_flag(item,OptExit::FAIL_NOCHANGE)) v.emplace_back("FAIL_NOCHANGE");
        if(has_flag(item,OptExit::FAIL_WORSENED)) v.emplace_back("FAIL_WORSENED");
        if(has_flag(item,OptExit::FAIL_ERROR))    v.emplace_back("FAIL_ERROR");
        if(has_flag(item,OptExit::NONE))          v.emplace_back("NONE");
        if(v.empty()) return "SUCCESS";
    }
    if constexpr(std::is_same_v<T,OptWhen>){
        if(has_flag(item,OptWhen::PREV_FAIL_GRADIENT)) v.emplace_back("PREV_FAIL_GRADIENT");
        if(has_flag(item,OptWhen::PREV_FAIL_RESIDUAL)) v.emplace_back("PREV_FAIL_RESIDUAL");
        if(has_flag(item,OptWhen::PREV_FAIL_OVERLAP))  v.emplace_back("PREV_FAIL_OVERLAP");
        if(has_flag(item,OptWhen::PREV_FAIL_NOCHANGE)) v.emplace_back("PREV_FAIL_NOCHANGE");
        if(has_flag(item,OptWhen::PREV_FAIL_WORSENED)) v.emplace_back("PREV_FAIL_WORSENED");
        if(has_flag(item,OptWhen::PREV_FAIL_ERROR))    v.emplace_back("PREV_FAIL_ERROR");
        if(has_flag(item,OptWhen::ALWAYS))             v.emplace_back("ALWAYS");
        if(v.empty()) return "NEVER";
    }
    return  std::accumulate(std::begin(v), std::end(v), std::string(),
                            [](const std::string &ss, const std::string &s)
                            {return ss.empty() ? s : ss + "|" + s;});

}


template<class T, class... Ts>
struct is_any : std::disjunction<std::is_same<T, Ts>...> {};
template<class T, class... Ts>
inline constexpr bool is_any_v = is_any<T, Ts...>::value;



template<typename T>
constexpr auto sv2enum(std::string_view item) {
    static_assert(is_any_v<T,
        AlgorithmType,
        AlgorithmStop,
        MultisiteMove,
        MultisiteWhen,
        SVDLibrary,
        UpdateWhen,
        GateMove,
        ModelType,
        EdgeStatus,
        StorageLevel,
        StorageReason,
        CopyPolicy,
        ResetReason,
        NormPolicy,
        FileCollisionPolicy,
        FileResumePolicy,
        LogPolicy,
        RandomizerMode,
        OptType,
        OptMode,
        OptSolver,
        OptRitz,
        OptWhen,
        OptExit,
        OptMark,
        OptInit,
        OptEigs,
        StateInitType,
        StateInit,
        fdmrg_task,
        xdmrg_task,
        flbit_task>);



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
    if constexpr(std::is_same_v<T, MultisiteWhen>) {
        if(item == "OFF")                                   return MultisiteWhen::OFF;
        if(item == "SATURATED")                             return MultisiteWhen::SATURATED;
        if(item == "ALWAYS")                                return MultisiteWhen::ALWAYS;
    }
    if constexpr(std::is_same_v<T, OptRitz>) {
        if(item == "SR")                                    return OptRitz::SR;
        if(item == "LR")                                    return OptRitz::LR;
        if(item == "SM")                                    return OptRitz::SM;
    }
    if constexpr(std::is_same_v<T, SVDLibrary>) {
        if(item == "EIGEN")                                 return SVDLibrary::EIGEN;
        if(item == "LAPACKE")                               return SVDLibrary::LAPACKE;
        if(item == "RSVD")                                  return SVDLibrary::RSVD;
    }
    if constexpr(std::is_same_v<T, UpdateWhen>) {
        if(item == "NEVER")                                 return UpdateWhen::NEVER;
        if(item == "TRUNCATED")                             return UpdateWhen::TRUNCATED;
        if(item == "SATURATED")                             return UpdateWhen::SATURATED;
        if(item == "ITERATION")                             return UpdateWhen::ITERATION;
    }
    if constexpr(std::is_same_v<T, GateMove>) {
        if(item == "OFF")                                   return GateMove::OFF;
        if(item == "ON")                                    return GateMove::ON;
        if(item == "AUTO")                                  return GateMove::AUTO;
    }
    if constexpr(std::is_same_v<T, ModelType>) {
        if(item == "ising_tf_rf")                           return ModelType::ising_tf_rf;
        if(item == "ising_sdual")                           return ModelType::ising_sdual;
        if(item == "ising_majorana")                        return ModelType::ising_majorana;
        if(item == "lbit")                                  return ModelType::lbit;
    }
    if constexpr(std::is_same_v<T, EdgeStatus>) {
        if(item == "STALE")                                 return EdgeStatus::STALE ;
        if(item == "FRESH")                                 return EdgeStatus::FRESH ;
    }
    if constexpr(std::is_same_v<T, AlgorithmStop>) {
        if(item == "SUCCESS")                               return AlgorithmStop::SUCCESS;
        if(item == "SATURATED")                             return AlgorithmStop::SATURATED;
        if(item == "MAX_ITERS")                             return AlgorithmStop::MAX_ITERS;
        if(item == "MAX_RESET")                             return AlgorithmStop::MAX_RESET;
        if(item == "RANDOMIZE")                             return AlgorithmStop::RANDOMIZE;
        if(item == "NONE")                                  return AlgorithmStop::NONE;
    }
    if constexpr(std::is_same_v<T, ResetReason>) {
        if(item == "INIT")                                  return ResetReason::INIT;
        if(item == "FIND_WINDOW")                           return ResetReason::FIND_WINDOW;
        if(item == "SATURATED")                             return ResetReason::SATURATED;
        if(item == "NEW_STATE")                             return ResetReason::NEW_STATE;
        if(item == "BOND_UPDATE")                           return ResetReason::BOND_UPDATE;
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
        if(item == "PROJ_STATE")                            return StorageReason::PROJ_STATE;
        if(item == "INIT_STATE")                            return StorageReason::INIT_STATE;
        if(item == "EMIN_STATE")                            return StorageReason::EMIN_STATE;
        if(item == "EMAX_STATE")                            return StorageReason::EMAX_STATE;
        if(item == "MODEL")                                 return StorageReason::MODEL;
        if(item == "BOND_INCREASE")                         return StorageReason::BOND_INCREASE;
        if(item == "TRNC_DECREASE")                         return StorageReason::TRNC_DECREASE;
        if(item == "FES")                                   return StorageReason::FES;

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
    if constexpr(std::is_same_v<T, FileCollisionPolicy>) {
        if(item == "RESUME")                                return FileCollisionPolicy::RESUME;
        if(item == "RENAME")                                return FileCollisionPolicy::RENAME;
        if(item == "BACKUP")                                return FileCollisionPolicy::BACKUP;
        if(item == "REVIVE")                                return FileCollisionPolicy::REVIVE;
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
        if(item == "INIT_BOND_LIMITS")                      return fdmrg_task::INIT_BOND_LIMITS;
        if(item == "INIT_TRNC_LIMITS")                      return fdmrg_task::INIT_TRNC_LIMITS;
        if(item == "INIT_WRITE_MODEL")                      return fdmrg_task::INIT_WRITE_MODEL;
        if(item == "INIT_CLEAR_STATUS")                     return fdmrg_task::INIT_CLEAR_STATUS;
        if(item == "INIT_CLEAR_CONVERGENCE")                return fdmrg_task::INIT_CLEAR_CONVERGENCE;
        if(item == "INIT_DEFAULT")                          return fdmrg_task::INIT_DEFAULT;
        if(item == "FIND_GROUND_STATE")                     return fdmrg_task::FIND_GROUND_STATE;
        if(item == "FIND_HIGHEST_STATE")                    return fdmrg_task::FIND_HIGHEST_STATE;
        if(item == "POST_WRITE_RESULT")                     return fdmrg_task::POST_WRITE_RESULT;
        if(item == "POST_PRINT_RESULT")                     return fdmrg_task::POST_PRINT_RESULT;
        if(item == "POST_PRINT_TIMERS")                     return fdmrg_task::POST_PRINT_TIMERS;
        if(item == "POST_FES_ANALYSIS")                     return fdmrg_task::POST_FES_ANALYSIS;
        if(item == "POST_DEFAULT")                          return fdmrg_task::POST_DEFAULT;
        if(item == "TIMER_RESET")                           return fdmrg_task::TIMER_RESET;
    }

    if constexpr(std::is_same_v<T,flbit_task>){
        if(item == "INIT_RANDOMIZE_MODEL")                  return flbit_task::INIT_RANDOMIZE_MODEL;
        if(item == "INIT_RANDOMIZE_INTO_PRODUCT_STATE")     return flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE;
        if(item == "INIT_RANDOMIZE_INTO_ENTANGLED_STATE")   return flbit_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE;
        if(item == "INIT_BOND_LIMITS")                      return flbit_task::INIT_BOND_LIMITS;
        if(item == "INIT_TRNC_LIMITS")                      return flbit_task::INIT_TRNC_LIMITS;
        if(item == "INIT_WRITE_MODEL")                      return flbit_task::INIT_WRITE_MODEL;
        if(item == "INIT_CLEAR_STATUS")                     return flbit_task::INIT_CLEAR_STATUS;
        if(item == "INIT_CLEAR_CONVERGENCE")                return flbit_task::INIT_CLEAR_CONVERGENCE;
        if(item == "INIT_DEFAULT")                          return flbit_task::INIT_DEFAULT;
        if(item == "INIT_GATES")                            return flbit_task::INIT_GATES;
        if(item == "INIT_TIME")                             return flbit_task::INIT_TIME;
        if(item == "TRANSFORM_TO_LBIT")                     return flbit_task::TRANSFORM_TO_LBIT;
        if(item == "TRANSFORM_TO_REAL")                     return flbit_task::TRANSFORM_TO_REAL;
        if(item == "TIME_EVOLVE")                           return flbit_task::TIME_EVOLVE;
        if(item == "POST_WRITE_RESULT")                     return flbit_task::POST_WRITE_RESULT;
        if(item == "POST_PRINT_RESULT")                     return flbit_task::POST_PRINT_RESULT;
        if(item == "POST_PRINT_TIMERS")                     return flbit_task::POST_PRINT_TIMERS;
        if(item == "POST_FES_ANALYSIS")                     return flbit_task::POST_FES_ANALYSIS;
        if(item == "POST_DEFAULT")                          return flbit_task::POST_DEFAULT;
        if(item == "TIMER_RESET")                           return flbit_task::TIMER_RESET;
    }
    if constexpr(std::is_same_v<T,xdmrg_task>){
        if(item == "INIT_RANDOMIZE_MODEL")                  return xdmrg_task::INIT_RANDOMIZE_MODEL;
        if(item == "INIT_RANDOMIZE_INTO_PRODUCT_STATE")     return xdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE;
        if(item == "INIT_RANDOMIZE_INTO_ENTANGLED_STATE")   return xdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE;
        if(item == "INIT_RANDOMIZE_INTO_STATE_IN_WIN")      return xdmrg_task::INIT_RANDOMIZE_INTO_STATE_IN_WIN;
        if(item == "INIT_RANDOMIZE_FROM_CURRENT_STATE")     return xdmrg_task::INIT_RANDOMIZE_FROM_CURRENT_STATE;
        if(item == "INIT_BOND_LIMITS")                      return xdmrg_task::INIT_BOND_LIMITS;
        if(item == "INIT_TRNC_LIMITS")                      return xdmrg_task::INIT_TRNC_LIMITS;
        if(item == "INIT_ENERGY_LIMITS")                    return xdmrg_task::INIT_ENERGY_LIMITS;
        if(item == "INIT_WRITE_MODEL")                      return xdmrg_task::INIT_WRITE_MODEL;
        if(item == "INIT_CLEAR_STATUS")                     return xdmrg_task::INIT_CLEAR_STATUS;
        if(item == "INIT_CLEAR_CONVERGENCE")                return xdmrg_task::INIT_CLEAR_CONVERGENCE;
        if(item == "INIT_DEFAULT")                          return xdmrg_task::INIT_DEFAULT;
        if(item == "NEXT_RANDOMIZE_INTO_PRODUCT_STATE")     return xdmrg_task::NEXT_RANDOMIZE_INTO_PRODUCT_STATE;
        if(item == "NEXT_RANDOMIZE_INTO_ENTANGLED_STATE")   return xdmrg_task::NEXT_RANDOMIZE_INTO_ENTANGLED_STATE;
        if(item == "NEXT_RANDOMIZE_PREVIOUS_STATE")         return xdmrg_task::NEXT_RANDOMIZE_PREVIOUS_STATE;
        if(item == "NEXT_RANDOMIZE_INTO_STATE_IN_WIN")      return xdmrg_task::NEXT_RANDOMIZE_INTO_STATE_IN_WIN;
        if(item == "FIND_ENERGY_RANGE")                     return xdmrg_task::FIND_ENERGY_RANGE;
        if(item == "FIND_EXCITED_STATE")                    return xdmrg_task::FIND_EXCITED_STATE;
        if(item == "POST_WRITE_RESULT")                     return xdmrg_task::POST_WRITE_RESULT;
        if(item == "POST_PRINT_RESULT")                     return xdmrg_task::POST_PRINT_RESULT;
        if(item == "POST_PRINT_TIMERS")                     return xdmrg_task::POST_PRINT_TIMERS;
        if(item == "POST_FES_ANALYSIS")                     return xdmrg_task::POST_FES_ANALYSIS;
        if(item == "POST_DEFAULT")                          return xdmrg_task::POST_DEFAULT;
        if(item == "TIMER_RESET")                           return xdmrg_task::TIMER_RESET;
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
        if(item == "ENERGY")                                return OptMode::ENERGY;
        if(item == "VARIANCE")                              return OptMode::VARIANCE;
        if(item == "OVERLAP")                               return OptMode::OVERLAP;
        if(item == "SUBSPACE")                              return OptMode::SUBSPACE;
    }
    if constexpr(std::is_same_v<T,OptSolver>){
        if(item == "EIGS")                                  return OptSolver::EIGS;
        if(item == "BFGS")                                  return OptSolver::BFGS;
    }
    if constexpr(std::is_same_v<T,OptWhen>){
        if(item == "NEVER")                                 return OptWhen::NEVER;
        if(item == "PREV_FAIL_GRADIENT")                    return OptWhen::PREV_FAIL_GRADIENT;
        if(item == "PREV_FAIL_RESIDUAL")                    return OptWhen::PREV_FAIL_RESIDUAL;
        if(item == "PREV_FAIL_OVERLAP")                     return OptWhen::PREV_FAIL_OVERLAP;
        if(item == "PREV_FAIL_NOCHANGE")                    return OptWhen::PREV_FAIL_NOCHANGE;
        if(item == "PREV_FAIL_WORSENED")                    return OptWhen::PREV_FAIL_WORSENED;
        if(item == "PREV_FAIL_ERROR")                       return OptWhen::PREV_FAIL_ERROR;
        if(item == "ALWAYS")                                return OptWhen::ALWAYS;
    }
    if constexpr(std::is_same_v<T,OptEigs>){
        if(item == "ALWAYS")                                return OptEigs::ALWAYS;
        if(item == "WHEN_SATURATED")                        return OptEigs::WHEN_SATURATED;
    }
    if constexpr(std::is_same_v<T,OptMark>){
        if(item == "PASS")                                  return OptMark::PASS;
        if(item == "FAIL")                                  return OptMark::FAIL;
    }
    if constexpr(std::is_same_v<T,OptInit>){
        if(item == "CURRENT_STATE")                         return OptInit::CURRENT_STATE;
        if(item == "LAST_RESULT")                           return OptInit::LAST_RESULT;
    }
    if constexpr(std::is_same_v<T,OptExit>){
        if(item == "SUCCESS")                               return OptExit::SUCCESS;
        if(item == "FAIL_GRADIENT")                         return OptExit::FAIL_GRADIENT;
        if(item == "FAIL_RESIDUAL")                         return OptExit::FAIL_RESIDUAL;
        if(item == "FAIL_OVERLAP")                          return OptExit::FAIL_OVERLAP;
        if(item == "FAIL_NOCHANGE")                         return OptExit::FAIL_NOCHANGE;
        if(item == "FAIL_WORSENED")                         return OptExit::FAIL_WORSENED;
        if(item == "FAIL_ERROR")                            return OptExit::FAIL_ERROR;
        if(item == "NONE")                                  return OptExit::NONE;
    }
    throw std::runtime_error("sv2enum given invalid string item: " + std::string(item));
}
/* clang-format on */

template<typename T, auto num>
using enumarray_t = std::array<std::pair<std::string, T>, num>;

template<typename T, typename... Args>
constexpr auto mapStr2Enum(Args... names) {
    constexpr auto num     = sizeof...(names);
    auto           pairgen = [](const std::string &name) -> std::pair<std::string, T> {
        return {name, sv2enum<T>(name)};
    };
    return enumarray_t<T, num>{pairgen(names)...};
}

template<typename T, typename... Args>
constexpr auto mapEnum2Str(Args... enums) {
    constexpr auto num     = sizeof...(enums);
    auto           pairgen = [](const T &e) -> std::pair<std::string, T> {
        return {std::string(enum2sv(e)), e};
    };
    return enumarray_t<T, num>{pairgen(enums)...};
}

inline auto ModelType_s2e = mapEnum2Str<ModelType>(ModelType::ising_tf_rf, ModelType::ising_sdual, ModelType::ising_majorana, ModelType::lbit);