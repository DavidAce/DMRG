#pragma once
#include <array>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

enum class AlgorithmType : int { iDMRG, fDMRG, xDMRG, iTEBD, fLBIT, ANY };
enum class AlgorithmStop : int { SUCCESS, SATURATED, MAX_ITERS, NONE };
enum class MeanType { ARITHMETIC, GEOMETRIC };
enum class MultisiteMove { ONE, MID, MAX };
enum class MultisiteWhen { NEVER, STUCK, SATURATED, ALWAYS };
enum class MultisiteGrow {
    OFF, /*!< Use the maximum number of sites immediately when conditions allow */
    ON   /*!< Slowly ramp up to the maximum number of sites as the number of stuck/saturated iterations increase */
};
enum class SVDLibrary { EIGEN, LAPACKE, RSVD };
enum class UpdateWhen { NEVER, TRUNCATED, STUCK, SATURATED, ITERATION };
enum class GateMove { OFF, ON, AUTO };
enum class GateOp { NONE, CNJ, ADJ, TRN };
enum class CircuitOp { NONE, ADJ, TRN };
enum class ModelType { ising_tf_rf, ising_sdual, ising_majorana, lbit, xxz };
enum class EdgeStatus { STALE, FRESH };
enum class TimeScale { LINSPACED, LOGSPACED };

/*! \brief The type of Hermitian matrix M_i used in the 2-site gates u_i of the random unitary circuit for our l-bit model: u_i = exp(-i f w_i M_i)
 */
enum class LbitCircuitGateMatrixKind : int {
    MATRIX_V1, /*!< Below, θᵢ, RE(c), and IM(c) are gaussian N(0,1) and λ is a constant
                * \verbatim
                * M_i = θ₀/4 (1 + σz[i] + σz[i+1] + λ σz[i] σz[i+1])
                *     + θ₁/4 (1 + σz[i] - σz[i+1] - λ σz[i] σz[i+1])
                *     + θ₂/4 (1 - σz[i] + σz[i+1] - λ σz[i] σz[i+1])
                *     + θ₃/4 (1 - σz[i] - σz[i+1] + λ σz[i] σz[i+1])
                *     + c σ+[i]σ-[i+1] + c^* σ-[i]σ+[i+1]
                * \endverbatim */
    MATRIX_V2, /*!< Below, θᵢ, RE(c), and IM(c) are gaussian N(0,1) and λ is a constant
                * \verbatim
                * M_i = θ₀/2 σz[i]
                *     + θ₁/2 σz[i+1]
                *     + θ₂/2 λ σz[i]σz[i+1]
                *     + c σ+[i]σ-[i+1] + c^* σ-[i]σ+[i+1]
                * \endverbatim */
    MATRIX_V3, /*!< Below, θᵢ, RE(c), and IM(c) are gaussian N(0,1) and λ is a constant
                * \verbatim
                * M_i = θ₀/2 σz[i]
                *     + θ₁/2 σz[i+1]
                *     + λ σz[i]σz[i+1]
                *     + c σ+[i]σ-[i+1] + c^* σ-[i]σ+[i+1]
                * \endverbatim */
};

/*! The type of weights w_i used in each 2-site gate u_i of the random unitary circuit for our l-bit model: u_i = exp(-i f w_i M_i) */
enum class LbitCircuitGateWeightKind : int {
    IDENTITY, /*!< w_i = 1 (i.e. disables weighs) */
    EXPDECAY, /*!< w_i = exp(-2|h[i] - h[i+1]|), where h[i] are on-site fields of the l-bit Hamiltonian */
};

/*! The reason that we are invoking a storage call */
enum class StorageEvent : int {
    NONE        = 0,    /*!< No event */
    INIT        = 1,    /*!< the initial state was defined */
    MODEL       = 2,    /*!< the model was defined */
    EMIN        = 4,    /*!< the ground state was found (e.g. before xDMRG) */
    EMAX        = 8,    /*!< the highest energy eigenstate was found (e.g. before xDMRG) */
    PROJECTION  = 16,   /*!< a projection to a spin parity sector was made */
    BOND_UPDATE = 32,   /*!< the bond dimension limit was updated */
    TRNC_UPDATE = 64,   /*!< the truncation error threshold for SVD was updated */
    RBDS_STEP   = 128,  /*!< reverse bond dimension scaling step was made */
    RTES_STEP   = 256,  /*!< reverse truncation error scaling step was made */
    ITERATION   = 512,  /*!< an iteration finished */
    FINISHED    = 1024, /*!< a simulation has finished */
};

/*! Determines in which cases we calculate and store to file
 *  Note that many flags can be set simultaneously
 */
enum class StoragePolicy : int {
    NONE    = 0,    /*!< Never store */
    INIT    = 1,    /*!< Store only once during initialization e.g. model (usually in preprocessing) */
    ITER    = 2,    /*!< Store after every iteration */
    EMIN    = 4,    /*!< Store after finding the ground state (e.g. before xDMRG) */
    EMAX    = 8,    /*!< Store after finding the highest energy eigenstate (e.g. before xDMRG) */
    PROJ    = 16,   /*!< Store after projections */
    BOND    = 32,   /*!< Store after bond updates */
    TRNC    = 64,   /*!< Store after truncation error limit updates */
    FAILURE = 128,  /*!< Store only if the simulation did not succeed (usually for debugging) */
    SUCCESS = 256,  /*!< Store only if the simulation succeeded */
    FINISH  = 512,  /*!< Store when the simulation has finished (regardless of failure or success) */
    ALWAYS  = 1024, /*!< Store every chance you get */
    REPLACE = 2048, /*!< Keep only the last event (i.e. replace previous events when possible) */
    allow_bitops
};

enum class CopyPolicy { FORCE, TRY, OFF };
enum class ResetReason { INIT, FIND_WINDOW, SATURATED, NEW_STATE, BOND_UPDATE };
enum class EnvExpandMode { ENE, VAR };
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
enum class OptFunc {
    /*! Choose the objective function for the optimization.
     * Below, H and |psi> are the effective Hamiltonian and state corresponding to
     * a one or more lattice sites.
     */
    ENERGY,   /*!< Keep an energy eigenstate H|psi> = E|psi> (ground state or excited state)   */
    VARIANCE, /*!< Keep the smallest eigenstate of (H-E)^2|psi> = Var(H) |psi>, where E is a target energy */
    OVERLAP   /*!< Among all eigenstates of H, keep the state which maximizes the overlap with the current state |<psi'|psi>| -> 1 */
};

/*! Choose the algorithm for the optimization.
 *
 * Note that H and |psi> are the effective Hamiltonian and state corresponding to
 * a one or more lattice sites.
 */
enum class OptAlgo {
    DIRECT,   /*!< Apply the eigensolver directly on H|psi> or (H-E)^2|psi> */
    SUBSPACE, /*!< Find a state |psi> = λ_0|psi_0> + λ_1|psi_1>... by finding the ground state of of (H-E)^2_ij = <psi_i|(H-E)^2|psi_j> , where psi_i,psi_j are
                 eigenstates of H which span the current state */
    SHIFTINV, /*!< Find a mid-spectrum eigenstate of H using shift-invert (only compatible with OptFunc::ENERGY).  */
    MPSEIGS   /*!< Apply all MPO's onto all MPS (all L sites) in the matrix-vector product of a matrix-free eigenvalue solver  */
};            //
enum class OptSolver {
    EIG, /*!< Apply the eigensolver directly on H|psi> or (H-E)^2|psi> */
    EIGS /*!< Apply the eigensolver directly on H|psi> or (H-E)^2|psi> */
};

/*! Choose the target energy eigenpair */
enum class OptRitz {
    NONE, /*!< No eigenpair is targeted (e.g. time evolution) */
    LR, /*!< Largest Real eigenvalue */
    SR, /*!< Smallest Real eigenvalue */
    SM, /*!< Smallest magnitude eigenvalue. MPO² Energy shift == 0. Use this to find an eigenstate with energy closest to 0) */
    IS, /*!< Initial State energy. MPO² Energy shift == Initial state energy. Targets an eigenstate with energy near that of the initial state */
    TE  /*!< Target Energy density in normalized units [0,1]. MPO² Energy shift == settings::xdmrg::energy_density_target * (EMIN+EMAX) + EMIN.  */
};

enum class OptWhen : int {
    NEVER              = 0,
    PREV_FAIL_GRADIENT = 1,
    PREV_FAIL_RESIDUAL = 2,
    PREV_FAIL_OVERLAP  = 4,
    PREV_FAIL_NOCHANGE = 8,
    PREV_FAIL_WORSENED = 16,
    PREV_FAIL_ERROR    = 32,
    ALWAYS             = 64,
    allow_bitops
};

enum class OptExit : int {
    SUCCESS       = 0,
    FAIL_GRADIENT = 1,
    FAIL_RESIDUAL = 2,
    FAIL_OVERLAP  = 4,
    FAIL_NOCHANGE = 8,
    FAIL_WORSENED = 16,
    FAIL_ERROR    = 32,
    NONE          = 64,
    allow_bitops
};
enum class OptMark {
    PASS,
    FAIL,
};
enum class StateInitType { REAL, CPLX };
enum class StateInit {
    RANDOM_PRODUCT_STATE,
    RANDOM_ENTANGLED_STATE,
    RANDOMIZE_PREVIOUS_STATE,
    MIDCHAIN_SINGLET_NEEL_STATE,
    PRODUCT_STATE_ALIGNED,
    PRODUCT_STATE_DOMAIN_WALL,
    PRODUCT_STATE_NEEL,
    PRODUCT_STATE_NEEL_SHUFFLED,
    PRODUCT_STATE_NEEL_DISLOCATED,
    PRODUCT_STATE_PATTERN,
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
    POST_RBDS_ANALYSIS,
    POST_RTES_ANALYSIS,
    POST_DEFAULT,
    TIMER_RESET,
};

enum class flbit_task {
    INIT_RANDOMIZE_MODEL,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_SHUFFLED,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_DISLOCATED,
    INIT_RANDOMIZE_INTO_ENTANGLED_STATE,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE_PATTERN,
    INIT_RANDOMIZE_INTO_MIDCHAIN_SINGLET_NEEL_STATE,
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
    POST_DEFAULT,
    TIMER_RESET,
};

enum class xdmrg_task {
    INIT_RANDOMIZE_MODEL,
    INIT_RANDOMIZE_INTO_PRODUCT_STATE,
    INIT_RANDOMIZE_INTO_ENTANGLED_STATE,
    INIT_RANDOMIZE_FROM_CURRENT_STATE,
    INIT_BOND_LIMITS,
    INIT_TRNC_LIMITS,
    INIT_ENERGY_TARGET,
    INIT_WRITE_MODEL,
    INIT_CLEAR_STATUS,
    INIT_CLEAR_CONVERGENCE,
    INIT_DEFAULT,
    FIND_ENERGY_RANGE,
    FIND_EXCITED_STATE,
    POST_WRITE_RESULT,
    POST_PRINT_RESULT,
    POST_PRINT_TIMERS,
    POST_RBDS_ANALYSIS,
    POST_RTES_ANALYSIS,
    POST_DEFAULT,
    TIMER_RESET,
};

template<typename T, typename = std::void_t<>>
struct enum_is_bitflag : std::false_type {};
template<typename T>
struct enum_is_bitflag<T, std::void_t<decltype(T::allow_bitops)>> : public std::true_type {};
template<typename T>
constexpr bool enum_is_bitflag_v = enum_is_bitflag<T>();

template<typename E>
constexpr auto operator|(E lhs, E rhs) noexcept -> decltype(E::allow_bitops) {
    using U = std::underlying_type_t<E>;
    return static_cast<E>(static_cast<U>(lhs) | static_cast<U>(rhs));
}
template<typename E>
constexpr auto operator&(E lhs, E rhs) noexcept -> decltype(E::allow_bitops) {
    using U = std::underlying_type_t<E>;
    return static_cast<E>(static_cast<U>(lhs) & static_cast<U>(rhs));
}
template<typename E>
constexpr auto operator|=(E &lhs, E rhs) noexcept -> decltype(E::allow_bitops) {
    lhs = lhs | rhs;
    return lhs;
}

template<typename E>
inline bool has_flag(E target, E check) noexcept {
    using U = std::underlying_type_t<E>;
    return (static_cast<U>(target) & static_cast<U>(check)) == static_cast<U>(check);
}

template<typename E1, typename E2>
inline bool have_common(E1 lhs, E2 rhs) noexcept {
    using U1 = std::underlying_type_t<E1>;
    using U2 = std::underlying_type_t<E2>;
    static_assert(std::is_same_v<U1, U2>);
    return (static_cast<U1>(lhs) & static_cast<U2>(rhs)) != 0;
}

namespace enum_sfinae {
    template<class T, class... Ts>
    struct is_any : std::disjunction<std::is_same<T, Ts>...> {};
    template<class T, class... Ts>
    inline constexpr bool is_any_v = is_any<T, Ts...>::value;
}

template<typename T>
std::vector<std::string_view> enum2sv(const std::vector<T> &items) noexcept {
    auto v = std::vector<std::string_view>();
    v.reserve(items.size());
    for(const auto &item : items) v.emplace_back(enum2sv(item));
    return v;
}

template<typename T>
std::string flag2str(const T &item) noexcept {
    static_assert(std::is_same_v<T, OptExit> or //
                  std::is_same_v<T, OptWhen> or //
                  std::is_same_v<T, StoragePolicy>);

    std::vector<std::string> v;
    if constexpr(std::is_same_v<T, OptExit>) {
        if(has_flag(item, OptExit::FAIL_GRADIENT)) v.emplace_back("FAIL_GRADIENT");
        if(has_flag(item, OptExit::FAIL_RESIDUAL)) v.emplace_back("FAIL_RESIDUAL");
        if(has_flag(item, OptExit::FAIL_OVERLAP)) v.emplace_back("FAIL_OVERLAP");
        if(has_flag(item, OptExit::FAIL_NOCHANGE)) v.emplace_back("FAIL_NOCHANGE");
        if(has_flag(item, OptExit::FAIL_WORSENED)) v.emplace_back("FAIL_WORSENED");
        if(has_flag(item, OptExit::FAIL_ERROR)) v.emplace_back("FAIL_ERROR");
        if(has_flag(item, OptExit::NONE)) v.emplace_back("NONE");
        if(v.empty()) return "SUCCESS";
    }
    if constexpr(std::is_same_v<T, OptWhen>) {
        if(has_flag(item, OptWhen::PREV_FAIL_GRADIENT)) v.emplace_back("PREV_FAIL_GRADIENT");
        if(has_flag(item, OptWhen::PREV_FAIL_RESIDUAL)) v.emplace_back("PREV_FAIL_RESIDUAL");
        if(has_flag(item, OptWhen::PREV_FAIL_OVERLAP)) v.emplace_back("PREV_FAIL_OVERLAP");
        if(has_flag(item, OptWhen::PREV_FAIL_NOCHANGE)) v.emplace_back("PREV_FAIL_NOCHANGE");
        if(has_flag(item, OptWhen::PREV_FAIL_WORSENED)) v.emplace_back("PREV_FAIL_WORSENED");
        if(has_flag(item, OptWhen::PREV_FAIL_ERROR)) v.emplace_back("PREV_FAIL_ERROR");
        if(has_flag(item, OptWhen::ALWAYS)) v.emplace_back("ALWAYS");
        if(v.empty()) return "NEVER";
    }
    if constexpr(std::is_same_v<T, StoragePolicy>) {
        if(has_flag(item, StoragePolicy::INIT)) v.emplace_back("INIT");
        if(has_flag(item, StoragePolicy::ITER)) v.emplace_back("ITER");
        if(has_flag(item, StoragePolicy::EMIN)) v.emplace_back("EMIN");
        if(has_flag(item, StoragePolicy::EMAX)) v.emplace_back("EMAX");
        if(has_flag(item, StoragePolicy::PROJ)) v.emplace_back("PROJ");
        if(has_flag(item, StoragePolicy::BOND)) v.emplace_back("BOND");
        if(has_flag(item, StoragePolicy::TRNC)) v.emplace_back("TRNC");
        if(has_flag(item, StoragePolicy::FAILURE)) v.emplace_back("FAILURE");
        if(has_flag(item, StoragePolicy::SUCCESS)) v.emplace_back("SUCCESS");
        if(has_flag(item, StoragePolicy::FINISH)) v.emplace_back("FINISH");
        if(has_flag(item, StoragePolicy::ALWAYS)) v.emplace_back("ALWAYS");
        if(has_flag(item, StoragePolicy::REPLACE)) v.emplace_back("REPLACE");
        if(v.empty()) return "NONE";
    }
    return std::accumulate(std::begin(v), std::end(v), std::string(),
                           [](const std::string &ss, const std::string &s) { return ss.empty() ? s : ss + "|" + s; });
}

template<typename T>
constexpr std::string_view enum2sv(const T item) noexcept {
    /* clang-format off */
    static_assert(std::is_enum_v<T> and "enum2sv<T>: T must be an enum");
    static_assert(enum_sfinae::is_any_v<T,
        AlgorithmType,
        AlgorithmStop,
        MeanType,
        MultisiteMove,
        MultisiteWhen,
        MultisiteGrow,
        SVDLibrary,
        UpdateWhen,
        GateMove,
        GateOp,
        CircuitOp,
        LbitCircuitGateMatrixKind,
        LbitCircuitGateWeightKind,
        ModelType,
        EdgeStatus,
        TimeScale,
        StorageEvent,
        StoragePolicy,
        CopyPolicy,
        ResetReason,
        NormPolicy,
        FileCollisionPolicy,
        FileResumePolicy,
        LogPolicy,
        RandomizerMode,
        OptType,
        OptFunc,
        OptAlgo,
        OptSolver,
        OptRitz,
        OptWhen,
        OptExit,
        OptMark,
        StateInitType,
        StateInit,
        fdmrg_task,
        xdmrg_task,
        flbit_task>);


    if constexpr(std::is_same_v<T, AlgorithmType>) switch(item) {
        case AlgorithmType::iDMRG:                                      return "iDMRG";
        case AlgorithmType::fDMRG:                                      return "fDMRG";
        case AlgorithmType::fLBIT:                                      return "fLBIT";
        case AlgorithmType::xDMRG:                                      return "xDMRG";
        case AlgorithmType::iTEBD:                                      return "iTEBD";
        case AlgorithmType::ANY:                                        return "ANY";
        default: return "AlgorithmType::UNDEFINED";
    }
    if constexpr(std::is_same_v<T, MultisiteMove>) switch(item){
        case MultisiteMove::ONE :                                       return "ONE";
        case MultisiteMove::MID :                                       return "MID";
        case MultisiteMove::MAX :                                       return "MAX";
        default: return "MultisiteMove::UNDEFINED";
    }
    if constexpr(std::is_same_v<T, MultisiteWhen>) switch(item){
        case MultisiteWhen::NEVER :                                     return "NEVER";
        case MultisiteWhen::STUCK :                                     return "STUCK";
        case MultisiteWhen::SATURATED :                                 return "SATURATED";
        case MultisiteWhen::ALWAYS :                                    return "ALWAYS";
        default: return "MultisiteWhen::UNDEFINED";
    }
    if constexpr(std::is_same_v<T, MultisiteGrow>) switch(item){
        case MultisiteGrow::OFF :                                       return "NEVER";
        case MultisiteGrow::ON :                                        return "ON";
        default: return "MultisiteGrow::UNDEFINED";
    }
    if constexpr(std::is_same_v<T, OptRitz>) {
        if(item == OptRitz::NONE)                                       return "NONE";
        if(item == OptRitz::SR)                                         return "SR";
        if(item == OptRitz::LR)                                         return "LR";
        if(item == OptRitz::SM)                                         return "SM";
        if(item == OptRitz::IS)                                         return "IS";
        if(item == OptRitz::TE)                                         return "TE";
    }
    if constexpr(std::is_same_v<T, SVDLibrary>) {
        if(item == SVDLibrary::EIGEN)                                   return "EIGEN";
        if(item == SVDLibrary::LAPACKE)                                 return "LAPACKE";
        if(item == SVDLibrary::RSVD)                                    return "RSVD";
    }
    if constexpr(std::is_same_v<T, UpdateWhen>) {
        if(item == UpdateWhen::NEVER)                                   return "NEVER";
        if(item == UpdateWhen::TRUNCATED)                               return "TRUNCATED";
        if(item == UpdateWhen::STUCK)                                   return "STUCK";
        if(item == UpdateWhen::SATURATED)                               return "SATURATED";
        if(item == UpdateWhen::ITERATION)                               return "ITERATION";
    }
    if constexpr(std::is_same_v<T, GateMove>) {
        if(item == GateMove::OFF)                                       return "OFF";
        if(item == GateMove::ON)                                        return "ON";
        if(item == GateMove::AUTO)                                      return "AUTO";
    }
    if constexpr(std::is_same_v<T, GateOp>) {
        if(item == GateOp::NONE)                                        return "NONE";
        if(item == GateOp::CNJ)                                         return "CNJ";
        if(item == GateOp::ADJ)                                         return "ADJ";
        if(item == GateOp::TRN)                                         return "TRN";
    }
    if constexpr(std::is_same_v<T, CircuitOp>) {
        if(item == CircuitOp::NONE)                                     return "NONE";
        if(item == CircuitOp::ADJ)                                      return "ADJ";
        if(item == CircuitOp::TRN)                                      return "TRN";
    }
    if constexpr(std::is_same_v<T, LbitCircuitGateMatrixKind>) {
        if(item == LbitCircuitGateMatrixKind::MATRIX_V1)                    return "MATRIX_V1";
        if(item == LbitCircuitGateMatrixKind::MATRIX_V2)                    return "MATRIX_V2";
        if(item == LbitCircuitGateMatrixKind::MATRIX_V3)                    return "MATRIX_V3";
    }
    if constexpr(std::is_same_v<T, LbitCircuitGateWeightKind>) {
        if(item == LbitCircuitGateWeightKind::IDENTITY)                     return "IDENTITY";
        if(item == LbitCircuitGateWeightKind::EXPDECAY)                     return "EXPDECAY";
    }
    if constexpr(std::is_same_v<T, ModelType>) {
        if(item == ModelType::ising_tf_rf)                              return "ising_tf_rf";
        if(item == ModelType::ising_sdual)                              return "ising_sdual";
        if(item == ModelType::ising_majorana)                           return "ising_majorana";
        if(item == ModelType::lbit)                                     return "lbit";
        if(item == ModelType::xxz)                                      return "xxz";
    }
    if constexpr(std::is_same_v<T, EdgeStatus>) {
        if(item == EdgeStatus::STALE)                                   return "STALE";
        if(item == EdgeStatus::FRESH)                                   return "FRESH";
    }
    if constexpr(std::is_same_v<T, TimeScale>) {
        if(item == TimeScale::LINSPACED)                                return "LINSPACED";
        if(item == TimeScale::LOGSPACED)                                return "LOGSPACED";
    }
    if constexpr(std::is_same_v<T, AlgorithmStop>) {
        if(item == AlgorithmStop::SUCCESS)                              return "SUCCESS";
        if(item == AlgorithmStop::SATURATED)                            return "SATURATED";
        if(item == AlgorithmStop::MAX_ITERS)                            return "MAX_ITERS";
        if(item == AlgorithmStop::NONE)                                 return "NONE";
    }
    if constexpr(std::is_same_v<T, MeanType>) {
        if(item == MeanType::ARITHMETIC)                                return "ARITHMETIC";
        if(item == MeanType::GEOMETRIC)                                 return "GEOMETRIC";
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
    if constexpr(std::is_same_v<T, StorageEvent>) {
        if(item == StorageEvent::NONE)                                  return "NONE";
        if(item == StorageEvent::INIT)                                  return "INIT";
        if(item == StorageEvent::MODEL)                                 return "MODEL";
        if(item == StorageEvent::EMIN)                                  return "EMIN";
        if(item == StorageEvent::EMAX)                                  return "EMAX";
        if(item == StorageEvent::PROJECTION)                            return "PROJECTION";
        if(item == StorageEvent::BOND_UPDATE)                           return "BOND_UPDATE";
        if(item == StorageEvent::TRNC_UPDATE)                           return "TRNC_UPDATE";
        if(item == StorageEvent::RBDS_STEP)                             return "RBDS_STEP";
        if(item == StorageEvent::RTES_STEP)                             return "RTES_STEP";
        if(item == StorageEvent::ITERATION)                             return "ITERATION";
        if(item == StorageEvent::FINISHED)                              return "FINISHED";
    }
    if constexpr(std::is_same_v<T, StoragePolicy>) {
        if(item == StoragePolicy::NONE)                                 return "NONE";
        if(item == StoragePolicy::INIT)                                 return "INIT";
        if(item == StoragePolicy::ITER)                                 return "ITER";
        if(item == StoragePolicy::EMIN)                                 return "EMIN";
        if(item == StoragePolicy::EMAX)                                 return "EMAX";
        if(item == StoragePolicy::PROJ)                                 return "PROJ";
        if(item == StoragePolicy::BOND)                                 return "BOND";
        if(item == StoragePolicy::TRNC)                                 return "TRNC";
        if(item == StoragePolicy::FAILURE)                              return "FAILURE";
        if(item == StoragePolicy::SUCCESS)                              return "SUCCESS";
        if(item == StoragePolicy::FINISH)                               return "FINISH";
        if(item == StoragePolicy::ALWAYS)                               return "ALWAYS";
        return "BITFLAG";
    }
    if constexpr(std::is_same_v<T, StateInit>) {
        if(item == StateInit::RANDOM_PRODUCT_STATE)                     return "RANDOM_PRODUCT_STATE";
        if(item == StateInit::RANDOM_ENTANGLED_STATE)                   return "RANDOM_ENTANGLED_STATE";
        if(item == StateInit::RANDOMIZE_PREVIOUS_STATE)                 return "RANDOMIZE_PREVIOUS_STATE";
        if(item == StateInit::MIDCHAIN_SINGLET_NEEL_STATE)              return "MIDCHAIN_SINGLET_NEEL_STATE";
        if(item == StateInit::PRODUCT_STATE_DOMAIN_WALL)                return "PRODUCT_STATE_DOMAIN_WALL";
        if(item == StateInit::PRODUCT_STATE_ALIGNED)                    return "PRODUCT_STATE_ALIGNED";
        if(item == StateInit::PRODUCT_STATE_NEEL)                       return "PRODUCT_STATE_NEEL";
        if(item == StateInit::PRODUCT_STATE_NEEL_SHUFFLED)              return "PRODUCT_STATE_NEEL_SHUFFLED";
        if(item == StateInit::PRODUCT_STATE_NEEL_DISLOCATED)            return "PRODUCT_STATE_NEEL_DISLOCATED";
        if(item == StateInit::PRODUCT_STATE_PATTERN)                    return "PRODUCT_STATE_PATTERN";
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
        if(item == flbit_task::INIT_RANDOMIZE_MODEL)                                return "INIT_RANDOMIZE_MODEL";
        if(item == flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE)                   return "INIT_RANDOMIZE_INTO_PRODUCT_STATE";
        if(item == flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_SHUFFLED)     return "INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_SHUFFLED";
        if(item == flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_DISLOCATED)   return "INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_DISLOCATED";
        if(item == flbit_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE)                 return "INIT_RANDOMIZE_INTO_ENTANGLED_STATE";
        if(item == flbit_task::INIT_BOND_LIMITS)                                    return "INIT_BOND_LIMITS";
        if(item == flbit_task::INIT_TRNC_LIMITS)                                    return "INIT_TRNC_LIMITS";
        if(item == flbit_task::INIT_WRITE_MODEL)                                    return "INIT_WRITE_MODEL";
        if(item == flbit_task::INIT_CLEAR_STATUS)                                   return "INIT_CLEAR_STATUS";
        if(item == flbit_task::INIT_CLEAR_CONVERGENCE)                              return "INIT_CLEAR_CONVERGENCE";
        if(item == flbit_task::INIT_DEFAULT)                                        return "INIT_DEFAULT";
        if(item == flbit_task::INIT_GATES)                                          return "INIT_GATES";
        if(item == flbit_task::INIT_TIME)                                           return "INIT_TIME";
        if(item == flbit_task::TRANSFORM_TO_LBIT)                                   return "TRANSFORM_TO_LBIT";
        if(item == flbit_task::TRANSFORM_TO_REAL)                                   return "TRANSFORM_TO_REAL";
        if(item == flbit_task::TIME_EVOLVE)                                         return "TIME_EVOLVE";
        if(item == flbit_task::POST_WRITE_RESULT)                                   return "POST_WRITE_RESULT";
        if(item == flbit_task::POST_PRINT_RESULT)                                   return "POST_PRINT_RESULT";
        if(item == flbit_task::POST_PRINT_TIMERS)                                   return "POST_PRINT_TIMERS";
        if(item == flbit_task::POST_DEFAULT)                                        return "POST_DEFAULT";
        if(item == flbit_task::TIMER_RESET)                                         return "TIMER_RESET";
    }
    if constexpr(std::is_same_v<T,xdmrg_task>){
        if(item == xdmrg_task::INIT_RANDOMIZE_MODEL)                   return "INIT_RANDOMIZE_MODEL";
        if(item == xdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE)      return "INIT_RANDOMIZE_INTO_PRODUCT_STATE";
        if(item == xdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE)    return "INIT_RANDOMIZE_INTO_ENTANGLED_STATE";
        if(item == xdmrg_task::INIT_RANDOMIZE_FROM_CURRENT_STATE)      return "INIT_RANDOMIZE_FROM_CURRENT_STATE";
        if(item == xdmrg_task::INIT_BOND_LIMITS)                       return "INIT_BOND_LIMITS";
        if(item == xdmrg_task::INIT_TRNC_LIMITS)                       return "INIT_TRNC_LIMITS";
        if(item == xdmrg_task::INIT_ENERGY_TARGET)                     return "INIT_ENERGY_TARGET";
        if(item == xdmrg_task::INIT_WRITE_MODEL)                       return "INIT_WRITE_MODEL";
        if(item == xdmrg_task::INIT_CLEAR_STATUS)                      return "INIT_CLEAR_STATUS";
        if(item == xdmrg_task::INIT_CLEAR_CONVERGENCE)                 return "INIT_CLEAR_CONVERGENCE";
        if(item == xdmrg_task::INIT_DEFAULT)                           return "INIT_DEFAULT";
        if(item == xdmrg_task::FIND_ENERGY_RANGE)                      return "FIND_ENERGY_RANGE";
        if(item == xdmrg_task::FIND_EXCITED_STATE)                     return "FIND_EXCITED_STATE";
        if(item == xdmrg_task::POST_WRITE_RESULT)                      return "POST_WRITE_RESULT";
        if(item == xdmrg_task::POST_PRINT_RESULT)                      return "POST_PRINT_RESULT";
        if(item == xdmrg_task::POST_PRINT_TIMERS)                      return "POST_PRINT_TIMERS";
        if(item == xdmrg_task::POST_RBDS_ANALYSIS)                     return "POST_RBDS_ANALYSIS";
        if(item == xdmrg_task::POST_RTES_ANALYSIS)                     return "POST_RTES_ANALYSIS";
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
    if constexpr(std::is_same_v<T,OptFunc>){
        if(item == OptFunc::ENERGY  )                                  return "ENERGY";
        if(item == OptFunc::VARIANCE)                                  return "VARIANCE";
        if(item == OptFunc::OVERLAP)                                   return "OVERLAP";
    }
    if constexpr(std::is_same_v<T,OptAlgo>){
        if(item == OptAlgo::DIRECT)                                    return "DIRECT";
        if(item == OptAlgo::SUBSPACE)                                  return "SUBSPACE";
        if(item == OptAlgo::SHIFTINV)                                  return "SHIFTINV";
        if(item == OptAlgo::MPSEIGS)                                   return "MPSEIGS";
    }
    if constexpr(std::is_same_v<T,OptSolver>){
        if(item == OptSolver::EIG)                                     return "EIG";
        if(item == OptSolver::EIGS)                                    return "EIGS";
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
    if constexpr(std::is_same_v<T,OptMark>){
        if(item == OptMark::PASS)                                      return "PASS";
        if(item == OptMark::FAIL)                                      return "FAIL";
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
    return "UNRECOGNIZED ENUM";
}




template<typename T>
constexpr auto sv2enum(std::string_view item) {
    static_assert(enum_sfinae::is_any_v<T,
        AlgorithmType,
        AlgorithmStop,
        MeanType,
        MultisiteMove,
        MultisiteWhen,
        MultisiteGrow,
        SVDLibrary,
        UpdateWhen,
        GateMove,
        GateOp,
        CircuitOp,
        LbitCircuitGateMatrixKind,
        LbitCircuitGateWeightKind,
        ModelType,
        EdgeStatus,
        TimeScale,
        StorageEvent,
        StoragePolicy,
        CopyPolicy,
        ResetReason,
        NormPolicy,
        FileCollisionPolicy,
        FileResumePolicy,
        LogPolicy,
        RandomizerMode,
        OptType,
        OptFunc,
        OptSolver,
        OptRitz,
        OptWhen,
        OptExit,
        OptMark,
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
        return "AlgorithmType::UNDEFINED";
    }
    if constexpr(std::is_same_v<T, MultisiteMove>) {
        if(item == "ONE")                                   return MultisiteMove::ONE;
        if(item == "MID")                                   return MultisiteMove::MID;
        if(item == "MAX")                                   return MultisiteMove::MAX;
    }
    if constexpr(std::is_same_v<T, MultisiteWhen>) {
        if(item == "NEVER")                                 return MultisiteWhen::NEVER;
        if(item == "STUCK")                                 return MultisiteWhen::STUCK;
        if(item == "SATURATED")                             return MultisiteWhen::SATURATED;
        if(item == "ALWAYS")                                return MultisiteWhen::ALWAYS;
    }
    if constexpr(std::is_same_v<T, MultisiteGrow>) {
        if(item == "OFF")                                   return MultisiteGrow::OFF;
        if(item == "ON")                                    return MultisiteGrow::ON;
    }
    if constexpr(std::is_same_v<T, OptRitz>) {
        if(item == "NONE")                                  return OptRitz::NONE;
        if(item == "SR")                                    return OptRitz::SR;
        if(item == "LR")                                    return OptRitz::LR;
        if(item == "SM")                                    return OptRitz::SM;
        if(item == "IS")                                    return OptRitz::IS;
        if(item == "TE")                                    return OptRitz::TE;
    }
    if constexpr(std::is_same_v<T, SVDLibrary>) {
        if(item == "EIGEN")                                 return SVDLibrary::EIGEN;
        if(item == "LAPACKE")                               return SVDLibrary::LAPACKE;
        if(item == "RSVD")                                  return SVDLibrary::RSVD;
    }
    if constexpr(std::is_same_v<T, UpdateWhen>) {
        if(item == "NEVER")                                 return UpdateWhen::NEVER;
        if(item == "TRUNCATED")                             return UpdateWhen::TRUNCATED;
        if(item == "STUCK")                                 return UpdateWhen::STUCK;
        if(item == "SATURATED")                             return UpdateWhen::SATURATED;
        if(item == "ITERATION")                             return UpdateWhen::ITERATION;
    }
    if constexpr(std::is_same_v<T, GateMove>) {
        if(item == "OFF")                                   return GateMove::OFF;
        if(item == "ON")                                    return GateMove::ON;
        if(item == "AUTO")                                  return GateMove::AUTO;
    }
    if constexpr(std::is_same_v<T, GateOp>) {
        if(item == "NONE")                                  return GateOp::NONE;
        if(item == "CNJ")                                   return GateOp::CNJ;
        if(item == "ADJ")                                   return GateOp::ADJ;
        if(item == "TRN")                                   return GateOp::TRN;
    }
    if constexpr(std::is_same_v<T, CircuitOp>) {
        if(item == "NONE")                                  return CircuitOp::NONE;
        if(item == "ADJ")                                   return CircuitOp::ADJ;
        if(item == "TRN")                                   return CircuitOp::TRN;
    }
    if constexpr(std::is_same_v<T, LbitCircuitGateMatrixKind>) {
        if(item == "MATRIX_V1")                             return LbitCircuitGateMatrixKind::MATRIX_V1;
        if(item == "MATRIX_V2")                             return LbitCircuitGateMatrixKind::MATRIX_V2;
        if(item == "MATRIX_V3")                             return LbitCircuitGateMatrixKind::MATRIX_V3;

    }
    if constexpr(std::is_same_v<T, LbitCircuitGateWeightKind>) {
        if(item == "IDENTITY")                              return LbitCircuitGateWeightKind::IDENTITY;
        if(item == "EXPDECAY")                              return LbitCircuitGateWeightKind::EXPDECAY;
    }
    if constexpr(std::is_same_v<T, ModelType>) {
        if(item == "ising_tf_rf")                           return ModelType::ising_tf_rf;
        if(item == "ising_sdual")                           return ModelType::ising_sdual;
        if(item == "ising_majorana")                        return ModelType::ising_majorana;
        if(item == "lbit")                                  return ModelType::lbit;
        if(item == "xxz")                                   return ModelType::xxz;
    }
    if constexpr(std::is_same_v<T, EdgeStatus>) {
        if(item == "STALE")                                 return EdgeStatus::STALE ;
        if(item == "FRESH")                                 return EdgeStatus::FRESH ;
    }
    if constexpr(std::is_same_v<T, TimeScale>) {
        if(item == "LINSPACED")                             return TimeScale::LINSPACED ;
        if(item == "LOGSPACED")                             return TimeScale::LOGSPACED ;
    }
    if constexpr(std::is_same_v<T, AlgorithmStop>) {
        if(item == "SUCCESS")                               return AlgorithmStop::SUCCESS;
        if(item == "SATURATED")                             return AlgorithmStop::SATURATED;
        if(item == "MAX_ITERS")                             return AlgorithmStop::MAX_ITERS;
        if(item == "NONE")                                  return AlgorithmStop::NONE;
    }
    if constexpr(std::is_same_v<T, MeanType>) {
        if(item == "ARITHMETIC")                            return MeanType::ARITHMETIC;
        if(item == "GEOMETRIC")                             return MeanType::GEOMETRIC;
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
    if constexpr(std::is_same_v<T, StorageEvent>) {
        if(item == "NONE")                                  return  StorageEvent::NONE;
        if(item == "INIT")                                  return  StorageEvent::INIT;
        if(item == "MODEL")                                 return  StorageEvent::MODEL;
        if(item == "EMIN")                                  return  StorageEvent::EMIN;
        if(item == "EMAX")                                  return  StorageEvent::EMAX;
        if(item == "PROJECTION")                            return  StorageEvent::PROJECTION;
        if(item == "BOND_UPDATE")                           return  StorageEvent::BOND_UPDATE;
        if(item == "TRNC_UPDATE")                           return  StorageEvent::TRNC_UPDATE;
        if(item == "RBDS_STEP")                             return  StorageEvent::RBDS_STEP;
        if(item == "RTES_STEP")                             return  StorageEvent::RTES_STEP;
        if(item == "ITERATION")                             return  StorageEvent::ITERATION;
        if(item == "FINISHED")                              return  StorageEvent::FINISHED;
    }
    if constexpr(std::is_same_v<T, StoragePolicy>) {
        auto policy = StoragePolicy::NONE;
        if(item.find("INIT") != std::string_view::npos) policy |= StoragePolicy::INIT;
        if(item.find("ITER") != std::string_view::npos) policy |= StoragePolicy::ITER;
        if(item.find("EMIN") != std::string_view::npos) policy |= StoragePolicy::EMIN;
        if(item.find("EMAX") != std::string_view::npos) policy |= StoragePolicy::EMAX;
        if(item.find("PROJ") != std::string_view::npos) policy |= StoragePolicy::PROJ;
        if(item.find("BOND") != std::string_view::npos) policy |= StoragePolicy::BOND;
        if(item.find("TRNC") != std::string_view::npos) policy |= StoragePolicy::TRNC;
        if(item.find("FAIL") != std::string_view::npos) policy |= StoragePolicy::FAILURE;
        if(item.find("SUCC") != std::string_view::npos) policy |= StoragePolicy::SUCCESS;
        if(item.find("FINI") != std::string_view::npos) policy |= StoragePolicy::FINISH;
        if(item.find("ALWA") != std::string_view::npos) policy |= StoragePolicy::ALWAYS;
        if(item.find("REPL") != std::string_view::npos) policy |= StoragePolicy::REPLACE;
        return policy;
    }
    if constexpr(std::is_same_v<T, StateInit>) {
        if(item == "RANDOM_PRODUCT_STATE")                  return StateInit::RANDOM_PRODUCT_STATE;
        if(item == "RANDOM_ENTANGLED_STATE")                return StateInit::RANDOM_ENTANGLED_STATE;
        if(item == "RANDOMIZE_PREVIOUS_STATE")              return StateInit::RANDOMIZE_PREVIOUS_STATE;
        if(item == "MIDCHAIN_SINGLET_NEEL_STATE")           return StateInit::MIDCHAIN_SINGLET_NEEL_STATE;
        if(item == "PRODUCT_STATE_DOMAIN_WALL")             return StateInit::PRODUCT_STATE_DOMAIN_WALL;
        if(item == "PRODUCT_STATE_ALIGNED")                 return StateInit::PRODUCT_STATE_ALIGNED;
        if(item == "PRODUCT_STATE_NEEL")                    return StateInit::PRODUCT_STATE_NEEL;
        if(item == "PRODUCT_STATE_NEEL_SHUFFLED")           return StateInit::PRODUCT_STATE_NEEL_SHUFFLED;
        if(item == "PRODUCT_STATE_NEEL_DISLOCATED")         return StateInit::PRODUCT_STATE_NEEL_DISLOCATED;
        if(item == "PRODUCT_STATE_PATTERN")                 return StateInit::PRODUCT_STATE_PATTERN;
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
        if(item == "POST_RBDS_ANALYSIS")                    return fdmrg_task::POST_RBDS_ANALYSIS;
        if(item == "POST_RTES_ANALYSIS")                    return fdmrg_task::POST_RTES_ANALYSIS;
        if(item == "POST_DEFAULT")                          return fdmrg_task::POST_DEFAULT;
        if(item == "TIMER_RESET")                           return fdmrg_task::TIMER_RESET;
    }

    if constexpr(std::is_same_v<T,flbit_task>){
        if(item == "INIT_RANDOMIZE_MODEL")                              return flbit_task::INIT_RANDOMIZE_MODEL;
        if(item == "INIT_RANDOMIZE_INTO_PRODUCT_STATE")                 return flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE;
        if(item == "INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_SHUFFLED")   return flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_SHUFFLED;
        if(item == "INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_DISLOCATED") return flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE_NEEL_DISLOCATED;
        if(item == "INIT_RANDOMIZE_INTO_ENTANGLED_STATE")               return flbit_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE;
        if(item == "INIT_RANDOMIZE_INTO_PRODUCT_STATE_PATTERN")         return flbit_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE_PATTERN;
        if(item == "INIT_RANDOMIZE_INTO_MIDCHAIN_SINGLET_NEEL_STATE")   return flbit_task::INIT_RANDOMIZE_INTO_MIDCHAIN_SINGLET_NEEL_STATE;
        if(item == "INIT_BOND_LIMITS")                                  return flbit_task::INIT_BOND_LIMITS;
        if(item == "INIT_TRNC_LIMITS")                                  return flbit_task::INIT_TRNC_LIMITS;
        if(item == "INIT_WRITE_MODEL")                                  return flbit_task::INIT_WRITE_MODEL;
        if(item == "INIT_CLEAR_STATUS")                                 return flbit_task::INIT_CLEAR_STATUS;
        if(item == "INIT_CLEAR_CONVERGENCE")                            return flbit_task::INIT_CLEAR_CONVERGENCE;
        if(item == "INIT_DEFAULT")                                      return flbit_task::INIT_DEFAULT;
        if(item == "INIT_GATES")                                        return flbit_task::INIT_GATES;
        if(item == "INIT_TIME")                                         return flbit_task::INIT_TIME;
        if(item == "TRANSFORM_TO_LBIT")                                 return flbit_task::TRANSFORM_TO_LBIT;
        if(item == "TRANSFORM_TO_REAL")                                 return flbit_task::TRANSFORM_TO_REAL;
        if(item == "TIME_EVOLVE")                                       return flbit_task::TIME_EVOLVE;
        if(item == "POST_WRITE_RESULT")                                 return flbit_task::POST_WRITE_RESULT;
        if(item == "POST_PRINT_RESULT")                                 return flbit_task::POST_PRINT_RESULT;
        if(item == "POST_PRINT_TIMERS")                                 return flbit_task::POST_PRINT_TIMERS;
        if(item == "POST_DEFAULT")                                      return flbit_task::POST_DEFAULT;
        if(item == "TIMER_RESET")                                       return flbit_task::TIMER_RESET;
    }
    if constexpr(std::is_same_v<T,xdmrg_task>){
        if(item == "INIT_RANDOMIZE_MODEL")                  return xdmrg_task::INIT_RANDOMIZE_MODEL;
        if(item == "INIT_RANDOMIZE_INTO_PRODUCT_STATE")     return xdmrg_task::INIT_RANDOMIZE_INTO_PRODUCT_STATE;
        if(item == "INIT_RANDOMIZE_INTO_ENTANGLED_STATE")   return xdmrg_task::INIT_RANDOMIZE_INTO_ENTANGLED_STATE;
        if(item == "INIT_RANDOMIZE_FROM_CURRENT_STATE")     return xdmrg_task::INIT_RANDOMIZE_FROM_CURRENT_STATE;
        if(item == "INIT_BOND_LIMITS")                      return xdmrg_task::INIT_BOND_LIMITS;
        if(item == "INIT_TRNC_LIMITS")                      return xdmrg_task::INIT_TRNC_LIMITS;
        if(item == "INIT_ENERGY_TARGET")                    return xdmrg_task::INIT_ENERGY_TARGET;
        if(item == "INIT_WRITE_MODEL")                      return xdmrg_task::INIT_WRITE_MODEL;
        if(item == "INIT_CLEAR_STATUS")                     return xdmrg_task::INIT_CLEAR_STATUS;
        if(item == "INIT_CLEAR_CONVERGENCE")                return xdmrg_task::INIT_CLEAR_CONVERGENCE;
        if(item == "INIT_DEFAULT")                          return xdmrg_task::INIT_DEFAULT;
        if(item == "FIND_ENERGY_RANGE")                     return xdmrg_task::FIND_ENERGY_RANGE;
        if(item == "FIND_EXCITED_STATE")                    return xdmrg_task::FIND_EXCITED_STATE;
        if(item == "POST_WRITE_RESULT")                     return xdmrg_task::POST_WRITE_RESULT;
        if(item == "POST_PRINT_RESULT")                     return xdmrg_task::POST_PRINT_RESULT;
        if(item == "POST_PRINT_TIMERS")                     return xdmrg_task::POST_PRINT_TIMERS;
        if(item == "POST_RBDS_ANALYSIS")                    return xdmrg_task::POST_RBDS_ANALYSIS;
        if(item == "POST_RTES_ANALYSIS")                    return xdmrg_task::POST_RTES_ANALYSIS;
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
    if constexpr(std::is_same_v<T,OptFunc>){
        if(item == "ENERGY")                                return OptFunc::ENERGY;
        if(item == "VARIANCE")                              return OptFunc::VARIANCE;
        if(item == "OVERLAP")                               return OptFunc::OVERLAP;
    }
    if constexpr(std::is_same_v<T,OptAlgo>){
        if(item == "DIRECT")                                return OptAlgo::DIRECT;
        if(item == "SUBSPACE")                              return OptAlgo::SUBSPACE;
        if(item == "SHIFTINV")                              return OptAlgo::SHIFTINV;
        if(item == "MPSEIGS")                               return OptAlgo::MPSEIGS;
    }
    if constexpr(std::is_same_v<T,OptSolver>){
        if(item == "EIG")                                  return OptSolver::EIG;
        if(item == "EIGS")                                 return OptSolver::EIGS;
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
    if constexpr(std::is_same_v<T,OptMark>){
        if(item == "PASS")                                  return OptMark::PASS;
        if(item == "FAIL")                                  return OptMark::FAIL;
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
    auto           pairgen = [](const std::string &name) -> std::pair<std::string, T> { return {name, sv2enum<T>(name)}; };
    return enumarray_t<T, num>{pairgen(names)...};
}

template<typename T, typename... Args>
constexpr auto mapEnum2Str(Args... enums) {
    constexpr auto num     = sizeof...(enums);
    auto           pairgen = [](const T &e) -> std::pair<std::string, T> { return {std::string(enum2sv(e)), e}; };
    return enumarray_t<T, num>{pairgen(enums)...};
}

// inline auto ModelType_s2e = mapEnum2Str<ModelType>(ModelType::ising_tf_rf, ModelType::ising_sdual, ModelType::ising_majorana, ModelType::lbit);
//  inline auto ModelType_s2e = mapEnum2Str<ModelType>(ModelType::ising_tf_rf, ModelType::ising_sdual, ModelType::ising_majorana, ModelType::lbit);
