#pragma once
#include <array>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

enum class AlgorithmStop : int { SUCCESS, SATURATED, MAX_ITERS, NONE };
enum class AlgorithmType : int { iDMRG, fDMRG, xDMRG, iTEBD, fLBIT, ANY };
enum class ModelType { ising_tf_rf, ising_sdual, ising_majorana, lbit, xxz };
enum class MpoCompress {
    NONE, /*!< Do not compress */
    SVD,  /*!< Use SVD on each mpo */
    DPL   /*!< Deparallelization: removes parallel columns/rows from each mpo */
};
enum class MposWithEdges { OFF, ON };
enum class MeanType { ARITHMETIC, GEOMETRIC };

/*! \brief Controls the dynamical adjustments of the dmrg block size (number of sites in a dmrg optimization step).
 *  The default block size is always settings::strategy::dmrg_min_blocksize.
 *  Combine these bitflags to increase the block size on a given status.
 *
 *  Examples:
 *      - MIN                   : Set the block size to settings::strategy::dmrg_min_blocksize at all steps
 *      - MAX                   : Set the block size to settings::strategy::dmrg_max_blocksize at all steps
 *      - MAX|SATURATED         : Set the block size to settings::strategy::dmrg_max_blocksize if the algorithm has saturated
 *      - ICOM|STUCK|FIN_TRNC : Set the block size to the ceil of "information_center_of_mass" if the algorithm is stuck and the trnc_lim ==
 * settings::precision::svd_truncation_min
 *
 *  Note that settings::strategy::dmrg_min_blocksize and settings::strategy::dmrg_max_blocksize are hard limits on the block size.
 *  This means that "ICOM" sets a block size within those limits.
 */

/*! \brief Why an mps merge event is invoked
 * Currently, this is used to check if an optimized multisite_mps was passed, in which case we need to store the truncation errors (otherwise they are
 * discarded)
 */
enum class MergeEvent {
    MOVE, /*!< Moved the center position */
    NORM, /*!< Normalized the state */
    SWAP, /*!< Swapped sites (e.g. in fLBIT) */
    GATE, /*!< Applied a gate, e.g. during time evolution */
    OPT,  /*!< Optimized some sites, e.g. during the main DMRG update */
    EXP,  /*!< Subspace expansion (aka perurbation, noise, enrichment) to increase the bond dimension */
};

enum class BlockSizePolicy {
    MIN       = 0,             /*!< Set the block size to dmrg_max_blocksize on special status (set by additional bitflags) */
    MAX       = 1,             /*!< Set the block size to dmrg_max_blocksize on special status (set by additional bitflags) */
    ICOM      = 2,             /*!< Set the block size to the ceil of "information_center_of_mass" on special status (set by additional bitflags) */
    ICOMPLUS1 = 4,             /*!< Set the block size to the ceil of "information_center_of_mass + 1" on special status (set by additional bitflags) */
    ICOM150   = 8,             /*!< Set the block size to the ceil of "information_center_of_mass + 1" on special status (set by additional bitflags) */
    ICOM200   = 16,            /*!< Set the block size to the ceil of "information_center_of_mass + 1" on special status (set by additional bitflags) */
    SAT_ENT   = 32,            /*!< Set the block size when the entanglement entropy has saturated */
    SAT_ICOM  = 64,            /*!< Set the block size when the information center of mass has saturated */
    SAT_VAR   = 128,           /*!< Set the block size when the energy variance has saturated */
    SAT_ALGO  = 256,           /*!< Set the block size when the algorithm status == saturated (implies all saturations) */
    STK_ALGO  = 512,           /*!< Set the block size when the algorithm status == stuck ( */
    FIN_BOND  = 1024,          /*!< Require that the bond dimension has reached its final (maximum) value before increasing the block size */
    FIN_TRNC  = 2048,          /*!< Require that the truncation error has reached its final (minimum) value before increasing the block size   */
    EXP       = 4096,          /*!< Use the enlarged blocksize during environment expansion  */
    OPT       = 8192,          /*!< Use the enlarged blocksize during the optimization step */
    DEFAULT   = ICOM150 | EXP, /*!< The default choice usually works well */
    allow_bitops
};

/*! \brief Policy that determines when a quantity should increase or improve
 *  For example, when the number of eigs iterations should increase, or when
 *  the mps block size should increase, or jacobi preconditioner block size
 */
enum class GainPolicy : int {

    NEVER     = 0,  /*!< Never increase eigs iterations  */
    HALFSWEEP = 1,  /*!< every halfsweep    */
    FULLSWEEP = 2,  /*!< every halfsweep    */
    SAT_VAR   = 4,  /*!< when the energy variance has saturated  */
    SAT_ALGO  = 8,  /*!< when the algorithm has saturated  */
    STK_ALGO  = 16, /*!< when the algorithm is stuck  */
    FIN_BOND  = 32, /*!< only when the bond dimension has reached its maximum   */
    FIN_TRNC  = 64, /*!< only when the truncation error has reached its minimum   */
    allow_bitops
};

enum class UpdatePolicy {
    NEVER     = 0,  /*!< Never update */
    WARMUP    = 1,  /*!< Update during warmup */
    HALFSWEEP = 2,  /*!< Update every iteration */
    FULLSWEEP = 4,  /*!< Update every second iteration (left to right + right to left sweep) */
    TRUNCATED = 8,  /*!< Update whenever the state is truncated */
    SAT_VAR   = 16, /*!< when the algorithm has saturated  <*/
    SAT_ALGO  = 32, /*!< Update when the algorithm is saturated */
    STK_ALGO  = 64, /*!< Update when the algorithm is stuck */
    allow_bitops
};

/*! \brief Controls how saturation is checked.
 *  Given a series
 *      x = [x0,x1,x2...xN],
 *  we generate a new sequence
 *      y = [f0, f1, f2... fN],
 *  where fk = function(xk,...,xN).
 *  That is, we apply a function f on
 *      x,
 *      x with first element removed,
 *      x with first and second elements removed,
 *  and so on.
 *  The function f is either
 *       val: the standard deviation on xk...xN
 *       mov: the standard deviation on the moving average of xk...xN
 *       min: the standard deviation on the minimum of xk...xN
 *       max: the standard deviation on the maximum of xk...xN
 *       mid: the standard deviation on the midpoint between min and max
 */
enum class SaturationPolicy {
    val = 0,   /*!< Check the standard deviation on xk...xN */
    avg = 1,   /*!< Check the standard deviation on xk...xN */
    med = 2,   /*!< Check the standard deviation on xk...xN */
    mov = 4,   /*!< Check the standard deviation on the moving average of xk...xN */
    min = 8,   /*!< Check the standard deviation on the minimum of xk...xN */
    max = 16,  /*!< Check the standard deviation on the maximum of xk...xN */
    mid = 32,  /*!< Check the standard deviation on the midpoint between min and max */
    dif = 64,  /*!< Check the standard deviation on the midpoint between min and max */
    log = 128, /*!< Transform x -> log(x) first */
    allow_bitops
};

enum class SVDLibrary { EIGEN, LAPACKE, RSVD };
enum class GateMove { OFF, ON, AUTO };
enum class GateOp { NONE, CNJ, ADJ, TRN };
enum class CircuitOp { NONE, ADJ, TRN };
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

/*! When to do projections to a symmetry sector */
enum class ProjectionPolicy {
    NEVER     = 0,                        /*!< Never project */
    INIT      = 1,                        /*!< Project after initializing the state  */
    WARMUP    = 2,                        /*!< Project during warmup  */
    STUCK     = 4,                        /*!< Project when the algorithm is stuck */
    ITER      = 8,                        /*!< Project every iteration */
    CONVERGED = 16,                       /*!< Project on converged iterations */
    FINISHED  = 32,                       /*!< Project the finished state during postprocessing */
    FORCE     = 64,                       /*!< Project even if not needed */
    DEFAULT   = INIT | STUCK | CONVERGED, /*!< Default policy for multisite dmrg steps */
    allow_bitops
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
    allow_bitops
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
    RBDS    = 4096, /*!< Store rbds steps */
    RTES    = 8192, /*!< Store rtes steps */
    allow_bitops
};
enum class CopyPolicy { FORCE, TRY, OFF };
enum class ResetReason { INIT, FIND_WINDOW, SATURATED, NEW_STATE, BOND_UPDATE };

enum class EnvExpandMode : int {
    NONE     = 0,  /*!< No environment expansion */
    H1       = 1,  /*!< environment expansion using H (valid with SSITE, ignored otherwise) */
    H2       = 2,  /*!< environment expansion using H²  (valid with SSITE, ignored otherwise)  */
    SSITE    = 4,  /*!< strictly single site expansion */
    NSITE    = 8,  /*!< nsite expansion (given by the local info scale "icom") */
    BACKWARD = 16, /*!< Expand the trailing environment as in the original DMRG3S algorithm */
    FORWARD  = 32, /*!< Expand the forward environment before optimization */
    DEFAULT  = H1 | H2 | NSITE | FORWARD,
    allow_bitops
};

enum class NormPolicy { ALWAYS, IFNEEDED }; // Rules of engagement
enum class ResumePolicy {
    IF_MAX_ITERS,
    IF_SATURATED,
    IF_UNSUCCESSFUL,
    IF_SUCCESSFUL,
    ALWAYS,
};
enum class FileCollisionPolicy {
    RESUME, /*!< If finished -> exit, else resume simulation from the latest "FULL" storage state. Throw if none is found. */
    BACKUP, /*!< Backup the existing file by appending .bak, then start with a new file. */
    RENAME, /*!< Rename the current file by appending .# to avoid collision with existing. */
    REVIVE, /*!< Try RESUME, but do REPLACE on error instead of throwing */
    REPLACE /*!< Just erase/truncate the existing file and start from the beginning. */
};

enum class FileResumePolicy { FULL, FAST };
enum class LogPolicy {
    SILENT, /*!< Never log */
    DEBUG,  /*!< Log on debug runs */
    VERBOSE /*!< Always log */
};
enum class RandomizerMode { SHUFFLE, SELECT1, ASIS };
enum class OptType { REAL, CPLX };

/*! Choose the algorithm for the optimization.
 *
 * In the descriptions below, H and |v> are the local effective Hamiltonians and states.
 */
enum class OptAlgo {
    DMRG,         /*!< Plain DMRG that solves the eigenvalue problem H|v>=E|v> for some given  */
    DMRGX,        /*!< Find an eigenvector of H that maximizes the overlap:  |v_new> = max_k <v_old | v_k> */
    HYBRID_DMRGX, /*!< Minimize the variance of |v_new> = sum_{k in K} λ_k|v_k>..., where K is a basis of H_eff eigenvectors satisfying <v_old | v_k> != 0 */
    XDMRG,        /*!< Solve the folded eigenvalue problem (H-Eshift)²|v>  */
    GDMRG,        /*!< Solve the generalized shift-invert problem (H-Eshift)|v> = (1/E)(H-Eshift)²|v> */
};

enum class OptSolver {
    EIG, /*!< Apply the eigensolver directly on H|psi> or (H-E)^2|psi> */
    EIGS /*!< Apply the eigensolver directly on H|psi> or (H-E)^2|psi> */
};

/*! Choose the target eigenpair */
enum class OptRitz {
    NONE, /*!< No eigenpair is targeted (e.g. time evolution) */
    LR,   /*!< Largest Real eigenvalue */
    LM,   /*!< Largest Absolute eigenvalue */
    SR,   /*!< Smallest Real eigenvalue */
    SM,   /*!< Smallest magnitude eigenvalue. MPO² Energy shift == 0. Use this to find an eigenstate with energy closest to 0) */
    IS,   /*!< Initial State energy. Energy shift == Initial state energy. Targets an eigenstate with energy near that of the initial state */
    TE,   /*!< Target Energy density in normalized units [0,1]. Energy shift == settings::xdmrg::energy_density_target * (EMIN+EMAX) + EMIN.  */
    // CE,    /*!< Closest to the current energy  */
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

template<typename T>
concept enum_is_bitflag_v = requires(T m) {
    { T::allow_bitops };
};

template<typename E>
requires enum_is_bitflag_v<E>
constexpr auto operator|(E lhs, E rhs) noexcept -> decltype(E::allow_bitops) {
    using U = std::underlying_type_t<E>;
    return static_cast<E>(static_cast<U>(lhs) | static_cast<U>(rhs));
}
template<typename E>
requires enum_is_bitflag_v<E>
constexpr auto operator&(E lhs, E rhs) noexcept -> decltype(E::allow_bitops) {
    using U = std::underlying_type_t<E>;
    return static_cast<E>(static_cast<U>(lhs) & static_cast<U>(rhs));
}
template<typename E>
requires enum_is_bitflag_v<E>
constexpr auto operator|=(E &lhs, E rhs) noexcept -> decltype(E::allow_bitops) {
    lhs = lhs | rhs;
    return lhs;
}

template<typename E>
inline bool has_flag(E target, E check) noexcept {
    using U = std::underlying_type_t<E>;
    return (static_cast<U>(target) & static_cast<U>(check)) == static_cast<U>(check);
}

template<typename E, typename... Args>
requires std::conjunction_v<std::is_same<E, Args>...>
bool has_any_flags(E target, Args &&...check) {
    using U = std::underlying_type_t<E>;
    return (((static_cast<U>(target) & static_cast<U>(check)) == static_cast<U>(check)) || ...);
}
template<typename E, typename... Args>
requires std::conjunction_v<std::is_same<E, Args>...>
bool has_all_flags(E target, Args &&...check) {
    using U = std::underlying_type_t<E>;
    return (((static_cast<U>(target) & static_cast<U>(check)) == static_cast<U>(check)) && ...);
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
    static_assert(enum_sfinae::is_any_v<T, OptExit, OptWhen, GainPolicy, StoragePolicy, BlockSizePolicy, UpdatePolicy, SaturationPolicy, ProjectionPolicy,
                                        EnvExpandMode>);

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
        if(has_flag(item, StoragePolicy::RBDS)) v.emplace_back("RBDS");
        if(has_flag(item, StoragePolicy::RTES)) v.emplace_back("RTES");
        if(v.empty()) return "NONE";
    }
    if constexpr(std::is_same_v<T, BlockSizePolicy>) {
        if(has_flag(item, BlockSizePolicy::MAX)) v.emplace_back("MAX");
        if(has_flag(item, BlockSizePolicy::ICOM)) v.emplace_back("ICOM");
        if(has_flag(item, BlockSizePolicy::ICOMPLUS1)) v.emplace_back("ICOMPLUS1");
        if(has_flag(item, BlockSizePolicy::ICOM150)) v.emplace_back("ICOM150");
        if(has_flag(item, BlockSizePolicy::ICOM200)) v.emplace_back("ICOM200");
        if(has_flag(item, BlockSizePolicy::SAT_ENT)) v.emplace_back("SAT_ENT");
        if(has_flag(item, BlockSizePolicy::SAT_ICOM)) v.emplace_back("SAT_ICOM");
        if(has_flag(item, BlockSizePolicy::SAT_VAR)) v.emplace_back("SAT_VAR");
        if(has_flag(item, BlockSizePolicy::SAT_ALGO)) v.emplace_back("SAT_ALGO");
        if(has_flag(item, BlockSizePolicy::STK_ALGO)) v.emplace_back("STK_ALGO");
        if(has_flag(item, BlockSizePolicy::FIN_BOND)) v.emplace_back("FIN_BOND");
        if(has_flag(item, BlockSizePolicy::FIN_TRNC)) v.emplace_back("FIN_TRNC");
        if(has_flag(item, BlockSizePolicy::EXP)) v.emplace_back("EXP");
        if(has_flag(item, BlockSizePolicy::OPT)) v.emplace_back("OPT");
        if(has_flag(item, BlockSizePolicy::DEFAULT)) v.emplace_back("DEFAULT");
        if(v.empty()) return "MIN";
    }
    if constexpr(std::is_same_v<T, UpdatePolicy>) {
        if(has_flag(item, UpdatePolicy::WARMUP)) v.emplace_back("WARMUP");
        if(has_flag(item, UpdatePolicy::HALFSWEEP)) v.emplace_back("ITERATION");
        if(has_flag(item, UpdatePolicy::FULLSWEEP)) v.emplace_back("FULLSWEEP");
        if(has_flag(item, UpdatePolicy::TRUNCATED)) v.emplace_back("TRUNCATED");
        if(has_flag(item, UpdatePolicy::SAT_ALGO)) v.emplace_back("SATURATED");
        if(has_flag(item, UpdatePolicy::STK_ALGO)) v.emplace_back("STUCK");
        if(v.empty()) return "NEVER";
    }
    if constexpr(std::is_same_v<T, SaturationPolicy>) {
        if(has_flag(item, SaturationPolicy::avg)) v.emplace_back("avg");
        if(has_flag(item, SaturationPolicy::med)) v.emplace_back("med");
        if(has_flag(item, SaturationPolicy::mov)) v.emplace_back("mov");
        if(has_flag(item, SaturationPolicy::min)) v.emplace_back("min");
        if(has_flag(item, SaturationPolicy::max)) v.emplace_back("max");
        if(has_flag(item, SaturationPolicy::mid)) v.emplace_back("mid");
        if(has_flag(item, SaturationPolicy::dif)) v.emplace_back("dif");
        if(has_flag(item, SaturationPolicy::log)) v.emplace_back("log");
        if(v.empty()) return "val";
    }
    if constexpr(std::is_same_v<T, ProjectionPolicy>) {
        if(has_flag(item, ProjectionPolicy::INIT)) v.emplace_back("INIT");
        if(has_flag(item, ProjectionPolicy::WARMUP)) v.emplace_back("WARMUP");
        if(has_flag(item, ProjectionPolicy::STUCK)) v.emplace_back("STUCK");
        if(has_flag(item, ProjectionPolicy::ITER)) v.emplace_back("ITER");
        if(has_flag(item, ProjectionPolicy::CONVERGED)) v.emplace_back("CONVERGED");
        if(has_flag(item, ProjectionPolicy::FINISHED)) v.emplace_back("FINISHED");
        if(has_flag(item, ProjectionPolicy::FORCE)) v.emplace_back("FORCE");
        if(v.empty()) return "NEVER";
    }
    if constexpr(std::is_same_v<T, EnvExpandMode>) {
        if(has_flag(item, EnvExpandMode::H1)) v.emplace_back("H1");
        if(has_flag(item, EnvExpandMode::H2)) v.emplace_back("H2");
        if(has_flag(item, EnvExpandMode::SSITE)) v.emplace_back("SSITE");
        if(has_flag(item, EnvExpandMode::NSITE)) v.emplace_back("NSITE");
        if(has_flag(item, EnvExpandMode::BACKWARD)) v.emplace_back("BACKWARD");
        if(has_flag(item, EnvExpandMode::FORWARD)) v.emplace_back("FORWARD");
        if(has_flag(item, EnvExpandMode::DEFAULT)) v.emplace_back("DEFAULT");
        if(v.empty()) return "NONE";
    }
    if constexpr(std::is_same_v<T, GainPolicy>) {
        if(has_flag(item, GainPolicy::HALFSWEEP)) v.emplace_back("HALFSWEEP");
        if(has_flag(item, GainPolicy::FULLSWEEP)) v.emplace_back("FULLSWEEP");
        if(has_flag(item, GainPolicy::SAT_VAR)) v.emplace_back("SAT_VAR");
        if(has_flag(item, GainPolicy::SAT_ALGO)) v.emplace_back("SAT_ALGO");
        if(has_flag(item, GainPolicy::STK_ALGO)) v.emplace_back("STK_ALGO");
        if(has_flag(item, GainPolicy::FIN_BOND)) v.emplace_back("FIN_BOND");
        if(has_flag(item, GainPolicy::FIN_TRNC)) v.emplace_back("FIN_TRNC");
        if(v.empty()) return "NEVER";
    }
    return std::accumulate(std::begin(v), std::end(v), std::string(),
                           [](const std::string &ss, const std::string &s) { return ss.empty() ? s : ss + "|" + s; });
}

template<typename T>
constexpr std::string_view enum2sv(const T item) noexcept {
    /* clang-format off */
    static_assert(std::is_enum_v<T> and "enum2sv<T>: T must be an enum");
    static_assert(enum_sfinae::is_any_v<T,
        AlgorithmStop,
        AlgorithmType,
        MpoCompress,
        MposWithEdges,
        MeanType,
        MergeEvent,
        BlockSizePolicy,
        SVDLibrary,
        UpdatePolicy,
        SaturationPolicy,
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
        EnvExpandMode,
        NormPolicy,
        ResumePolicy,
        FileCollisionPolicy,
        FileResumePolicy,
        LogPolicy,
        RandomizerMode,
        OptType,
        OptAlgo,
        OptSolver,
        OptRitz,
        OptWhen,
        OptExit,
        OptMark,
        GainPolicy,
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
    if constexpr(std::is_same_v<T, MpoCompress>) switch(item) {
        case MpoCompress::NONE:                                  return "NONE";
        case MpoCompress::SVD:                                   return "SVD";
        case MpoCompress::DPL:                                   return "DPL";
        default: return "MpoCompress::UNDEFINED";
    }
    if constexpr(std::is_same_v<T, MposWithEdges>) switch(item) {
        case MposWithEdges::OFF:                                 return "OFF";
        case MposWithEdges::ON:                                  return "ON";
        default: return "MposWithEdges::UNDEFINED";
    }
    if constexpr(std::is_same_v<T, BlockSizePolicy>) switch(item){
        case BlockSizePolicy::MIN :                              return "MIN";
        case BlockSizePolicy::MAX :                              return "MAX";
        case BlockSizePolicy::ICOM :                             return "ICOM";
        case BlockSizePolicy::ICOMPLUS1 :                        return "ICOMPLUS1";
        case BlockSizePolicy::ICOM150 :                          return "ICOM150";
        case BlockSizePolicy::ICOM200 :                          return "ICOM200";
        case BlockSizePolicy::SAT_ENT :                          return "SAT_ENT";
        case BlockSizePolicy::SAT_ICOM :                         return "SAT_ICOM";
        case BlockSizePolicy::SAT_VAR :                          return "SAT_VAR";
        case BlockSizePolicy::SAT_ALGO:                          return "SAT_ALGO";
        case BlockSizePolicy::STK_ALGO:                          return "STK_ALGO";
        case BlockSizePolicy::FIN_BOND:                          return "FIN_BOND";
        case BlockSizePolicy::FIN_TRNC:                          return "FIN_TRNC";
        case BlockSizePolicy::EXP:                               return "EXP";
        case BlockSizePolicy::OPT:                               return "OPT";
        case BlockSizePolicy::DEFAULT:                           return "DEFAULT";
        default: return "BlockSizePolicy::BITFLAG";
    }
    if constexpr(std::is_same_v<T, OptRitz>) {
        if(item == OptRitz::NONE)                                return "NONE";
        if(item == OptRitz::LR)                                  return "LR";
        if(item == OptRitz::LM)                                  return "LM";
        if(item == OptRitz::SM)                                  return "SM";
        if(item == OptRitz::SM)                                  return "SM";
        if(item == OptRitz::IS)                                  return "IS";
        if(item == OptRitz::TE)                                  return "TE";
    }
    if constexpr(std::is_same_v<T, SVDLibrary>) {
        if(item == SVDLibrary::EIGEN)                            return "EIGEN";
        if(item == SVDLibrary::LAPACKE)                          return "LAPACKE";
        if(item == SVDLibrary::RSVD)                             return "RSVD";
    }
    if constexpr(std::is_same_v<T, UpdatePolicy>) {
        if(item == UpdatePolicy::NEVER)                          return "NEVER";
        if(item == UpdatePolicy::WARMUP)                         return "WARMUP";
        if(item == UpdatePolicy::HALFSWEEP)                      return "ITERATION";
        if(item == UpdatePolicy::FULLSWEEP)                      return "FULLSWEEP";
        if(item == UpdatePolicy::TRUNCATED)                      return "TRUNCATED";
        if(item == UpdatePolicy::SAT_ALGO)                       return "SATURATED";
        if(item == UpdatePolicy::STK_ALGO)                       return "STUCK";
    }
    if constexpr(std::is_same_v<T, SaturationPolicy>) {
        if(item == SaturationPolicy::val)                        return "val";
        if(item == SaturationPolicy::avg)                        return "avg";
        if(item == SaturationPolicy::med)                        return "med";
        if(item == SaturationPolicy::mov)                        return "mov";
        if(item == SaturationPolicy::min)                        return "min";
        if(item == SaturationPolicy::max)                        return "max";
        if(item == SaturationPolicy::mid)                        return "mid";
        if(item == SaturationPolicy::dif)                        return "dif";
        if(item == SaturationPolicy::log)                        return "log";
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
    if constexpr(std::is_same_v<T, ProjectionPolicy>) {
        if(item == ProjectionPolicy::NEVER)                             return "NEVER";
        if(item == ProjectionPolicy::INIT)                              return "INIT";
        if(item == ProjectionPolicy::WARMUP)                            return "WARMUP";
        if(item == ProjectionPolicy::STUCK)                             return "STUCK";
        if(item == ProjectionPolicy::ITER)                              return "ITER";
        if(item == ProjectionPolicy::CONVERGED)                         return "CONVERGED";
        if(item == ProjectionPolicy::FINISHED)                          return "FINISHED";
        if(item == ProjectionPolicy::FORCE)                             return "FORCE";
        if(item == ProjectionPolicy::DEFAULT)                           return "DEFAULT";
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
    if constexpr(std::is_same_v<T, MergeEvent>) {
        if(item == MergeEvent::MOVE)                                    return "MOVE";
        if(item == MergeEvent::NORM)                                    return "NORM";
        if(item == MergeEvent::SWAP)                                    return "SWAP";
        if(item == MergeEvent::GATE)                                    return "GATE";
        if(item == MergeEvent::OPT)                                     return "OPT";
    }
    if constexpr(std::is_same_v<T, ResetReason>) {
        if(item == ResetReason::INIT)                                   return "INIT";
        if(item == ResetReason::FIND_WINDOW)                            return "FIND_WINDOW";
        if(item == ResetReason::SATURATED)                              return "SATURATED";
        if(item == ResetReason::NEW_STATE)                              return "NEW_STATE";
        if(item == ResetReason::BOND_UPDATE)                            return "BOND_UPDATE";
    }
    if constexpr(std::is_same_v<T, EnvExpandMode>) {
        if(item == EnvExpandMode::NONE)                                 return "NONE";
        if(item == EnvExpandMode::H1)                                   return "H1";
        if(item == EnvExpandMode::H2)                                   return "H2";
        if(item == EnvExpandMode::SSITE)                                return "SSITE";
        if(item == EnvExpandMode::NSITE)                                return "NSITE";
        if(item == EnvExpandMode::BACKWARD)                             return "BACKWARD";
        if(item == EnvExpandMode::FORWARD)                              return "FORWARD";
        if(item == EnvExpandMode::DEFAULT)                              return "DEFAULT";
    }
    if constexpr(std::is_same_v<T, NormPolicy>) {
        if(item == NormPolicy::ALWAYS)                                  return "ALWAYS";
        if(item == NormPolicy::IFNEEDED)                                return "IFNEEDED";
    }
    if constexpr(std::is_same_v<T, LogPolicy>) {
        if(item == LogPolicy::SILENT)                                   return "SILENT";
        if(item == LogPolicy::DEBUG)                                    return "DEBUG";
        if(item == LogPolicy::VERBOSE)                                  return "VERBOSE";
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
        if(item == StoragePolicy::REPLACE)                              return "REPLACE";
        if(item == StoragePolicy::RBDS)                                 return "RBDS";
        if(item == StoragePolicy::RTES)                                 return "RTES";
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
    if constexpr(std::is_same_v<T, ResumePolicy>) {
        if(item == ResumePolicy::IF_MAX_ITERS)                          return "IF_MAX_ITERS";
        if(item == ResumePolicy::IF_SATURATED)                          return "IF_SATURATED";
        if(item == ResumePolicy::IF_UNSUCCESSFUL)                       return "IF_UNSUCCESSFUL";
        if(item == ResumePolicy::IF_SUCCESSFUL)                         return "IF_SUCCESSFUL";
        if(item == ResumePolicy::ALWAYS)                                return "ALWAYS";
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
    if constexpr(std::is_same_v<T,OptAlgo>){
        if(item == OptAlgo::DMRG)                                      return "DMRG";
        if(item == OptAlgo::DMRGX)                                     return "DMRGX";
        if(item == OptAlgo::HYBRID_DMRGX)                              return "HYBRID_DMRGX";
        if(item == OptAlgo::XDMRG)                                     return "XDMRG";
        if(item == OptAlgo::GDMRG)                                     return "GDMRG";
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
    if constexpr(std::is_same_v<T,GainPolicy>){
        if(item == GainPolicy::NEVER)                                  return "NEVER";
        if(item == GainPolicy::HALFSWEEP)                              return "HALFSWEEP";
        if(item == GainPolicy::FULLSWEEP)                              return "FULLSWEEP";
        if(item == GainPolicy::SAT_VAR)                                return "SAT_VAR";
        if(item == GainPolicy::SAT_ALGO)                               return "SAT_ALGO";
        if(item == GainPolicy::STK_ALGO)                               return "STK_ALGO";
        if(item == GainPolicy::FIN_BOND)                               return "FIN_BOND";
        if(item == GainPolicy::FIN_TRNC)                               return "FIN_TRNC";
    }
    return "UNRECOGNIZED ENUM";
}




template<typename T>
constexpr auto sv2enum(std::string_view item) {
    static_assert(enum_sfinae::is_any_v<T,
        AlgorithmStop,
        AlgorithmType,
        MpoCompress,
        MposWithEdges,
        MeanType,
        MergeEvent,
        BlockSizePolicy,
        SVDLibrary,
        UpdatePolicy,
        SaturationPolicy,
        GateMove,
        GateOp,
        CircuitOp,
        LbitCircuitGateMatrixKind,
        LbitCircuitGateWeightKind,
        ProjectionPolicy,
        ModelType,
        EdgeStatus,
        TimeScale,
        StorageEvent,
        StoragePolicy,
        CopyPolicy,
        ResetReason,
        EnvExpandMode,
        NormPolicy,
        ResumePolicy,
        FileCollisionPolicy,
        FileResumePolicy,
        LogPolicy,
        RandomizerMode,
        OptAlgo,
        OptType,
        OptSolver,
        OptRitz,
        OptWhen,
        OptExit,
        OptMark,
        GainPolicy,
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
    if constexpr(std::is_same_v<T, MpoCompress>){
        if(item == "NONE")                                  return  MpoCompress::NONE;
        if(item == "SVD")                                   return  MpoCompress::SVD;
        if(item == "DPL")                                   return  MpoCompress::DPL;
    }
    if constexpr(std::is_same_v<T, MposWithEdges>){
        if(item == "OFF")                                   return MposWithEdges::OFF;
        if(item == "ON")                                    return MposWithEdges::ON;
    }
    if constexpr(std::is_same_v<T, BlockSizePolicy>) {
        auto policy = BlockSizePolicy::MIN;
        if(item.find("MAX")         != std::string_view::npos)  policy |= BlockSizePolicy::MAX;
        if(item.find("ICOM")        != std::string_view::npos)  policy |= BlockSizePolicy::ICOM;
        if(item.find("ICOMPLUS1")   != std::string_view::npos)  policy |= BlockSizePolicy::ICOMPLUS1;
        if(item.find("ICOM150")     != std::string_view::npos)  policy |= BlockSizePolicy::ICOM150;
        if(item.find("ICOM200")     != std::string_view::npos)  policy |= BlockSizePolicy::ICOM200;
        if(item.find("SAT_ENT")     != std::string_view::npos)  policy |= BlockSizePolicy::SAT_ENT;
        if(item.find("SAT_ICOM")    != std::string_view::npos)  policy |= BlockSizePolicy::SAT_ICOM;
        if(item.find("SAT_VAR")     != std::string_view::npos)  policy |= BlockSizePolicy::SAT_VAR;
        if(item.find("SAT_ALGO")    != std::string_view::npos)  policy |= BlockSizePolicy::SAT_ALGO;
        if(item.find("STK_ALGO")    != std::string_view::npos)  policy |= BlockSizePolicy::STK_ALGO;
        if(item.find("FIN_BOND")    != std::string_view::npos)  policy |= BlockSizePolicy::FIN_BOND;
        if(item.find("FIN_TRNC")    != std::string_view::npos)  policy |= BlockSizePolicy::FIN_TRNC;
        if(item.find("OPT")         != std::string_view::npos)  policy |= BlockSizePolicy::OPT;
        if(item.find("EXP")         != std::string_view::npos)  policy |= BlockSizePolicy::EXP;
        if(item.find("DEFAULT")     != std::string_view::npos)  policy |= BlockSizePolicy::DEFAULT;
        return policy;

    }
    if constexpr(std::is_same_v<T,GainPolicy>){
        auto policy = GainPolicy::NEVER;
        if(item.find("HALFSWEEP")   != std::string_view::npos) policy |= GainPolicy::HALFSWEEP;
        if(item.find("FULLSWEEP")   != std::string_view::npos) policy |= GainPolicy::FULLSWEEP;
        if(item.find("SAT_VAR")     != std::string_view::npos) policy |= GainPolicy::SAT_VAR;
        if(item.find("SAT_ALGO")    != std::string_view::npos) policy |= GainPolicy::SAT_ALGO;
        if(item.find("STK_ALGO")    != std::string_view::npos) policy |= GainPolicy::STK_ALGO;
        if(item.find("FIN_BOND")    != std::string_view::npos) policy |= GainPolicy::FIN_BOND;
        if(item.find("FIN_TRNC")    != std::string_view::npos) policy |= GainPolicy::FIN_TRNC;


        return policy;
    }
    if constexpr(std::is_same_v<T, OptRitz>) {
        if(item == "NONE")                                  return OptRitz::NONE;
        if(item == "LR")                                    return OptRitz::LR;
        if(item == "LM")                                    return OptRitz::LM;
        if(item == "SR")                                    return OptRitz::SR;
        if(item == "SM")                                    return OptRitz::SM;
        if(item == "IS")                                    return OptRitz::IS;
        if(item == "TE")                                    return OptRitz::TE;
    }
    if constexpr(std::is_same_v<T, SVDLibrary>) {
        if(item == "EIGEN")                                 return SVDLibrary::EIGEN;
        if(item == "LAPACKE")                               return SVDLibrary::LAPACKE;
        if(item == "RSVD")                                  return SVDLibrary::RSVD;
    }
    if constexpr(std::is_same_v<T, UpdatePolicy>) {
        auto policy = UpdatePolicy::NEVER;
        if(item.find("WARMUP")     != std::string_view::npos) policy |= UpdatePolicy::WARMUP;
        if(item.find("ITERATION")  != std::string_view::npos) policy |= UpdatePolicy::HALFSWEEP;
        if(item.find("FULLSWEEP")  != std::string_view::npos) policy |= UpdatePolicy::FULLSWEEP;
        if(item.find("TRUNCATED")  != std::string_view::npos) policy |= UpdatePolicy::TRUNCATED;
        if(item.find("SATURATED")  != std::string_view::npos) policy |= UpdatePolicy::SAT_ALGO;
        if(item.find("STUCK")      != std::string_view::npos) policy |= UpdatePolicy::STK_ALGO;
        return policy;
    }
    if constexpr(std::is_same_v<T, SaturationPolicy>) {
        auto policy = SaturationPolicy::val;
        if(item.find("val")  != std::string_view::npos) policy |= SaturationPolicy::val;
        if(item.find("avg")  != std::string_view::npos) policy |= SaturationPolicy::avg;
        if(item.find("med")  != std::string_view::npos) policy |= SaturationPolicy::med;
        if(item.find("mov")  != std::string_view::npos) policy |= SaturationPolicy::mov;
        if(item.find("min")  != std::string_view::npos) policy |= SaturationPolicy::min;
        if(item.find("max")  != std::string_view::npos) policy |= SaturationPolicy::max;
        if(item.find("mid")  != std::string_view::npos) policy |= SaturationPolicy::mid;
        if(item.find("dif")  != std::string_view::npos) policy |= SaturationPolicy::dif;
        if(item.find("log")  != std::string_view::npos) policy |= SaturationPolicy::log;
        return policy;
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
    if constexpr(std::is_same_v<T, ProjectionPolicy>) {
        auto policy = ProjectionPolicy::NEVER;
        if(item.find("INIT")        != std::string_view::npos) policy |= ProjectionPolicy::INIT;
        if(item.find("WARMUP")      != std::string_view::npos) policy |= ProjectionPolicy::WARMUP;
        if(item.find("STUCK")       != std::string_view::npos) policy |= ProjectionPolicy::STUCK;
        if(item.find("ITER")        != std::string_view::npos) policy |= ProjectionPolicy::ITER;
        if(item.find("CONVERGED")   != std::string_view::npos) policy |= ProjectionPolicy::CONVERGED;
        if(item.find("FINISHED")    != std::string_view::npos) policy |= ProjectionPolicy::FINISHED;
        if(item.find("FORCE")       != std::string_view::npos) policy |= ProjectionPolicy::FORCE;
        if(item.find("DEFAULT")     != std::string_view::npos) policy |= ProjectionPolicy::DEFAULT;
        return policy;
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
    if constexpr(std::is_same_v<T, MergeEvent>) {
        if(item == "MOVE")                                  return MergeEvent::MOVE;
        if(item == "NORM")                                  return MergeEvent::NORM;
        if(item == "SWAP")                                  return MergeEvent::SWAP;
        if(item == "GATE")                                  return MergeEvent::GATE;
        if(item == "OPT")                                   return MergeEvent::OPT;
    }
    if constexpr(std::is_same_v<T, ResetReason>) {
        if(item == "INIT")                                  return ResetReason::INIT;
        if(item == "FIND_WINDOW")                           return ResetReason::FIND_WINDOW;
        if(item == "SATURATED")                             return ResetReason::SATURATED;
        if(item == "NEW_STATE")                             return ResetReason::NEW_STATE;
        if(item == "BOND_UPDATE")                           return ResetReason::BOND_UPDATE;
    }
    if constexpr(std::is_same_v<T, EnvExpandMode>) {
        auto policy = EnvExpandMode::NONE;
        if(item.find("H1")       != std::string_view::npos) policy |= EnvExpandMode::H1;
        if(item.find("H2")       != std::string_view::npos) policy |= EnvExpandMode::H2;
        if(item.find("SSITE")    != std::string_view::npos) policy |= EnvExpandMode::SSITE;
        if(item.find("NSITE")    != std::string_view::npos) policy |= EnvExpandMode::NSITE;
        if(item.find("BACKWARD") != std::string_view::npos) policy |= EnvExpandMode::BACKWARD;
        if(item.find("FORWARD")  != std::string_view::npos) policy |= EnvExpandMode::FORWARD;
        if(item.find("DEFAULT")  != std::string_view::npos) policy |= EnvExpandMode::DEFAULT;
        return policy;
    }
    if constexpr(std::is_same_v<T, NormPolicy>) {
        if(item == "ALWAYS")                                return NormPolicy::ALWAYS;
        if(item == "IFNEEDED")                              return NormPolicy::IFNEEDED;
    }
    if constexpr(std::is_same_v<T, LogPolicy>) {
        if(item == "SILENT")                                return LogPolicy::SILENT;
        if(item == "DEBUG")                                 return LogPolicy::DEBUG;
        if(item == "VERBOSE")                               return LogPolicy::VERBOSE;
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
        if(item.find("RBDS") != std::string_view::npos) policy |= StoragePolicy::RBDS;
        if(item.find("RTES") != std::string_view::npos) policy |= StoragePolicy::RTES;
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
    if constexpr(std::is_same_v<T, ResumePolicy>) {
        if(item == "IF_MAX_ITERS")                          return  ResumePolicy::IF_MAX_ITERS;
        if(item == "IF_SATURATED")                          return  ResumePolicy::IF_SATURATED;
        if(item == "IF_UNSUCCESSFUL")                       return  ResumePolicy::IF_UNSUCCESSFUL;
        if(item == "IF_SUCCESSFUL")                         return  ResumePolicy::IF_SUCCESSFUL;
        if(item == "ALWAYS")                                return  ResumePolicy::ALWAYS;
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

    if constexpr(std::is_same_v<T,OptAlgo>){
        if(item == "DMRG")                                  return OptAlgo::DMRG;
        if(item == "DMRGX")                                 return OptAlgo::DMRGX;
        if(item == "HYBRID_DMRGX")                          return OptAlgo::HYBRID_DMRGX;
        if(item == "XDMRG")                                 return OptAlgo::XDMRG;
        if(item == "GDMRG")                                 return OptAlgo::GDMRG;
    }
    if constexpr(std::is_same_v<T,OptSolver>){
        if(item == "EIG")                                   return OptSolver::EIG;
        if(item == "EIGS")                                  return OptSolver::EIGS;
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
