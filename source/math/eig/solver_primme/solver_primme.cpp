#include "../log.h"
#include "../matvec/matvec_dense.h"
#include "../matvec/matvec_mpo.h"
#include "../matvec/matvec_mpos.h"
#include "../matvec/matvec_sparse.h"
#include "../matvec/matvec_zero.h"
#include "../sfinae.h"
#include "../solver.h"
#include "debug/exceptions.h"
#include "tid/tid.h"
#include <math/tenx.h>
#include <primme/primme.h>

inline primme_target stringToTarget(std::string_view ritzstring) {
    if(ritzstring == "LA") return primme_target::primme_largest;
    if(ritzstring == "SA") return primme_target::primme_smallest;
    if(ritzstring == "LM") return primme_target::primme_largest_abs;
    if(ritzstring == "SM") throw except::runtime_error("Primme Ritz SM (smallest magnitude) not implemented");
    if(ritzstring == "LR") return primme_target::primme_largest;
    if(ritzstring == "SR") return primme_target::primme_smallest;
    if(ritzstring == "LI") throw except::runtime_error("Primme Ritz LI (largest imaginary) not implemented");
    if(ritzstring == "SI") throw except::runtime_error("Primme Ritz SI (smallest imaginary) not implemented");
    if(ritzstring == "BE") throw except::runtime_error("Primme Ritz BE (both ends) not implemented");
    if(ritzstring == "primme_smallest") return primme_target::primme_smallest;
    if(ritzstring == "primme_largest") return primme_target::primme_largest;
    if(ritzstring == "primme_closest_geq") return primme_target::primme_closest_geq;
    if(ritzstring == "primme_closest_leq") return primme_target::primme_closest_leq;
    if(ritzstring == "primme_closest_abs") return primme_target::primme_closest_abs;
    if(ritzstring == "primme_largest_abs") return primme_target::primme_largest_abs;
    throw except::runtime_error("Wrong ritz string: " + std::string(ritzstring));
}

inline std::string_view TargetToString(primme_target target) {
    if(target == primme_target::primme_smallest) return "SA";
    if(target == primme_target::primme_largest) return "LA";
    if(target == primme_target::primme_closest_geq) return "GQ";
    if(target == primme_target::primme_closest_leq) return "LQ";
    if(target == primme_target::primme_closest_abs) return "SM";
    if(target == primme_target::primme_largest_abs) return "LM";
    throw std::logic_error("Could not match primme_target");
}

inline primme_target RitzToTarget(eig::Ritz ritz) {
    if(ritz == eig::Ritz::LA) return primme_target::primme_largest;
    if(ritz == eig::Ritz::SA) return primme_target::primme_smallest;
    if(ritz == eig::Ritz::LM) return primme_target::primme_largest_abs;
    if(ritz == eig::Ritz::SM) return primme_target::primme_closest_abs;
    if(ritz == eig::Ritz::LR) return primme_target::primme_largest;
    if(ritz == eig::Ritz::SR) return primme_target::primme_smallest;
    if(ritz == eig::Ritz::LI) throw except::runtime_error("Primme Ritz LI (largest imaginary) not implemented");
    if(ritz == eig::Ritz::SI) throw except::runtime_error("Primme Ritz SI (smallest imaginary) not implemented");
    if(ritz == eig::Ritz::BE) throw except::runtime_error("Primme Ritz BE (both ends) not implemented");
    if(ritz == eig::Ritz::primme_smallest) return primme_target::primme_smallest;
    if(ritz == eig::Ritz::primme_largest) return primme_target::primme_largest;
    if(ritz == eig::Ritz::primme_closest_geq) return primme_target::primme_closest_geq;
    if(ritz == eig::Ritz::primme_closest_leq) return primme_target::primme_closest_leq;
    if(ritz == eig::Ritz::primme_closest_abs) return primme_target::primme_closest_abs;
    if(ritz == eig::Ritz::primme_largest_abs) return primme_target::primme_largest_abs;
    throw except::runtime_error("Wrong ritz enum");
}

inline primme_preset_method stringToMethod(std::optional<std::string> methodstring) {
    if(not methodstring.has_value()) return PRIMME_DEFAULT_MIN_MATVECS;
    if(methodstring.value() == "PRIMME_DEFAULT_METHOD") return primme_preset_method::PRIMME_DEFAULT_METHOD;
    if(methodstring.value() == "PRIMME_DYNAMIC") return primme_preset_method::PRIMME_DYNAMIC;
    if(methodstring.value() == "PRIMME_DEFAULT_MIN_TIME") return primme_preset_method::PRIMME_DEFAULT_MIN_TIME;
    if(methodstring.value() == "PRIMME_DEFAULT_MIN_MATVECS") return primme_preset_method::PRIMME_DEFAULT_MIN_MATVECS;
    if(methodstring.value() == "PRIMME_Arnoldi") return primme_preset_method::PRIMME_Arnoldi;
    if(methodstring.value() == "PRIMME_GD") return primme_preset_method::PRIMME_GD;
    if(methodstring.value() == "PRIMME_GD_plusK") return primme_preset_method::PRIMME_GD_plusK;
    if(methodstring.value() == "PRIMME_GD_Olsen_plusK") return primme_preset_method::PRIMME_GD_Olsen_plusK;
    if(methodstring.value() == "PRIMME_JD_Olsen_plusK") return primme_preset_method::PRIMME_JD_Olsen_plusK;
    if(methodstring.value() == "PRIMME_RQI") return primme_preset_method::PRIMME_RQI;
    if(methodstring.value() == "PRIMME_JDQR") return primme_preset_method::PRIMME_JDQR;
    if(methodstring.value() == "PRIMME_JDQMR") return primme_preset_method::PRIMME_JDQMR;
    if(methodstring.value() == "PRIMME_JDQMR_ETol") return primme_preset_method::PRIMME_JDQMR_ETol;
    if(methodstring.value() == "PRIMME_STEEPEST_DESCENT") return primme_preset_method::PRIMME_STEEPEST_DESCENT;
    if(methodstring.value() == "PRIMME_LOBPCG_OrthoBasis") return primme_preset_method::PRIMME_LOBPCG_OrthoBasis;
    if(methodstring.value() == "PRIMME_LOBPCG_OrthoBasis_Window") return primme_preset_method::PRIMME_LOBPCG_OrthoBasis_Window;
    throw except::runtime_error("Wrong method string: " + methodstring.value());
}

inline primme_preset_method MethodAdapter(std::optional<eig::PrimmeMethod> method) {
    if(not method.has_value()) return primme_preset_method::PRIMME_DYNAMIC;
    if(method.value() == eig::PrimmeMethod::PRIMME_DEFAULT_METHOD) return primme_preset_method::PRIMME_DEFAULT_METHOD;
    if(method.value() == eig::PrimmeMethod::PRIMME_DYNAMIC) return primme_preset_method::PRIMME_DYNAMIC;
    if(method.value() == eig::PrimmeMethod::PRIMME_DEFAULT_MIN_TIME) return primme_preset_method::PRIMME_DEFAULT_MIN_TIME;
    if(method.value() == eig::PrimmeMethod::PRIMME_DEFAULT_MIN_MATVECS) return primme_preset_method::PRIMME_DEFAULT_MIN_MATVECS;
    if(method.value() == eig::PrimmeMethod::PRIMME_Arnoldi) return primme_preset_method::PRIMME_Arnoldi;
    if(method.value() == eig::PrimmeMethod::PRIMME_GD) return primme_preset_method::PRIMME_GD;
    if(method.value() == eig::PrimmeMethod::PRIMME_GD_plusK) return primme_preset_method::PRIMME_GD_plusK;
    if(method.value() == eig::PrimmeMethod::PRIMME_GD_Olsen_plusK) return primme_preset_method::PRIMME_GD_Olsen_plusK;
    if(method.value() == eig::PrimmeMethod::PRIMME_JD_Olsen_plusK) return primme_preset_method::PRIMME_JD_Olsen_plusK;
    if(method.value() == eig::PrimmeMethod::PRIMME_RQI) return primme_preset_method::PRIMME_RQI;
    if(method.value() == eig::PrimmeMethod::PRIMME_JDQR) return primme_preset_method::PRIMME_JDQR;
    if(method.value() == eig::PrimmeMethod::PRIMME_JDQMR) return primme_preset_method::PRIMME_JDQMR;
    if(method.value() == eig::PrimmeMethod::PRIMME_JDQMR_ETol) return primme_preset_method::PRIMME_JDQMR_ETol;
    if(method.value() == eig::PrimmeMethod::PRIMME_STEEPEST_DESCENT) return primme_preset_method::PRIMME_STEEPEST_DESCENT;
    if(method.value() == eig::PrimmeMethod::PRIMME_LOBPCG_OrthoBasis) return primme_preset_method::PRIMME_LOBPCG_OrthoBasis;
    if(method.value() == eig::PrimmeMethod::PRIMME_LOBPCG_OrthoBasis_Window) return primme_preset_method::PRIMME_LOBPCG_OrthoBasis_Window;
    return primme_preset_method::PRIMME_DYNAMIC;
}

inline primme_projection stringToProj(std::optional<std::string> projstring) {
    if(not projstring.has_value()) return primme_projection::primme_proj_default;
    if(projstring.value() == "primme_proj_RR") return primme_projection::primme_proj_RR;
    if(projstring.value() == "primme_proj_harmonic") return primme_projection::primme_proj_harmonic;
    if(projstring.value() == "primme_proj_refined") return primme_projection::primme_proj_refined;
    if(projstring.value() == "primme_proj_default") return primme_projection::primme_proj_default;
    throw except::runtime_error("Wrong projection string: " + projstring.value());
}
inline std::string_view projToString(primme_projection proj) {
    switch(proj) {
        case primme_projection::primme_proj_RR: return "primme_proj_RR";
        case primme_projection::primme_proj_harmonic: return "primme_proj_harmonic";
        case primme_projection::primme_proj_refined: return "primme_proj_refined";
        case primme_projection::primme_proj_default: return "primme_proj_default";
        default: return "primme_proj_UNKNOWN";
    }
}

template<typename MatrixProductType>
void eig::solver::MultAx_wrapper(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
    auto t_mv = tid::tic_scope("matvec", tid::level::highest);
    if(primme->matrix == nullptr) throw std::logic_error("primme->matrix == nullptr");
    auto matrix_ptr = static_cast<MatrixProductType *>(primme->matrix);
    matrix_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
}

template<typename MatrixProductType>
void eig::solver::MultOPv_wrapper(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
    auto t_iv = tid::tic_scope("invvec", tid::level::highest);
    if(primme->matrix == nullptr) throw std::logic_error("primme->matrix == nullptr");
    auto matrix_ptr = static_cast<MatrixProductType *>(primme->matrix);
    matrix_ptr->MultOPv(x, ldx, y, ldy, blockSize, primme, ierr);
}

std::string getLogMessage(struct primme_params *primme) {
    if(primme->monitor == nullptr) {
        return fmt::format("iter {:>6} | mv {:>6} | size {} | λ {:19.16f} | time {:9.3e}s | {:8.2e} it/s | {:8.2e} mv/s", primme->stats.numOuterIterations,
                           primme->stats.numMatvecs, primme->n, primme->stats.estimateMinEVal, primme->stats.elapsedTime,
                           primme->stats.numOuterIterations / primme->stats.elapsedTime, primme->stats.numMatvecs / primme->stats.timeMatvec);
    }
    auto &solver     = *static_cast<eig::solver *>(primme->monitor);
    auto &result     = solver.result;
    auto &eigvals    = result.get_eigvals<eig::Form::SYMM>();
    auto  get_eigval = [&](long idx) -> double {
        auto eigvmap = Eigen::Map<const Eigen::VectorXd>(eigvals.data(), eigvals.size());
        // auto basismap = Eigen::Map<const Eigen::VectorXd>(primme->stats.data(), eigvals.size());
        if(solver.config.shift_invert == eig::Shinv::ON) {
            return 1.0 / (eigvmap[idx]) + std::real(solver.config.sigma.value_or(0.0));
        } else {
            switch(primme->target) {
                case primme_target::primme_closest_abs: [[fallthrough]];
                case primme_target::primme_closest_geq: [[fallthrough]];
                case primme_target::primme_closest_leq: return result.meta.last_basis_eval;
                case primme_target::primme_smallest: return primme->stats.estimateMinEVal;
                case primme_target::primme_largest: return primme->stats.estimateMaxEVal;
                case primme_target::primme_largest_abs: return std::max(std::abs(primme->stats.estimateMaxEVal), std::abs(primme->stats.estimateMaxEVal));
                default: return primme->stats.estimateMinEVal;
            }
        }
    };
    std::string msg_diff = eigvals.size() >= 2 ? fmt::format(" | λ1-λ0 {:19.16f}", std::abs(get_eigval(0) - get_eigval(1))) : "";
    return fmt::format("iter {:>6} | mv {:>6}/{}| size {} | λ {:19.16f}{} | rnorm {:8.2e} | time {:9.3e}s | {:8.2e} "
                       "it/s | {:8.2e} mv/s | {} | {}",
                       primme->stats.numOuterIterations, primme->stats.numMatvecs, primme->maxMatvecs, primme->n, get_eigval(0), msg_diff,
                       std::max(primme->stats.maxConvTol, solver.result.meta.last_res_norm), primme->stats.elapsedTime,
                       primme->stats.numOuterIterations / primme->stats.elapsedTime, static_cast<double>(primme->stats.numMatvecs) / primme->stats.timeMatvec,
                       eig::TypeToString(solver.config.type), eig::MethodToString(solver.config.primme_method));
}

void monitorFun([[maybe_unused]] void *basisEvals, [[maybe_unused]] int *basisSize, [[maybe_unused]] int *basisFlags, [[maybe_unused]] int *iblock,
                [[maybe_unused]] int *blockSize, [[maybe_unused]] void *basisNorms, [[maybe_unused]] int *numConverged, [[maybe_unused]] void *lockedEvals,
                [[maybe_unused]] int *numLocked, [[maybe_unused]] int *lockedFlags, [[maybe_unused]] void *lockedNorms, [[maybe_unused]] int *inner_its,
                [[maybe_unused]] void *LSRes, [[maybe_unused]] const char *msg, [[maybe_unused]] double *time, [[maybe_unused]] primme_event *event,
                [[maybe_unused]] struct primme_params *primme, [[maybe_unused]] int *ierr) {
    if(event == nullptr) return;
    std::string eventMessage;
    /* clang-format off */
    switch(*event){
        case primme_event_outer_iteration : eventMessage = "event_outer_iteration"; break;
        case primme_event_inner_iteration : eventMessage = "event_inner_iteration"; break;
        case primme_event_restart         : return;
        case primme_event_reset           : return;
        case primme_event_converged       : return;
        case primme_event_locked          : return;
        case primme_event_message         : return;
        case primme_event_profile         : return;
    }
    /* clang-format on */
    std::string basisMessage = basisSize != nullptr ? fmt::format(" | ncv {:3}", *basisSize) : "";
    std::string nlockMessage = numLocked != nullptr and primme->locking == 1 ? fmt::format(" | nlk {:3}", *numLocked) : "";

    if(primme->monitor != nullptr) {
        auto &solver              = *static_cast<eig::solver *>(primme->monitor);
        auto &config              = solver.config;
        auto &result              = solver.result;
        auto  level               = spdlog::level::trace;
        auto  nmvs_since_last_log = std::abs(primme->stats.numMatvecs - result.meta.last_log_iter);
        auto  time_since_last_log = std::abs(primme->stats.elapsedTime - result.meta.last_log_time);
        if(*event == primme_event_outer_iteration or *event == primme_event_inner_iteration) {
            if(config.logTime and config.logTime.value() <= time_since_last_log) level = eig::log->level();
            if(config.logIter and config.logIter.value() <= nmvs_since_last_log) level = eig::log->level();
        } else if(*event == primme_event_converged)
            level = spdlog::level::debug;

        // Terminate if it's taking too long
        if(config.maxTime.has_value() and primme->stats.elapsedTime > config.maxTime.value()) {
            if(primme->maxMatvecs > 0 and primme->maxOuterIterations > 0) eig::log->warn("primme: max time has been exeeded: {:.2f}", config.maxTime.value());
            primme->maxMatvecs         = 0;
            primme->maxOuterIterations = 0;

            // *ierr = 60;
        }
        if(eig::log->level() <= level) {
            // double mineval = 0;
            if(basisNorms != nullptr and basisSize != nullptr) {
                auto norms                  = Eigen::Map<const Eigen::VectorXd>(static_cast<double *>(basisNorms), *basisSize);
                auto evals                  = Eigen::Map<const Eigen::VectorXd>(static_cast<double *>(basisEvals), *basisSize);
                result.meta.last_res_norm   = norms.maxCoeff();
                result.meta.last_basis_eval = evals[0];
            }
            eig::log->log(level, "{}{}{} | {}", getLogMessage(primme), basisMessage, nlockMessage, eventMessage);
            result.meta.last_log_time = primme->stats.elapsedTime;
            result.meta.last_log_iter = primme->stats.numMatvecs;
        }
    } else {
        eig::log->trace("{} | {}", getLogMessage(primme), eventMessage);
    }
}

template<typename MatrixProductType>
int eig::solver::eigs_primme(MatrixProductType &matrix) {
    using Scalar  = typename MatrixProductType::Scalar;
    auto t_primme = tid::tic_scope("primme", tid::level::higher);
    auto t_prep   = tid::tic_scope("prep");
    if constexpr(MatrixProductType::can_shift) {
        if(config.sigma) {
            auto t_shift = tid::tic_scope("shift", tid::level::highest);
            eig::log->debug("Setting shift with sigma = {}", std::real(config.sigma.value()));
            matrix.set_shift(config.sigma.value());
            if constexpr(MatrixProductType::can_shift_invert) {
                if(config.shift_invert == Shinv::ON) {
                    eig::log->debug("Enabling shift-invert mode");
                    auto t_inv = tid::tic_scope("invert", tid::level::highest);
                    matrix.FactorOP();
                }
            } else if(config.shift_invert == Shinv::ON)
                throw std::runtime_error("Tried to shift-invert an incompatible matrix");
        }
    } else if(config.sigma) {
        throw std::runtime_error("Tried to apply shift on an incompatible matrix");
    }

    primme_params primme;

    /* Set default values in PRIMME configuration struct */
    primme_initialize(&primme);

    /* Set problem matrix
     * We need to bind the matrix-vector product to a void pointer, then call it from a wrapper.
     * This is a delicate step.
     */
    primme.matrix       = &matrix;
    primme.matrixMatvec = eig::solver::MultAx_wrapper<MatrixProductType>;

    if(matrix.isReadyFactorOp()) { primme.matrixMatvec = eig::solver::MultOPv_wrapper<MatrixProductType>; }

    primme.monitorFun = monitorFun; // Set the function which prints log output
    primme.monitor    = this;       // Make eig objects visible from within monitorFun

    // Set custom convergence test if given
    if(config.primme_convTestFun) {
        primme.convTestFun = config.primme_convTestFun.value();
        primme.convtest    = this;
    }
    // Set custom preconditioner if given
    if(config.primme_preconditioner) {
        primme.applyPreconditioner           = config.primme_preconditioner.value();
        primme.preconditioner                = this;
        primme.correctionParams.precondition = true;
    }

    /* Set problem parameters */
    // https://www.cs.wm.edu/~andreas/software/doc/appendix.html#c.primme_params.eps
    primme.n                           = matrix.rows();                                                                /*!< set problem dimension */
    primme.numEvals                    = config.maxNev.value_or(1);                                                    /*!<  Number of desired eigenpairs */
    primme.eps                         = std::max(config.tol.value_or(1e-12), std::numeric_limits<double>::epsilon()); /*!< 1e-12 is good, see link above. */
    primme.target                      = RitzToTarget(config.ritz.value_or(Ritz::primme_smallest));
    primme.projectionParams.projection = stringToProj(config.primme_projection.value_or("primme_proj_default"));
    primme.maxOuterIterations          = config.maxIter.value_or(primme.maxOuterIterations);
    primme.correctionParams.maxInnerIterations = -1; // Up to the remaining matvecs
    primme.maxMatvecs                          = config.maxIter.value_or(primme.maxMatvecs);
    primme.maxBlockSize                        = config.primme_maxBlockSize.value_or(1);
    primme.maxBasisSize                        = getBasisSize(primme.n, primme.numEvals, config.maxNcv);
    primme.minRestartSize                      = config.primme_minRestartSize.value_or(primme.maxBasisSize / 2);
    // Make sure the basis is bigger than minRestartSize + maxPrevRetain, where the latter will be set to max(2,maxBlockSize) in primme_set_method
    primme.maxBasisSize = std::clamp(primme.maxBasisSize, primme.minRestartSize + std::max(2, primme.maxBlockSize + 1), primme.n);

    // Shifts
    switch(primme.target) {
        case primme_target::primme_largest:
        case primme_target::primme_smallest: {
            primme.numTargetShifts = 0;
            break;
        }
        case primme_target::primme_largest_abs: {
            // According to the manual, this is required when the target is primme_largest_abs
            config.primme_targetShifts = {0.0};
            primme.numTargetShifts     = 1;
            primme.targetShifts        = config.primme_targetShifts.data();
            break;
        }
        case primme_target::primme_closest_geq:
        case primme_target::primme_closest_leq:
        case primme_target::primme_closest_abs: {
            if(not config.primme_targetShifts.empty()) {
                eig::log->debug("Setting target shifts: {:.8f}", fmt::join(config.primme_targetShifts, ", "));
                primme.numTargetShifts = safe_cast<int>(config.primme_targetShifts.size());
                primme.targetShifts    = config.primme_targetShifts.data();
            } else if(config.sigma and not matrix.isReadyShift()) {
                // We handle shifts by applying them directly on the matrix if possible. Else here:
                eig::log->debug("Setting target shift: {:.8f}", std::real(config.sigma.value()));
                config.primme_targetShifts = {std::real(config.sigma.value())};
                primme.numTargetShifts     = safe_cast<int>(config.primme_targetShifts.size());
                primme.targetShifts        = config.primme_targetShifts.data();
            } else
                throw std::runtime_error("primme_target:primme_closest_???: no target shift given");

            break;
        }
    }

    /* Set method to solve the problem
     * DYNAMIC uses a runtime heuristic to choose the fastest method between
     * PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
     * set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */
    primme_set_method(MethodAdapter(config.primme_method), &primme);

    eig::log->debug("numEvals {} | maxMatvecs {} | eps {:.2e} | maxBasisSize {} | minRestartSize {} | maxBlockSize {} | maxPrevRetain {} | {}", primme.numEvals,
                    primme.maxMatvecs, primme.eps, primme.maxBasisSize, primme.minRestartSize, primme.maxBlockSize, primme.restartingParams.maxPrevRetain,
                    projToString(primme.projectionParams.projection));
    // Allocate space
    auto &eigvals = result.get_eigvals<eig::Form::SYMM>();
    auto &eigvecs = result.get_eigvecs<Scalar, eig::Form::SYMM>();
    eigvals.resize(static_cast<size_t>(primme.numEvals));
    eigvecs.resize(static_cast<size_t>(primme.numEvals) * static_cast<size_t>(primme.n));
    result.meta.residual_norms.resize(static_cast<size_t>(primme.numEvals));

    {
        // Copy initial guess
        primme.initSize = 0;
        for(const auto &ig : config.initial_guess) {
            if(ig.idx < primme.numEvals and ig.idx >= 0) {
                auto source = static_cast<Scalar *>(ig.ptr);
                auto target = eigvecs.data() + static_cast<size_t>(ig.idx * primme.n);
                std::copy_n(source, primme.n, target);
                primme.initSize = safe_cast<int>(ig.idx) + 1;
            }
        }
    }

    t_prep.toc();
    /* Call primme  */
    /* clang-format off */
    int info = 0;
    if constexpr(std::is_same_v<Scalar, real>) {
        try {
            auto t_dprimme = tid::tic_scope("dprimme");
            info = dprimme(eigvals.data(), eigvecs.data(), result.meta.residual_norms.data(), &primme);
        } catch(const std::exception &ex) { eig::log->error("dprimme has exited with error: info {} | msg: {}", info, ex.what()); } catch(...) {
            eig::log->error("dprimme has exited with error: info {}", info);
        }
    } else if constexpr(std::is_same_v<typename MatrixProductType::Scalar, cplx>) {
        try {
            auto t_zprimme = tid::tic_scope("zprimme");
            info = zprimme(eigvals.data(), eigvecs.data(), result.meta.residual_norms.data(), &primme);
        } catch(const std::exception &ex) { eig::log->error("zprimme has exited with error: info {} | msg: {}", info, ex.what()); } catch(...) {
            eig::log->error("zprimme has exited with error: info {}", info);
        }
    } else {
        throw std::runtime_error("Primme: type not implemented {}", eig::sfinae::type_name<typename MatrixProductType::Scalar>());
    }
    /* clang-format on */

    if(info == 0) {
        switch(primme.dynamicMethodSwitch) {
            case -1: eig::log->debug("Recommended method for next run: DEFAULT_MIN_MATVECS"); break;
            case -2: eig::log->debug("Recommended method for next run: DEFAULT_MIN_TIME"); break;
            case -3: eig::log->debug("Recommended method for next run: DYNAMIC (close call)"); break;
            default: break;
        }
    } else {
        switch(info) {
            case -1: eig::log->error("PRIMME_UNEXPECTED_FAILURE (exit {}): set printLevel > 0 to see the call stack", info); break;
            case -2: eig::log->error("PRIMME_MALLOC_FAILURE (exit {}): eith<er CPU or GPU", info); break;
            case -3: eig::log->debug("{}", getLogMessage(&primme)); break;
            case -4: eig::log->error("PRIMME_ARGUMENT_IS_NULL (exit {})", info); break;
            case -5: eig::log->error("PRIMME_INVALID_ARG n < 0 or nLocal < 0 or nLocal > n (exit {})", info); break;
            case -6: eig::log->error("PRIMME_INVALID_ARG numProcs < 1 (exit {})", info); break;
            case -25:
                eig::log->error("PRIMME_INVALID_ARG  maxPrevRetain({}) + minRestartSize({}) >= maxBasisSize({}), and n({}) > maxBasisSize({}) (exit {})",
                                primme.restartingParams.maxPrevRetain, primme.minRestartSize, primme.maxBasisSize, primme.n, primme.maxBasisSize, info);
                break;
            case -40: eig::log->error("PRIMME_LAPACK_FAILURE (exit {})", info); break;
            case -41: eig::log->error("PRIMME_USER_FAILURE (exit {})", info); break;
            case -42: eig::log->error("PRIMME_ORTHO_CONST_FAILURE (exit {})", info); break;
            case -43: eig::log->error("PRIMME_PARALLEL_FAILURE (exit {})", info); break;
            case -44: eig::log->error("PRIMME_FUNCTION_UNAVAILABLE (exit {})", info); break;
            case -60: eig::log->debug("PRIMME_MAX_TIME_EXCEEDED (exit {})", info); break;
            default: eig::log->error("Unknown error code: {}.\n Go to http://www.cs.wm.edu/~andreas/software/doc/appendix.html", info);
        }
    }
    result.meta.eigvecsR_found = true; // We can use partial results
    result.meta.eigvals_found  = true; // We can use partial results
    result.meta.rows           = primme.n;
    result.meta.cols           = primme.numEvals;
    result.meta.nev            = primme.numEvals;
    result.meta.nev_converged  = primme.initSize; // Holds the number of converged eigenpairs on exit
    result.meta.ncv            = primme.maxBasisSize;
    result.meta.tol            = primme.eps;
    result.meta.iter           = primme.stats.numOuterIterations;
    result.meta.num_mv         = primme.stats.numMatvecs;
    result.meta.num_pc         = primme.stats.numPreconds;
    result.meta.num_op         = matrix.num_op;
    result.meta.time_mv        = primme.stats.timeMatvec;
    result.meta.time_pc        = primme.stats.timePrecond;
    result.meta.time_op        = matrix.t_multOPv->get_time();
    result.meta.n              = primme.n;
    result.meta.tag            = config.tag;
    result.meta.ritz           = TargetToString(primme.target);
    result.meta.form           = config.form.value();
    result.meta.type           = config.type.value();
    result.meta.time_total     = primme.stats.elapsedTime;
    result.meta.last_res_norm  = *std::max_element(result.meta.residual_norms.begin(), result.meta.residual_norms.end());
    primme_free(&primme);
    return info;
}

template int eig::solver::eigs_primme(MatVecMPO<real> &matrix);
template int eig::solver::eigs_primme(MatVecMPO<cplx> &matrix);
template int eig::solver::eigs_primme(MatVecMPOS<real> &matrix);
template int eig::solver::eigs_primme(MatVecMPOS<cplx> &matrix);
template int eig::solver::eigs_primme(MatVecZero<real> &matrix);
template int eig::solver::eigs_primme(MatVecZero<cplx> &matrix);
template int eig::solver::eigs_primme(MatVecDense<real> &matrix);
template int eig::solver::eigs_primme(MatVecDense<cplx> &matrix);
template int eig::solver::eigs_primme(MatVecSparse<real> &matrix);
template int eig::solver::eigs_primme(MatVecSparse<cplx> &matrix);
