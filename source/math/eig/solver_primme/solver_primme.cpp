#include "../log.h"
#include "../matvec/matvec_dense.h"
#include "../matvec/matvec_mpo.h"
#include "../matvec/matvec_sparse.h"
#include "../sfinae.h"
#include "../solver.h"
#include <primme/primme.h>
#include <tid/tid.h>

inline primme_target stringToTarget(std::string_view ritzstring) {
    if(ritzstring == "LA") return primme_target::primme_largest;
    if(ritzstring == "SA") return primme_target::primme_smallest;
    if(ritzstring == "LM") return primme_target::primme_largest_abs;
    if(ritzstring == "SM") throw std::runtime_error("Primme Ritz SM (smallest magnitude) not implemented");
    if(ritzstring == "LR") return primme_target::primme_largest;
    if(ritzstring == "SR") return primme_target::primme_smallest;
    if(ritzstring == "LI") throw std::runtime_error("Primme Ritz LI (largest imaginary) not implemented");
    if(ritzstring == "SI") throw std::runtime_error("Primme Ritz SI (smallest imaginary) not implemented");
    if(ritzstring == "BE") throw std::runtime_error("Primme Ritz BE (both ends) not implemented");
    if(ritzstring == "primme_smallest") return primme_target::primme_smallest;
    if(ritzstring == "primme_largest") return primme_target::primme_largest;
    if(ritzstring == "primme_closest_geq") return primme_target::primme_closest_geq;
    if(ritzstring == "primme_closest_leq") return primme_target::primme_closest_leq;
    if(ritzstring == "primme_closest_abs") return primme_target::primme_closest_abs;
    if(ritzstring == "primme_largest_abs") return primme_target::primme_largest_abs;
    throw std::runtime_error("Wrong ritz string: " + std::string(ritzstring));
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
    if(ritz == eig::Ritz::LI) throw std::runtime_error("Primme Ritz LI (largest imaginary) not implemented");
    if(ritz == eig::Ritz::SI) throw std::runtime_error("Primme Ritz SI (smallest imaginary) not implemented");
    if(ritz == eig::Ritz::BE) throw std::runtime_error("Primme Ritz BE (both ends) not implemented");
    if(ritz == eig::Ritz::primme_smallest) return primme_target::primme_smallest;
    if(ritz == eig::Ritz::primme_largest) return primme_target::primme_largest;
    if(ritz == eig::Ritz::primme_closest_geq) return primme_target::primme_closest_geq;
    if(ritz == eig::Ritz::primme_closest_leq) return primme_target::primme_closest_leq;
    if(ritz == eig::Ritz::primme_closest_abs) return primme_target::primme_closest_abs;
    if(ritz == eig::Ritz::primme_largest_abs) return primme_target::primme_largest_abs;
    throw std::runtime_error("Wrong ritz enum");
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
    throw std::runtime_error("Wrong method string: " + methodstring.value());
}

inline primme_preset_method MethodAdapter(std::optional<eig::PrimmeMethod> method) {
    if(not method.has_value()) return primme_preset_method::PRIMME_DEFAULT_MIN_MATVECS;
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
    return primme_preset_method::PRIMME_DEFAULT_MIN_MATVECS;
}

inline primme_projection stringToProj(std::optional<std::string> projstring) {
    if(not projstring.has_value()) return primme_projection::primme_proj_default;
    if(projstring.value() == "primme_proj_RR") return primme_projection::primme_proj_RR;
    if(projstring.value() == "primme_proj_harmonic") return primme_projection::primme_proj_harmonic;
    if(projstring.value() == "primme_proj_refined") return primme_projection::primme_proj_refined;
    if(projstring.value() == "primme_proj_default") return primme_projection::primme_proj_default;
    throw std::runtime_error("Wrong projection string: " + projstring.value());
}

template<typename MatrixProductType>
void eig::solver::MultAx_wrapper(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
    if(primme->matrix == nullptr) throw std::logic_error("primme->matrix == nullptr");
    auto matrix_ptr = static_cast<MatrixProductType *>(primme->matrix);
    return matrix_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
}

template<typename MatrixProductType>
void eig::solver::MultOPv_wrapper(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
    if(primme->matrix == nullptr) throw std::logic_error("primme->matrix == nullptr");
    auto matrix_ptr = static_cast<MatrixProductType *>(primme->matrix);
    matrix_ptr->MultOPv(x, ldx, y, ldy, blockSize, primme, ierr);
}

std::string getLogMessage(struct primme_params *primme) {
    if(primme->monitor == nullptr) {
        return fmt::format(FMT_STRING("mv {:<5} | iter {:<4} | size {} | f {:20.16f} | time {:8.2f} s | {:8.2e} s/it | {:8.2e} s/mv"), primme->stats.numMatvecs,
                           primme->stats.numOuterIterations, primme->n, primme->stats.estimateMinEVal, primme->stats.elapsedTime,
                           primme->stats.elapsedTime / primme->stats.numOuterIterations, primme->stats.timeMatvec / primme->stats.numMatvecs);
    }
    auto       &solver     = *static_cast<eig::solver *>(primme->monitor);
    auto       &result     = solver.result;
    auto       &eigvals    = result.get_eigvals<eig::Form::SYMM>();
    std::string msg_diff   = eigvals.size() >= 2 ? fmt::format(" | f1-f0 {:20.16f}", std::abs(eigvals[0] - eigvals[1])) : "";
    std::string msg_grad   = primme->convTestFun != nullptr ? fmt::format(" | ∇fᵐᵃˣ {:8.2e}", result.meta.last_grad_max) : "";
    auto        res        = std::numeric_limits<double>::quiet_NaN();
    auto        max_res_it = std::max_element(result.meta.residual_norms.begin(), result.meta.residual_norms.end());
    if(max_res_it != result.meta.residual_norms.end()) {
        res                       = *max_res_it;
        result.meta.last_res_norm = res;
    }
    return fmt::format(FMT_STRING("mv {:<5} | iter {:<4} | size {} | res {:8.2e} | f {:20.16f}{}{} | time {:8.2f} s | {:8.2e} s/it | {:8.2e} s/mv"),
                       primme->stats.numMatvecs, primme->stats.numOuterIterations, primme->n, res, primme->stats.estimateMinEVal, msg_diff, msg_grad,
                       primme->stats.elapsedTime, primme->stats.elapsedTime / primme->stats.numOuterIterations,
                       primme->stats.timeMatvec / primme->stats.numMatvecs);
}

void monitorFun([[maybe_unused]] void *basisEvals, [[maybe_unused]] int *basisSize, [[maybe_unused]] int *basisFlags, [[maybe_unused]] int *iblock,
                [[maybe_unused]] int *blockSize, [[maybe_unused]] void *basisNorms, [[maybe_unused]] int *numConverged, [[maybe_unused]] void *lockedEvals,
                [[maybe_unused]] int *numLocked, [[maybe_unused]] int *lockedFlags, [[maybe_unused]] void *lockedNorms, [[maybe_unused]] int *inner_its,
                [[maybe_unused]] void *LSRes, [[maybe_unused]] const char *msg, [[maybe_unused]] double *time, [[maybe_unused]] primme_event *event,
                [[maybe_unused]] struct primme_params *primme, [[maybe_unused]] int *ierr) {
    if(event == nullptr) return;
    if(*event == primme_event_inner_iteration) return; // No need to log that often

    std::string eventMessage;
    /* clang-format off */
    if(*event == primme_event_outer_iteration)      eventMessage = "event_outer_iteration";
    else if(*event == primme_event_inner_iteration) eventMessage = "event_inner_iteration";
    else if(*event == primme_event_restart)         eventMessage = "event_restart";
    else if(*event == primme_event_reset)           eventMessage = "event_reset";
    else if(*event == primme_event_converged)       eventMessage = "event_converged";
    else if(*event == primme_event_locked)          eventMessage = "event_locked";
    else if(*event == primme_event_message)         eventMessage = "event_message";
    else if(*event == primme_event_profile)         eventMessage = "event_profile";
    /* clang-format on */
    std::string basisMessage = basisSize != nullptr ? fmt::format(" | ncv {:3}", *basisSize) : "";

    if(primme->monitor != nullptr) {
        auto &solver              = *static_cast<eig::solver *>(primme->monitor);
        auto &config              = solver.config;
        auto &result              = solver.result;
        auto  level               = spdlog::level::trace;
        auto  iter_since_last_log = std::abs(primme->stats.numOuterIterations - result.meta.last_log_iter);
        auto  time_since_last_log = std::abs(primme->stats.elapsedTime - result.meta.last_log_time);
        if(*event == primme_event_outer_iteration) {
            if(config.logTime and config.logTime.value() <= time_since_last_log) level = eig::log->level();
            if(config.logIter and config.logIter.value() <= iter_since_last_log) level = eig::log->level();
        } else if(*event == primme_event_converged)
            level = spdlog::level::debug;

        // Terminate if it's taking too long
        if(config.maxTime.has_value() and primme->stats.elapsedTime > config.maxTime.value()) {
            eig::log->warn("primme: max time has been exeeded: {:.2f}", config.maxTime.value());
            primme->maxMatvecs = 0;
        }
        if(eig::log->level() <= level) {
            eig::log->log(level, FMT_STRING("{}{} | {}"), getLogMessage(primme), basisMessage, eventMessage);
            result.meta.last_log_time = primme->stats.elapsedTime;
            result.meta.last_log_iter = primme->stats.numOuterIterations;
        }
    } else {
        eig::log->trace(FMT_STRING("{} | {}"), getLogMessage(primme), eventMessage);
    }
}

template<typename MatrixProductType>
int eig::solver::eigs_primme(MatrixProductType &matrix) {
    using Scalar = typename MatrixProductType::Scalar;

    tid::ur t_tot, t_pre;
    auto    t_tot_token = t_tot.tic_token();
    auto    t_pre_token = t_pre.tic_token();
    if constexpr(MatrixProductType::can_shift) {
        if(config.sigma) {
            eig::log->trace("Setting shift with sigma = {}", std::real(config.sigma.value()));
            matrix.set_shift(config.sigma.value());
            if constexpr(MatrixProductType::can_shift_invert) {
                if(config.shift_invert == Shinv::ON) {
                    eig::log->trace("Enabling shift-invert mode");
                    matrix.FactorOP();
                }
            } else if(config.shift_invert == Shinv::ON)
                throw std::runtime_error("Tried to shift-invert an incompatible matrix");
        }
    } else if(config.sigma) {
        throw std::runtime_error("Tried to apply shift on an incompatible matrix");
    }

    if constexpr(MatrixProductType::can_compress) {
        if(config.compress and config.compress.value()) matrix.compress();
    }
    t_pre_token.toc();

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
    primme.n = matrix.rows();                                                    /* set problem dimension */
    if(config.maxNev) primme.numEvals = static_cast<int>(config.maxNev.value()); /* Number of wanted eigenpairs */
    if(config.tol) primme.eps = config.tol.value();                              /* ||r|| <= eps * ||matrix|| */
    if(config.ritz) primme.target = RitzToTarget(config.ritz.value());
    if(config.maxNcv) primme.maxBasisSize = std::clamp<int>(static_cast<int>(config.maxNcv.value()), primme.numEvals + 1, primme.n);
    if(config.maxIter) primme.maxOuterIterations = config.maxIter.value();
    if(config.primme_projection) primme.projectionParams.projection = stringToProj(config.primme_projection);
    if(config.primme_locking) primme.locking = config.primme_locking.value();
    // Shifts
    switch(primme.target) {
        case primme_target::primme_largest:
        case primme_target::primme_smallest: {
            primme.numTargetShifts = 0;
            break;
        }
        case primme_target::primme_largest_abs: {
            // According to the manual, this is required when the target is primme_largest_abs
            config.primme_target_shifts = {0.0};
            primme.numTargetShifts      = 1;
            primme.targetShifts         = config.primme_target_shifts.data();
            break;
        }
        case primme_target::primme_closest_geq:
        case primme_target::primme_closest_leq:
        case primme_target::primme_closest_abs: {
            if(not config.primme_target_shifts.empty()) {
                eig::log->debug("Setting target shifts: {:.8f}", fmt::join(config.primme_target_shifts, ", "));
                primme.numTargetShifts = static_cast<int>(config.primme_target_shifts.size());
                primme.targetShifts    = config.primme_target_shifts.data();
            } else if(config.sigma and not matrix.isReadyShift()) {
                // We handle shifts by applying them directly on the matrix is possible. Else here:
                eig::log->debug("Setting target shift: {:.8f}", std::real(config.sigma.value()));
                config.primme_target_shifts = {std::real(config.sigma.value())};
                primme.numTargetShifts      = static_cast<int>(config.primme_target_shifts.size());
                primme.targetShifts         = config.primme_target_shifts.data();
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

    // Override some parameters
    if(config.primme_max_inner_iterations) primme.correctionParams.maxInnerIterations = config.primme_max_inner_iterations.value();
    //    primme.restartingParams.maxPrevRetain = 2;
    //    primme.minRestartSize = 2;
    //    primme.maxBlockSize = 1;

    // Allocate space
    auto &eigvals = result.get_eigvals<eig::Form::SYMM>();
    auto &eigvecs = result.get_eigvecs<Scalar, eig::Form::SYMM>();
    eigvals.resize(static_cast<size_t>(primme.numEvals));
    eigvecs.resize(static_cast<size_t>(primme.numEvals) * static_cast<size_t>(primme.n));
    result.meta.residual_norms.resize(static_cast<size_t>(primme.numEvals));

    // Copy initial guess
    primme.initSize = 0;
    for(const auto &ig : config.initial_guess) {
        if(ig.idx < primme.numEvals and ig.idx >= 0) {
            auto source = static_cast<Scalar *>(ig.ptr);
            auto target = eigvecs.data() + static_cast<size_t>(ig.idx * primme.n);
            std::copy_n(source, primme.n, target);
            primme.initSize = std::max(primme.initSize, static_cast<int>(ig.idx) + 1);
        }
    }

    /* Call primme  */
    int info = 0;
    if constexpr(std::is_same_v<Scalar, eig::real>) {
        try {
            info = dprimme(eigvals.data(), eigvecs.data(), result.meta.residual_norms.data(), &primme);
        } catch(const std::exception &ex) { eig::log->error("dprimme has exited with error: info {} | msg: {}", info, ex.what()); } catch(...) {
            eig::log->error("dprimme has exited with error: info {}", info);
        }
    } else if constexpr(std::is_same_v<typename MatrixProductType::Scalar, eig::cplx>) {
        try {
            info = zprimme(eigvals.data(), eigvecs.data(), result.meta.residual_norms.data(), &primme);
        } catch(const std::exception &ex) { eig::log->error("zprimme has exited with error: info {} | msg: {}", info, ex.what()); } catch(...) {
            eig::log->error("zprimme has exited with error: info {}", info);
        }
    } else {
        throw std::runtime_error("Primme: type not implemented {}", eig::sfinae::type_name<typename MatrixProductType::Scalar>());
    }

    t_tot_token.toc();
    if(info == 0) {
        switch(primme.dynamicMethodSwitch) {
            case -1: eig::log->debug("Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
            case -2: eig::log->debug("Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
            case -3: eig::log->debug("Recommended method for next run: DYNAMIC (close call)\n"); break;
        }
    } else {
        switch(info) {
            case -1: eig::log->error("PRIMME_UNEXPECTED_FAILURE (exit {}): set printLevel > 0 to see the call stack", info); break;
            case -2: eig::log->error("PRIMME_MALLOC_FAILURE (exit {}): eith<er CPU or GPU", info); break;
            case -3: eig::log->info("PRIMME_MAIN_ITER_FAILURE (exit {}): {}", info, getLogMessage(&primme)); break;
            case -4: eig::log->error("PRIMME_ARGUMENT_IS_NULL (exit {})", info); break;
            case -5: eig::log->error("PRIMME_INVALID_ARG n < 0 or nLocal < 0 or nLocal > n (exit {})", info); break;
            case -6: eig::log->error("PRIMME_INVALID_ARG numProcs < 1 (exit {})", info); break;
            case -40: eig::log->error("PRIMME_LAPACK_FAILURE (exit {})", info); break;
            case -41: eig::log->error("PRIMME_USER_FAILURE (exit {})", info); break;
            case -42: eig::log->error("PRIMME_ORTHO_CONST_FAILURE (exit {})", info); break;
            case -43: eig::log->error("PRIMME_PARALLEL_FAILURE (exit {})", info); break;
            case -44: eig::log->error("PRIMME_FUNCTION_UNAVAILABLE (exit {})", info); break;
            default: eig::log->error("Unknown error code: Go to http://www.cs.wm.edu/~andreas/software/doc/appendix.html?highlight=primme_main_iter_failure");
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
    result.meta.n              = primme.n;
    result.meta.tag            = config.tag;
    result.meta.ritz           = TargetToString(primme.target);
    result.meta.form           = config.form.value();
    result.meta.type           = config.type.value();
    result.meta.time_total     = t_tot.get_time();
    result.meta.last_res_norm  = *std::max_element(result.meta.residual_norms.begin(), result.meta.residual_norms.end());

    tid::get(config.tag) += t_tot.get_time();
    tid::get(config.tag).add_count(1ul);
    tid::get(fmt::format("{}.prep", config.tag)) += t_pre.get_time();
    tid::get(fmt::format("{}.prep", config.tag)).add_count(1ul);
    tid::get(fmt::format("{}.matvec", config.tag)) += primme.stats.timeMatvec;
    tid::get(fmt::format("{}.matvec", config.tag)).add_count(static_cast<size_t>(primme.stats.numMatvecs));
    if(config.primme_preconditioner) {
        tid::get(fmt::format("{}.precond", config.tag)) += primme.stats.timePrecond;
        tid::get(fmt::format("{}.precond", config.tag)).add_count(static_cast<size_t>(primme.stats.numPreconds));
    }
    primme_free(&primme);
    return info;
}

template int eig::solver::eigs_primme(MatVecMPO<eig::real> &matrix);
template int eig::solver::eigs_primme(MatVecMPO<eig::cplx> &matrix);
template int eig::solver::eigs_primme(MatVecDense<eig::real> &matrix);
template int eig::solver::eigs_primme(MatVecDense<eig::cplx> &matrix);
template int eig::solver::eigs_primme(MatVecSparse<eig::real> &matrix);
template int eig::solver::eigs_primme(MatVecSparse<eig::cplx> &matrix);
