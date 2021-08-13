#include "../log.h"
#include "../matvec/matvec_dense.h"
#include "../matvec/matvec_mps.h"
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
template void eig::solver::MultAx_wrapper<MatVecDense<eig::real>>(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr);
template void eig::solver::MultAx_wrapper<MatVecDense<eig::cplx>>(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr);
template void eig::solver::MultAx_wrapper<MatVecSparse<eig::real>>(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr);
template void eig::solver::MultAx_wrapper<MatVecSparse<eig::cplx>>(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr);
template void eig::solver::MultAx_wrapper<MatVecMps<eig::real>>(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr);
template void eig::solver::MultAx_wrapper<MatVecMps<eig::cplx>>(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr);

template<typename MatrixProductType>
void eig::solver::GradientConvTest(double *eval, void *evec, double *rNorm, int *isconv, struct primme_params *primme, int *ierr) {
    if(rNorm == nullptr) return;
    if(primme == nullptr) return;

    double problemNorm;
    if(not primme->massMatrixMatvec) {
        problemNorm = primme->aNorm > 0.0 ? primme->aNorm : primme->stats.estimateLargestSVal;
    } else {
        problemNorm = primme->aNorm > 0.0 && primme->invBNorm > 0.0 ? primme->aNorm * primme->invBNorm : primme->stats.estimateLargestSVal;
    }
    double prec         = primme->eps * problemNorm;
    int    default_crit = *rNorm < prec;

    *isconv = default_crit;
    *ierr   = 0;

    if(eval == nullptr) return;
    if(evec == nullptr) return;
    if(primme->convtest == nullptr) return;
    auto &solver = *static_cast<eig::solver *>(primme->convtest);
    auto &config = solver.config;
    auto &result = solver.result;

    // Store rNorm
    result.meta.last_res_norm = *rNorm;

    // Terminate if its taking too long
    if(config.maxTime.has_value() and primme->stats.elapsedTime > config.maxTime.value()) {
        eig::log->warn("primme: max time has been exeeded: {:.2f}", config.maxTime.value());
        primme->maxMatvecs = 0;
    }

    if constexpr(MatrixProductType::storage == eig::Storage::MPS) {
        if(config.primme_extra == nullptr) return;
        if(primme->matrix == nullptr) return;
        long   iter_since_last = primme->stats.numMatvecs - result.meta.last_iter_grad;
        double time_since_last = primme->stats.elapsedTime - result.meta.last_time_grad;

        double grad_tol  = config.primme_grad_tol ? config.primme_grad_tol.value() : 1e0;
        long   grad_iter = config.primme_grad_iter ? config.primme_grad_iter.value() : 100;
        double grad_time = config.primme_grad_time ? config.primme_grad_time.value() : 5.0;

        // If we do inner iterations we better check gradient convergence every outer iteration instead, so we override the
        // iter/time periods given in config.
        if(primme->correctionParams.maxInnerIterations != 0) {
            grad_iter = 0;
            grad_time = 0;
        }

        if(iter_since_last >= grad_iter or time_since_last >= grad_time) {
            auto H_ptr   = static_cast<MatrixProductType *>(config.primme_extra);
            auto H2_ptr  = static_cast<MatrixProductType *>(primme->matrix);
            using Scalar = typename MatrixProductType::Scalar;
            using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
            if(H2_ptr->rows() != H_ptr->rows()) throw std::logic_error("H2 and H size mismatch");
            H_ptr->compress();
            long mps_size = static_cast<long>(H2_ptr->rows());
            auto v        = Eigen::Map<Vector>(static_cast<Scalar *>(evec), mps_size);

            Vector Hv(mps_size), H2v(mps_size);

            H_ptr->MultAx(v.data(), Hv.data());
            H2_ptr->MultAx(v.data(), H2v.data());

            if(H_ptr->template get_shift<Scalar>() != 0.0) Hv += v * H_ptr->template get_shift<Scalar>();
            if(H2_ptr->template get_shift<Scalar>() != 0.0) H2v += v * H2_ptr->template get_shift<Scalar>();
            auto vHv      = v.dot(Hv);
            auto vH2v     = v.dot(H2v);
            auto var      = std::abs(vH2v - vHv * vHv);
            auto var_1    = 1.0 / var / std::log(10);
            auto norm_1   = 1.0 / v.norm();
            auto pref     = std::is_same_v<Scalar, cplx> ? 2.0 : 1.0; // Factor 2 for complex
            auto grad     = pref * var_1 * norm_1 * (H2v - 2.0 * vHv * Hv - (vH2v - 2.0 * vHv * vHv) * v);
            auto grad_max = grad.template lpNorm<Eigen::Infinity>();

            // Store gradient info
            result.meta.last_grad_max  = grad_max;
            result.meta.last_iter_grad = primme->stats.numMatvecs;
            result.meta.last_time_grad = primme->stats.elapsedTime;

            *isconv = std::max<int>(*isconv, grad_max < grad_tol);
            *ierr   = 0;
            //            if(*isconv != 0)
            eig::log->info(FMT_STRING("ops {:<5} | iter {:<4} | f {:20.16f} | ∇fₘₐₓ {:8.2e} | res {:8.2e} | time {:8.2f} s | dt {:8.2f} ms/op"),
                           primme->stats.numMatvecs, primme->stats.numOuterIterations, primme->stats.estimateMinEVal, grad_max, *rNorm,
                           primme->stats.elapsedTime, primme->stats.timeMatvec / primme->stats.numMatvecs * 1000);
        }
    }
}

void monitorFun([[maybe_unused]] void *basisEvals, [[maybe_unused]] int *basisSize, [[maybe_unused]] int *basisFlags, [[maybe_unused]] int *iblock,
                [[maybe_unused]] int *blockSize, [[maybe_unused]] void *basisNorms, [[maybe_unused]] int *numConverged, [[maybe_unused]] void *lockedEvals,
                [[maybe_unused]] int *numLocked, [[maybe_unused]] int *lockedFlags, [[maybe_unused]] void *lockedNorms, [[maybe_unused]] int *inner_its,
                [[maybe_unused]] void *LSRes, [[maybe_unused]] const char *msg, [[maybe_unused]] double *time, [[maybe_unused]] primme_event *event,
                [[maybe_unused]] struct primme_params *primme, [[maybe_unused]] int *ierr) {
    if(event == nullptr) return;
    if(*event == primme_event_inner_iteration) return; // No need to log that often

    std::string event_msg;
    /* clang-format off */
    if(*event == primme_event_outer_iteration)  event_msg = "event_outer_iteration";
    if(*event == primme_event_restart)          event_msg = "event_restart";
    if(*event == primme_event_reset)            event_msg = "event_reset";
    if(*event == primme_event_converged)        event_msg = "event_converged";
    if(*event == primme_event_locked)           event_msg = "event_locked";
    if(*event == primme_event_message)          event_msg = "event_message";
    if(*event == primme_event_profile)          event_msg = "event_profile";
    /* clang-format on */

    std::string basisMsg, blockMsg;
    if(basisSize != nullptr) basisMsg = fmt::format("| basis {:<3}", *basisSize);
    if(blockSize != nullptr) blockMsg = fmt::format("| block {:<3}", *blockSize);
    if(msg != nullptr) { event_msg += ": " + std::string(msg); }

    if(primme->monitor == nullptr) {
        eig::log->trace(FMT_STRING("ops {:<5} | iter {:<4} | f {:20.16f} | time {:8.2f} s | dt {:8.2f} ms/op {} {} | {}"), primme->stats.numMatvecs,
                        primme->stats.numOuterIterations, primme->stats.estimateMinEVal, primme->stats.elapsedTime,
                        primme->stats.timeMatvec / primme->stats.numMatvecs * 1000, basisMsg, blockMsg, event_msg);
        return;
    }

    auto &solver = *static_cast<eig::solver *>(primme->monitor);
    auto &config = solver.config;
    auto &result = solver.result;

    double primme_log_time = config.logTime ? config.logTime.value() : std::numeric_limits<double>::quiet_NaN();
    double time_since_last_log = std::abs(primme->stats.elapsedTime - result.meta.last_time_log);

    if(time_since_last_log > primme_log_time or *event == primme_event_converged) {
        auto level = (config.logTime and time_since_last_log > 60) ? spdlog::level::info : spdlog::level::debug;
        eig::log->log(level, FMT_STRING("ops {:<5} | iter {:<4} | f {:20.16f} | ∇fₘₐₓ {:8.2e} | res {:8.2e} | time {:8.2f} s | dt {:8.2f} ms/op {} {} | {}"),
                      primme->stats.numMatvecs, primme->stats.numOuterIterations, primme->stats.estimateMinEVal, result.meta.last_grad_max,
                      result.meta.last_res_norm, primme->stats.elapsedTime, primme->stats.timeMatvec / primme->stats.numMatvecs * 1000, basisMsg, blockMsg,
                      event_msg);
        result.meta.last_time_log = primme->stats.elapsedTime;
        result.meta.last_iter_log = primme->stats.numMatvecs;
    }
    // Terminate if its taking too long
    if(config.maxTime.has_value() and primme->stats.elapsedTime > config.maxTime.value()) {
        eig::log->warn("primme: max time has been exeeded: {:.2f}", config.maxTime.value());
        primme->maxMatvecs = 0;
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

    // Set the function which prints log output
    primme.monitorFun = monitorFun;

    // Make eig objects visible from within monitorFun
    primme.monitor = this;

    // Set the function to make convergence tests
    if constexpr(MatrixProductType::storage == eig::Storage::MPS) {
        primme.convTestFun = eig::solver::GradientConvTest<MatrixProductType>;
        primme.convtest    = this;
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
    std::vector<double> tgtShifts = {0.0};
    switch(primme.target) {
        case primme_target::primme_largest:
        case primme_target::primme_smallest: {
            primme.numTargetShifts = 0;
            break;
        }
        case primme_target::primme_largest_abs: {
            // According to the manual, this is required when the target is primme_largest_abs
            primme.numTargetShifts = 1;
            primme.targetShifts    = tgtShifts.data();
            break;
        }
        case primme_target::primme_closest_geq:
        case primme_target::primme_closest_leq:
        case primme_target::primme_closest_abs: {
            primme.numTargetShifts = 1;
            primme.targetShifts    = tgtShifts.data();
            if(config.sigma and not matrix.isReadyShift()) {
                // We handle shifts by applying them directly on the matrix is possible. Else here:
                primme.targetShifts[0] = std::real(config.sigma.value());
            }
            break;
        }
    }

    /* Set method to solve the problem
     * DYNAMIC uses a runtime heuristic to choose the fastest method between
     * PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
     * set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */
    primme_set_method(stringToMethod(config.primme_method), &primme);
    // Override some parameters
    //    primme.restartingParams.maxPrevRetain = 2;
    //    primme.minRestartSize = 2;
    //    primme.maxBlockSize = 32;
    if(config.primme_max_inner_iterations) primme.correctionParams.maxInnerIterations = config.primme_max_inner_iterations.value();

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
            auto target = eigvecs.data() + ig.idx * static_cast<size_t>(primme.n);
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
        std::string msg;
        switch(info) {
            case -1: msg = "PRIMME_UNEXPECTED_FAILURE: set printLevel > 0 to see the call stack"; break;
            case -2: msg = "PRIMME_MALLOC_FAILURE: either CPU or GPU"; break;
            case -3: msg = "PRIMME_MAIN_ITER_FAILURE: maxOuterIterations or maxMatvecs reached"; break;
            case -4: msg = "Argument primme is null"; break;
            case -5: msg = "n < 0 or nLocal < 0 or nLocal > n"; break;
            case -6: msg = "numProcs < 1"; break;
            case -40: msg = "PRIMME_LAPACK_FAILURE"; break;
            case -41: msg = "PRIMME_USER_FAILURE"; break;
            case -42: msg = "PRIMME_ORTHO_CONST_FAILURE"; break;
            case -43: msg = "PRIMME_PARALLEL_FAILURE"; break;
            case -44: msg = "PRIMME_FUNCTION_UNAVAILABLE"; break;
            default: msg = "Unknown error code: Go to http://www.cs.wm.edu/~andreas/software/doc/appendix.html?highlight=primme_main_iter_failure";
        }
        eig::log->warn("Primme returned with nonzero exit status {}: {}", info, msg);
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
    result.meta.matvecs        = primme.stats.numMatvecs;
    result.meta.n              = primme.n;
    result.meta.tag            = config.tag;
    result.meta.ritz           = TargetToString(primme.target);
    result.meta.form           = config.form.value();
    result.meta.type           = config.type.value();
    result.meta.time_total     = t_tot.get_time();

    tid::get("primme") += t_tot.get_time();
    tid::get("primme.prep") += t_pre.get_time();
    tid::get("primme.elapsed") += primme.stats.elapsedTime;
    tid::get("primme.matvec") += primme.stats.timeMatvec;

    primme_free(&primme);
    return info;
}

template int eig::solver::eigs_primme(MatVecMps<eig::real> &matrix);
template int eig::solver::eigs_primme(MatVecMps<eig::cplx> &matrix);
template int eig::solver::eigs_primme(MatVecDense<eig::real> &matrix);
template int eig::solver::eigs_primme(MatVecDense<eig::cplx> &matrix);
template int eig::solver::eigs_primme(MatVecSparse<eig::real> &matrix);
template int eig::solver::eigs_primme(MatVecSparse<eig::cplx> &matrix);
