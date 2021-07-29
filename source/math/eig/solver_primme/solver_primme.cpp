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
    if(ritzstring == "SM") return primme_target::primme_closest_abs;
    if(ritzstring == "LR") return primme_target::primme_largest;
    if(ritzstring == "SR") return primme_target::primme_smallest;
    if(ritzstring == "LI") throw std::runtime_error("Primme RITZ LI (largest imaginary) not implemented");
    if(ritzstring == "SI") throw std::runtime_error("Primme RITZ SI (smallest imaginary) not implemented");
    if(ritzstring == "BE") throw std::runtime_error("Primme RITZ BE (both ends) not implemented");
    throw std::runtime_error("Wrong ritz string: " + std::string(ritzstring));
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
    throw std::runtime_error("Wrong method string: " + std::string(methodstring.value()));
}

inline std::string_view TargetToString(primme_target target) {
    if(target == primme_target::primme_largest) return "LA";
    if(target == primme_target::primme_smallest) return "SA";
    if(target == primme_target::primme_largest_abs) return "LM";
    if(target == primme_target::primme_closest_abs) return "SM";
    if(target == primme_target::primme_largest) return "LR";
    if(target == primme_target::primme_smallest) return "SR";
    throw std::logic_error("Could not match primme_target");
}

inline primme_target RitzToTarget(eig::Ritz ritz) {
    if(ritz == eig::Ritz::LA) return primme_target::primme_largest;
    if(ritz == eig::Ritz::SA) return primme_target::primme_smallest;
    if(ritz == eig::Ritz::LM) return primme_target::primme_largest_abs;
    if(ritz == eig::Ritz::SM) return primme_target::primme_closest_abs;
    if(ritz == eig::Ritz::LR) return primme_target::primme_largest;
    if(ritz == eig::Ritz::SR) return primme_target::primme_smallest;
    if(ritz == eig::Ritz::LI) throw std::runtime_error("Primme RITZ LI (largest imaginary) not implemented");
    if(ritz == eig::Ritz::SI) throw std::runtime_error("Primme RITZ SI (smallest imaginary) not implemented");
    if(ritz == eig::Ritz::BE) throw std::runtime_error("Primme RITZ BE (both ends) not implemented");
    throw std::runtime_error("Wrong ritz enum");
}

static MatVecMPO<eig::real>    *matvecmpo_real_ptr    = nullptr;
static MatVecMPO<eig::cplx>    *matvecmpo_cplx_ptr    = nullptr;
static MatVecDense<eig::real>  *matvecdense_real_ptr  = nullptr;
static MatVecDense<eig::cplx>  *matvecdense_cplx_ptr  = nullptr;
static MatVecSparse<eig::real> *matvecsparse_real_ptr = nullptr;
static MatVecSparse<eig::cplx> *matvecsparse_cplx_ptr = nullptr;

static double primme_max_time = 60 * 60;

void MultAx_wrapper(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
    if(matvecmpo_real_ptr != nullptr) return matvecmpo_real_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecmpo_cplx_ptr != nullptr) return matvecmpo_cplx_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecdense_real_ptr != nullptr) return matvecdense_real_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecdense_cplx_ptr != nullptr) return matvecdense_cplx_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecsparse_real_ptr != nullptr) return matvecsparse_real_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecsparse_cplx_ptr != nullptr) return matvecsparse_cplx_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    throw std::runtime_error("No pointer was set");
}

void monitorFun([[maybe_unused]] void *basisEvals, [[maybe_unused]] int *basisSize, [[maybe_unused]] int *basisFlags, [[maybe_unused]] int *iblock,
                [[maybe_unused]] int *blockSize, [[maybe_unused]] void *basisNorms, [[maybe_unused]] int *numConverged, [[maybe_unused]] void *lockedEvals,
                [[maybe_unused]] int *numLocked, [[maybe_unused]] int *lockedFlags, [[maybe_unused]] void *lockedNorms, [[maybe_unused]] int *inner_its,
                [[maybe_unused]] void *LSRes, [[maybe_unused]] const char *msg, [[maybe_unused]] double *time, [[maybe_unused]] primme_event *event,
                [[maybe_unused]] class primme_params *primme, [[maybe_unused]] int *ierr) {
    std::string event_msg;
    switch(*event) {
        case primme_event_outer_iteration: {
            event_msg = "primme_outer_iteration";
            break;
        } /* report at every outer iteration        */
        case primme_event_inner_iteration: {
            event_msg = "primme_inner_iteration";
            break;
        } /* report at every QMR iteration          */
        case primme_event_restart: {
            event_msg = "primme_restart";
            break;
        } /* report at every basis restart          */
        case primme_event_reset: {
            event_msg = "primme_reset";
            break;
        } /* event launch if basis reset            */
        case primme_event_converged: {
            event_msg = "primme_converged";
            break;
        } /* report new pair marked as converged    */
        case primme_event_locked: {
            event_msg = "primme_locked";
            break;
        } /* report new pair marked as locked       */
        case primme_event_message: {
            event_msg = "primme_message";
            break;
        } /* report warning                         */
        case primme_event_profile: {
            event_msg = "primme_profile";
            break;
        } /* report time from consumed by a function*/
    }
    static double last_log_time = 0;
    if(last_log_time > primme->stats.elapsedTime) last_log_time = 0; // If this function has been run before

    bool log_time = std::abs(primme->stats.elapsedTime - last_log_time) > 5.0 or last_log_time == 0;
    bool log_init = primme->stats.numOuterIterations <= 1 or primme->stats.numMatvecs <= 1;

    if(log_time or log_init) {
        eig::log->info("iter {:<4} | nops {:<5} | min eval {:20.16f} | time {:8.2f} s | dt {:8.2f} ms/op | block {} | basis {} | {}",
                       primme->stats.numOuterIterations, primme->stats.numMatvecs, primme->stats.estimateMinEVal, primme->stats.elapsedTime,
                       primme->stats.timeMatvec / primme->stats.numMatvecs * 1000, *blockSize, *basisSize, event_msg);
        last_log_time = primme->stats.elapsedTime;
    }
    if(primme->stats.elapsedTime > primme_max_time) primme->maxMatvecs = 0;
}

template<typename MatrixProductType>
int eig::solver::eigs_primme(MatrixProductType &matrix, typename MatrixProductType::Scalar *residual) {
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
     * We need to bind the matrix-vector product operation using std::function.
     * This is a delicate step.
     */
    if constexpr(std::is_same_v<MatrixProductType, MatVecMPO<eig::real>>)
        matvecmpo_real_ptr = &matrix;
    else if constexpr(std::is_same_v<MatrixProductType, MatVecMPO<eig::cplx>>)
        matvecmpo_cplx_ptr = &matrix;
    else if constexpr(std::is_same_v<MatrixProductType, MatVecDense<eig::real>>)
        matvecdense_real_ptr = &matrix;
    else if constexpr(std::is_same_v<MatrixProductType, MatVecDense<eig::cplx>>)
        matvecdense_cplx_ptr = &matrix;
    else if constexpr(std::is_same_v<MatrixProductType, MatVecSparse<eig::real>>)
        matvecsparse_real_ptr = &matrix;
    else if constexpr(std::is_same_v<MatrixProductType, MatVecSparse<eig::cplx>>)
        matvecsparse_cplx_ptr = &matrix;
    else
        throw std::runtime_error("No matching MatVec operation for the wrapper was found");
    primme.matrixMatvec = MultAx_wrapper;
    primme.monitorFun   = monitorFun;

    /* Set problem parameters */
    primme.n        = matrix.rows();                           /* set problem dimension */
    primme.numEvals = static_cast<int>(config.maxNev.value()); /* Number of wanted eigenpairs */
    if(config.tol) primme.eps = config.tol.value();            /* ||r|| <= eps * ||matrix|| */
    primme.target       = RitzToTarget(config.ritz.value());
    primme.maxBasisSize = std::clamp<int>(static_cast<int>(config.maxNcv.value()), primme.numEvals + 1, primme.n);
    if(config.maxIter) primme.maxOuterIterations = config.maxIter.value();
    if(config.maxTime) primme_max_time = config.maxTime.value();

    // Shifts
    std::vector<double> tgtShifts = {0.0};
    if(primme.target == primme_largest_abs || primme.target == primme_closest_geq || primme.target == primme_closest_leq ||
       primme.target == primme_closest_abs) {
        if(config.sigma and not matrix.isReadyShift()) {
            // We handle shifts by applying them directly on the matrix is possible. Else here:
            primme.numTargetShifts = 1;
            tgtShifts              = {std::real(config.sigma.value())};
            primme.targetShifts    = tgtShifts.data();
        } else {
            // According to the manual, this is required when the target is primme_largest_abs
            primme.numTargetShifts = 1;
            primme.targetShifts    = tgtShifts.data();
        }

        if(primme.numTargetShifts <= 0)
            throw std::runtime_error(fmt::format("Primme ritz {} requires well defined numTargetShifts | got {}", config.ritz.value(), primme.numTargetShifts));
        else if(primme.targetShifts == nullptr)
            throw std::runtime_error(
                fmt::format("Primme ritz {} requires well defined targetShifts | got nullptr", config.ritz.value(), primme.numTargetShifts));
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

    int                    info = 0;
    std::vector<eig::real> rnorms(static_cast<size_t>(primme.numEvals));
    //    /* Allocate space for converged Ritz values and residual norms */
    if constexpr(std::is_same_v<typename MatrixProductType::Scalar, eig::real>) {
        auto &eigvals = result.get_eigvals<eig::Form::SYMM>();
        auto &eigvecs = result.get_eigvecs<eig::Form::SYMM, eig::Type::REAL>();
        eigvals.resize(static_cast<size_t>(primme.numEvals));
        eigvecs.resize(static_cast<size_t>(primme.numEvals) * static_cast<size_t>(primme.n));
        // Initial guess / residual
        if(residual != nullptr) {
            primme.initSize = 1;
            std::copy_n(residual, primme.n, eigvecs.data());
        }

        /* Call primme  */
        try {
            info = dprimme(eigvals.data(), eigvecs.data(), rnorms.data(), &primme);
        } catch(const std::exception &ex) { eig::log->error("dprimme has exited with error: info {} | msg: {}", info, ex.what()); } catch(...) {
            eig::log->error("dprimme has exited with error: info {}", info);
        }
    } else if constexpr(std::is_same_v<typename MatrixProductType::Scalar, eig::cplx>) {
        auto &eigvals = result.get_eigvals<eig::Form::SYMM>();
        auto &eigvecs = result.get_eigvecs<eig::Form::SYMM, eig::Type::CPLX>();
        eigvals.resize(static_cast<size_t>(primme.numEvals));
        eigvecs.resize(static_cast<size_t>(primme.numEvals) * static_cast<size_t>(primme.n));
        // Initial guess / residual
        if(residual != nullptr) {
            primme.initSize = 1;
            //            primme.numOrthoConst = 1;
            std::copy_n(residual, primme.n, eigvecs.data());
        }

        /* Call primme  */
        try {
            info = zprimme(eigvals.data(), eigvecs.data(), rnorms.data(), &primme);
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
    result.meta.counter        = primme.stats.numMatvecs;
    result.meta.n              = primme.n;
    result.meta.ritz           = TargetToString(primme.target);
    result.meta.form           = config.form.value();
    result.meta.type           = config.type.value();
    result.meta.time_total     = t_tot.get_time();

    tid::get("primme") += t_tot.get_time();
    tid::get("primme.prep") += t_pre.get_time();
    tid::get("primme.elapsed") += primme.stats.elapsedTime;
    tid::get("primme.matvec") += primme.stats.timeMatvec;

    primme_free(&primme);
    matvecmpo_real_ptr    = nullptr;
    matvecmpo_cplx_ptr    = nullptr;
    matvecdense_real_ptr  = nullptr;
    matvecdense_cplx_ptr  = nullptr;
    matvecsparse_real_ptr = nullptr;
    matvecsparse_cplx_ptr = nullptr;

    return info;
}

template int eig::solver::eigs_primme(MatVecMPO<eig::real> &matrix, eig::real *residual);
template int eig::solver::eigs_primme(MatVecMPO<eig::cplx> &matrix, eig::cplx *residual);
template int eig::solver::eigs_primme(MatVecDense<eig::real> &matrix, eig::real *residual);
template int eig::solver::eigs_primme(MatVecDense<eig::cplx> &matrix, eig::cplx *residual);
template int eig::solver::eigs_primme(MatVecSparse<eig::real> &matrix, eig::real *residual);
template int eig::solver::eigs_primme(MatVecSparse<eig::cplx> &matrix, eig::cplx *residual);
