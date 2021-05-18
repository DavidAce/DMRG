#include "math/eig/matvec/matvec_dense.h"
#include "math/eig/matvec/matvec_mpo.h"
#include "math/eig/matvec/matvec_sparse.h"
#include "sfinae.h"
#include "solver.h"
#include <general/class_tic_toc.h>
#include <primme/primme.h>
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

static MatVecMPO<eig::real> *   matvecmpo_real_ptr    = nullptr;
static MatVecMPO<eig::cplx> *   matvecmpo_cplx_ptr    = nullptr;
static MatVecDense<eig::real> * matvecdense_real_ptr  = nullptr;
static MatVecDense<eig::cplx> * matvecdense_cplx_ptr  = nullptr;
static MatVecSparse<eig::real> *matvecsparse_real_ptr = nullptr;
static MatVecSparse<eig::cplx> *matvecsparse_cplx_ptr = nullptr;

void MultAx_wrapper(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
    if(matvecmpo_real_ptr != nullptr) return matvecmpo_real_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecmpo_cplx_ptr != nullptr) return matvecmpo_cplx_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecdense_real_ptr != nullptr) return matvecdense_real_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecdense_cplx_ptr != nullptr) return matvecdense_cplx_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecsparse_real_ptr != nullptr) return matvecsparse_real_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    if(matvecsparse_cplx_ptr != nullptr) return matvecsparse_cplx_ptr->MultAx(x, ldx, y, ldy, blockSize, primme, ierr);
    throw std::runtime_error("No pointer was set");
}

template<typename MatrixProductType>
int eig::solver::eigs_primme(MatrixProductType &matrix) {
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
    class_tic_toc t_tot(true, 5, "");
    auto          t_tot_token = t_tot.tic_token();

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

    /* Set problem parameters */
    primme.n            = matrix.rows();                           /* set problem dimension */
    primme.numEvals     = static_cast<int>(config.maxNev.value()); /* Number of wanted eigenpairs */
    primme.eps          = config.tol.value();                      /* ||r|| <= eps * ||matrix|| */
    primme.target       = RitzToTarget(config.ritz.value());
    primme.maxBasisSize = std::clamp<int>(static_cast<int>(config.maxNcv.value()), static_cast<int>(config.maxNev.value()) + 1, matrix.rows());
    if(primme.target == primme_largest_abs || primme.target == primme_closest_geq || primme.target == primme_closest_leq ||
       primme.target == primme_closest_abs) {
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
    primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
    int                    info = 0;
    std::vector<eig::real> rnorms(static_cast<size_t>(primme.numEvals));
    //    /* Allocate space for converged Ritz values and residual norms */
    if constexpr(std::is_same_v<typename MatrixProductType::Scalar, eig::real>) {
        auto &eigvals = result.get_eigvals<eig::Form::SYMM>();
        auto &eigvecs = result.get_eigvecs<eig::Form::SYMM, eig::Type::REAL>();
        eigvals.resize(static_cast<size_t>(primme.numEvals));
        eigvecs.resize(static_cast<size_t>(primme.numEvals) * static_cast<size_t>(primme.n));

        /* Call primme  */
        info = dprimme(eigvals.data(), eigvecs.data(), rnorms.data(), &primme);
    } else if constexpr(std::is_same_v<typename MatrixProductType::Scalar, eig::cplx>) {
        auto &eigvals = result.get_eigvals<eig::Form::SYMM>();
        auto &eigvecs = result.get_eigvecs<eig::Form::SYMM, eig::Type::CPLX>();
        eigvals.resize(static_cast<size_t>(primme.numEvals));
        eigvecs.resize(static_cast<size_t>(primme.numEvals) * static_cast<size_t>(primme.n));

        /* Call primme  */
        info = zprimme(eigvals.data(), eigvecs.data(), rnorms.data(), &primme);
    } else {
        throw std::runtime_error("Primme: type not implemented {}", eig::sfinae::type_name<typename MatrixProductType::Scalar>());
    }

    t_tot_token.toc();
    if(info == 0) {
        result.meta.eigvecsR_found = true;
        result.meta.eigvals_found  = true;
        result.meta.rows           = primme.n;
        result.meta.cols           = primme.numEvals;
        result.meta.nev            = primme.numEvals;
        result.meta.nev_converged  = primme.numEvals;
        result.meta.ncv            = primme.maxBasisSize;
        result.meta.tol            = primme.eps;
        result.meta.iter           = primme.stats.numOuterIterations;
        result.meta.counter        = primme.stats.numMatvecs;
        result.meta.n              = primme.n;
        result.meta.ritz           = TargetToString(primme.target);
        result.meta.form           = config.form.value();
        result.meta.type           = config.type.value();
        result.meta.time_total     = t_tot.get_measured_time();

        switch(primme.dynamicMethodSwitch) {
            case -1: eig::log->debug("Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
            case -2: eig::log->debug("Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
            case -3: eig::log->debug("Recommended method for next run: DYNAMIC (close call)\n"); break;
        }

        primme_free(&primme);
    } else {
        primme_free(&primme);
        throw std::runtime_error("Primme returned with nonzero exit status: " + std::to_string(info));
    }

    matvecmpo_real_ptr    = nullptr;
    matvecmpo_cplx_ptr    = nullptr;
    matvecdense_real_ptr  = nullptr;
    matvecdense_cplx_ptr  = nullptr;
    matvecsparse_real_ptr = nullptr;
    matvecsparse_cplx_ptr = nullptr;

    return info;
}

template int eig::solver::eigs_primme(MatVecMPO<eig::real> &matrix);
template int eig::solver::eigs_primme(MatVecMPO<eig::cplx> &matrix);
template int eig::solver::eigs_primme(MatVecDense<eig::real> &matrix);
template int eig::solver::eigs_primme(MatVecDense<eig::cplx> &matrix);
template int eig::solver::eigs_primme(MatVecSparse<eig::real> &matrix);
template int eig::solver::eigs_primme(MatVecSparse<eig::cplx> &matrix);
