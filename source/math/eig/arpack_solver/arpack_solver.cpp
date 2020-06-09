//
// Created by david on 2018-10-30.
//

#include "arpack_solver.h"
#include "matrix_product_dense.h"
#include "matrix_product_hamiltonian.h"
#include "matrix_product_sparse.h"
#include <general/nmspc_sfinae.h>

#if defined(_MKL_LAPACK_H_)
    #error MKL_LAPACK IS NOT SUPPOSED TO BE DEFINED HERE
#endif

#if defined(LAPACK_H)
    #error LAPACK IS NOT SUPPOSED TO BE DEFINED HERE
#endif

#if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Woverloaded-virtual"
#elif defined(__GNUC__) || defined(__GNUG__)
#endif

#if __has_include(<arpackpp/arcomp.h>)
    #include <arpackpp/arcomp.h>
    #include <arpackpp/ardnsmat.h>
    #include <arpackpp/ardscomp.h>
    #include <arpackpp/arsnsym.h>
    #include <arpackpp/arssym.h>
#elif __has_include(<arpack++/arcomp.h>)
    #include <arpack++/ardnsmat.h>
    #include <arpack++/ardscomp.h>
    #include <arpack++/arsnsym.h>
    #include <arpack++/arssym.h>
#else
    #error Could not include arpack headers correctly
#endif

#if defined(__clang__)
    // turn the warnings back on
    #pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
#endif

namespace tc = sfinae;
using namespace eig;

template<typename MatrixType>
eig::arpack_solver<MatrixType>::arpack_solver(MatrixType &matrix_, eig::settings &config_, eig::solution &result_, Scalar *residual_)
    : matrix(matrix_), config(config_), result(result_), residual(residual_) {
    t_sol.set_properties(profile_arpack, 10, "Time iterating  ");
    t_get.set_properties(profile_arpack, 10, "Time getting result");
    t_sub.set_properties(profile_arpack, 10, "Time subtracting");
    t_all.set_properties(profile_arpack, 10, "Time doing all  ");
    nev_internal = std::min(matrix.rows() / 2, static_cast<int>(config.eigMaxNev.value()));
    ncv_internal = static_cast<int>(std::max(config.eigMaxNcv.value(), 2 + config.eigMaxNev.value()));
    ncv_internal = std::min(ncv_internal, matrix.rows());
}

template<typename MatrixType>
void eig::arpack_solver<MatrixType>::eigs() {
    auto loglevel = Logger::getLogLevel(eig::log);
    eig::log      = Logger::setLogger("eigs", loglevel);

    result.reset();
    nev_internal = std::min(matrix.rows() / 2, static_cast<int>(config.eigMaxNev.value()));
    ncv_internal = static_cast<int>(std::max(config.eigMaxNcv.value(), 2 + config.eigMaxNev.value()));
    ncv_internal = std::min(ncv_internal, matrix.rows());
    assert(nev_internal >= 1 and nev_internal <= matrix.rows() / 2);

    config.checkRitz();
    matrix.set_mode(config.form.value());
    matrix.set_side(config.side.value());

    // Dispatch to symmetric or nonsymmetric. If complex, there's only a nonsymmetric option available.
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        eigs_comp();
    } else {
        if(config.form == Form::SYMM) eigs_sym();
        else
            eigs_nsym();
    }
    eig::log = Logger::setLogger("eigs", loglevel);
}

template<typename MatrixType>
void eig::arpack_solver<MatrixType>::eigs_sym() {
    if constexpr(std::is_same<Scalar, double>::value) {
        assert(config.form == Form::SYMM and "ERROR: config not SYMMETRIC");
        assert(matrix.get_form() == Form::SYMM and "ERROR: matrix not SYMMETRIC");
        ARSymStdEig<double, MatrixType> solver(matrix.rows(), nev_internal, &matrix, &MatrixType::MultAx, config.get_ritz_string().c_str(), ncv_internal,
                                               config.eigThreshold.value(), static_cast<int>(config.eigMaxIter.value()), residual);

        if(config.sigma) {
            if(config.shift_invert == Shinv::ON) {
                if constexpr(MatrixType::can_shift_invert) {
                    eig::log->trace("Enabling shift-invert mode with sigma = {}", std::real(config.sigma.value()));
                    // Calculate shift-inverse mat-vec mult operator by LU decomposition
                    matrix.set_shift(config.sigma.value());
                    matrix.FactorOP();
                    if constexpr(std::is_same_v<Scalar, cplx>) solver.SetShiftInvertMode(config.sigma.value(), &matrix, &MatrixType::MultOPv);
                    if constexpr(std::is_same_v<Scalar, real>) solver.SetShiftInvertMode(std::real(config.sigma.value()), &matrix, &MatrixType::MultOPv);
                } else
                    throw std::runtime_error("Tried to shift-invert an incompatible matrix");
            } else {
                if constexpr(MatrixType::can_shift) {
                    eig::log->trace("Setting shift with sigma = {}", std::real(config.sigma.value()));
                    matrix.set_shift(config.sigma.value());
                } else
                    throw std::runtime_error("Tried to apply shift on an incompatible matrix");
            }
        }
        find_solution(solver, config.eigMaxNev.value());
        copy_solution(solver);

    } else {
        eig::log->critical("Called eigs_sym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
        throw std::runtime_error("Called eigs_sym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}

template<typename MatrixType>
void eig::arpack_solver<MatrixType>::eigs_nsym() {
    if constexpr(std::is_same<Scalar, double>::value) {
        assert(config.form == Form::NSYM and "ERROR: config not NSYM");
        assert(matrix.get_form() == Form::NSYM and "ERROR: matrix not NSYM");
        if(nev_internal == 1) { nev_internal++; }

        ARNonSymStdEig<double, MatrixType> solver(matrix.rows(), nev_internal, &matrix, &MatrixType::MultAx, config.get_ritz_string().c_str(), ncv_internal,
                                                  config.eigThreshold.value(), config.eigMaxIter.value(), residual);

        if(config.sigma) {
            if(config.shift_invert == Shinv::ON) {
                if constexpr(MatrixType::can_shift_invert) {
                    eig::log->trace("Enabling shift-invert mode with sigma = {}", std::real(config.sigma.value()));
                    // Calculate shift-inverse mat-vec mult operator by LU decomposition
                    matrix.set_shift(config.sigma.value());
                    matrix.FactorOP();
                    if constexpr(std::is_same_v<Scalar, cplx>) solver.SetShiftInvertMode(config.sigma.value(), &matrix, &MatrixType::MultOPv);
                    if constexpr(std::is_same_v<Scalar, real>) solver.SetShiftInvertMode(std::real(config.sigma.value()), &matrix, &MatrixType::MultOPv);
                } else
                    throw std::runtime_error("Tried to shift-invert an incompatible matrix");
            } else {
                if constexpr(MatrixType::can_shift) {
                    eig::log->trace("Setting shift with sigma = {}", std::real(config.sigma.value()));
                    matrix.set_shift(config.sigma.value());
                } else
                    throw std::runtime_error("Tried to apply shift on an incompatible matrix");
            }
        }

        find_solution(solver, config.eigMaxNev.value());
        copy_solution(solver);
    } else {
        throw std::runtime_error("Called eigs_nsym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}

template<typename MatrixType>
void eig::arpack_solver<MatrixType>::eigs_comp() {
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        ARCompStdEig<double, MatrixType> solver(matrix.rows(), nev_internal, &matrix, &MatrixType::MultAx, config.get_ritz_string().c_str(), ncv_internal,
                                                config.eigThreshold.value(), config.eigMaxIter.value(), residual);
        if(config.sigma) {
            if(config.shift_invert == Shinv::ON) {
                if constexpr(MatrixType::can_shift_invert) {
                    eig::log->trace("Enabling shift-invert mode with sigma = {}", std::real(config.sigma.value()));
                    // Calculate shift-inverse mat-vec mult operator by LU decomposition
                    matrix.set_shift(config.sigma.value());
                    matrix.FactorOP();
                    if constexpr(std::is_same_v<Scalar, cplx>) solver.SetShiftInvertMode(config.sigma.value(), &matrix, &MatrixType::MultOPv);
                    if constexpr(std::is_same_v<Scalar, real>) solver.SetShiftInvertMode(std::real(config.sigma.value()), &matrix, &MatrixType::MultOPv);
                } else
                    throw std::runtime_error("Tried to shift-invert an incompatible matrix");
            } else {
                if constexpr(MatrixType::can_shift) {
                    eig::log->trace("Setting shift with sigma = {}", std::real(config.sigma.value()));
                    matrix.set_shift(config.sigma.value());
                } else
                    throw std::runtime_error("Tried to apply shift on an incompatible matrix");
            }
        }
        find_solution(solver, config.eigMaxNev.value());
        copy_solution(solver);
    } else {
        throw std::runtime_error("Called eigs_nsym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}

//
//
template<typename MatrixType>
template<typename Derived>
void eig::arpack_solver<MatrixType>::find_solution(Derived &solver, eig::size_type nev) {
    if(config.compute_eigvecs) {
        eig::log->trace("Finding eigenvectors");
        solver.FindEigenvectors();
        if(config.side == eig::Side::L) result.meta.eigvecsL_found = solver.EigenvectorsFound();
        else
            result.meta.eigvecsR_found = solver.EigenvectorsFound(); // BOOL!

        result.meta.eigvals_found = solver.EigenvaluesFound(); // BOOL!
        result.meta.iter          = solver.GetIter();
        result.meta.n             = solver.GetN();
        result.meta.nev           = std::min(nev, static_cast<eig::size_type>(solver.GetNev()));
        result.meta.nev_converged = solver.ConvergedEigenvalues();
        result.meta.ncv_used      = solver.GetNcv();
        result.meta.rows          = solver.GetN();
        result.meta.cols          = result.meta.nev;
        result.meta.counter       = matrix.counter;
        result.meta.form          = config.form.value();
        result.meta.type          = config.type.value();
        /* clang-format off */
        eig::log->trace("- {:<30} = {}" ,"eigvals_found",  result.meta.eigvals_found);
        eig::log->trace("- {:<30} = {}" ,"eigvecsR_found", result.meta.eigvecsR_found);
        eig::log->trace("- {:<30} = {}" ,"eigvecsL_found", result.meta.eigvecsL_found);
        eig::log->trace("- {:<30} = {}" ,"iter",           result.meta.iter);
        eig::log->trace("- {:<30} = {}" ,"n",              result.meta.n);
        eig::log->trace("- {:<30} = {}" ,"nev",            result.meta.nev);
        eig::log->trace("- {:<30} = {}" ,"nev_converged",  result.meta.nev_converged);
        eig::log->trace("- {:<30} = {}" ,"ncv_used",       result.meta.ncv_used);
        eig::log->trace("- {:<30} = {}" ,"rows",           result.meta.rows);
        eig::log->trace("- {:<30} = {}" ,"cols",           result.meta.cols);
        eig::log->trace("- {:<30} = {}" ,"counter",        result.meta.counter);
        /* clang-format on */
        assert(result.meta.nev_converged >= result.meta.nev and "Not enough eigenvalues converged");
        assert(result.meta.eigvals_found and "Eigenvalues were not found");
        assert(solver.EigenvectorsFound() and "Eigenvectors were not found");

    } else {
        eig::log->trace("Finding eigenvalues");
        solver.FindEigenvalues();
        result.meta.eigvals_found = solver.EigenvaluesFound();
        result.meta.iter          = solver.GetIter();
        result.meta.n             = solver.GetN();
        result.meta.nev           = std::min(nev, static_cast<eig::size_type>(solver.GetNev()));
        result.meta.nev_converged = solver.ConvergedEigenvalues();
        result.meta.ncv_used      = solver.GetNcv();
        result.meta.rows          = solver.GetN();
        result.meta.cols          = result.meta.nev;
        result.meta.counter       = matrix.counter;
        result.meta.form          = config.form.value();
        result.meta.type          = config.type.value();
        /* clang-format off */
        eig::log->trace("- {:<30} = {}" ,"eigvals_found",  result.meta.eigvals_found);
        eig::log->trace("- {:<30} = {}" ,"eigvecsR_found", result.meta.eigvecsR_found);
        eig::log->trace("- {:<30} = {}" ,"eigvecsL_found", result.meta.eigvecsL_found);
        eig::log->trace("- {:<30} = {}" ,"iter",           result.meta.iter);
        eig::log->trace("- {:<30} = {}" ,"n",              result.meta.n);
        eig::log->trace("- {:<30} = {}" ,"nev",            result.meta.nev);
        eig::log->trace("- {:<30} = {}" ,"nev_converged",  result.meta.nev_converged);
        eig::log->trace("- {:<30} = {}" ,"ncv_used",       result.meta.ncv_used);
        eig::log->trace("- {:<30} = {}" ,"rows",           result.meta.rows);
        eig::log->trace("- {:<30} = {}" ,"cols",           result.meta.cols);
        eig::log->trace("- {:<30} = {}" ,"counter",        result.meta.counter);
        /* clang-format on */

        assert(result.meta.nev_converged >= result.meta.nev and "Not enough eigenvalues converged");
        assert(result.meta.eigvals_found and "Eigenvalues were not found");
    }
}

template class eig::arpack_solver<MatrixProductDense<real>>;
template class eig::arpack_solver<MatrixProductDense<cplx>>;
template class eig::arpack_solver<MatrixProductSparse<real>>;
template class eig::arpack_solver<MatrixProductSparse<cplx>>;
template class eig::arpack_solver<MatrixProductHamiltonian<real>>;
template class eig::arpack_solver<MatrixProductHamiltonian<cplx>>;
