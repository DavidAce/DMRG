#include "solver_arpack.h"
#include "../log.h"
#include "../matvec/matvec_dense.h"
#include "../matvec/matvec_mpo.h"
#include "../matvec/matvec_sparse.h"
#include <algorithm>
#include <general/sfinae.h>
#include <tid/tid.h>

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
eig::solver_arpack<MatrixType>::solver_arpack(MatrixType &matrix_, eig::settings &config_, eig::solution &result_, Scalar *residual_)
    : matrix(matrix_), config(config_), result(result_), residual(residual_) {
    t_tot        = std::make_unique<tid::ur>("total");
    t_mul        = std::make_unique<tid::ur>("matvec");
    t_fnd        = std::make_unique<tid::ur>("find_eigenpair");
    t_pre        = std::make_unique<tid::ur>("preparation");
    nev_internal = std::clamp<int>(static_cast<int>(config.maxNev.value()), 1, matrix.rows() / 2);
    ncv_internal = std::clamp<int>(static_cast<int>(config.maxNcv.value()), static_cast<int>(config.maxNev.value()) + 1, matrix.rows());
}

template<typename MatrixType>
void eig::solver_arpack<MatrixType>::eigs() {
    result.reset();
    nev_internal = std::clamp<int>(static_cast<int>(config.maxNev.value()), 1, matrix.rows() / 2);
    ncv_internal = std::clamp<int>(static_cast<int>(config.maxNcv.value()), static_cast<int>(config.maxNev.value()) + 1, matrix.rows());

    config.checkRitz();
    matrix.set_mode(config.form.value());
    matrix.set_side(config.side.value());

    // Dispatch to symmetric or nonsymmetric. If complex, there's only a nonsymmetric option available.
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        eigs_comp_rc();
    } else {
        if(config.form == Form::SYMM)
            eigs_sym_rc();
        else
            eigs_nsym_rc();
    }
}

template<typename MatrixType>
void eig::solver_arpack<MatrixType>::eigs_sym() {
    if constexpr(std::is_same<Scalar, double>::value) {
        if(config.form != Form::SYMM) throw std::runtime_error("ERROR: config not SYMMETRIC");
        if(matrix.get_form() != Form::SYMM) throw std::runtime_error("ERROR: matrix not SYMMETRIC");
        ARSymStdEig<double, MatrixType> solver(matrix.rows(), nev_internal, &matrix, &MatrixType::MultAx, config.get_ritz_string().data(), ncv_internal,
                                               config.tol.value(), static_cast<int>(config.maxIter.value()), residual);

        find_solution(solver, config.maxNev.value());
        copy_solution(solver);

    } else {
        eig::log->critical("Called eigs_sym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
        throw std::runtime_error("Called eigs_sym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}

template<typename MatrixType>
void eig::solver_arpack<MatrixType>::eigs_sym_rc() {
    if constexpr(std::is_same<Scalar, double>::value) {
        if(config.form != Form::SYMM) throw std::runtime_error("ERROR: config not SYMMETRIC");
        if(matrix.get_form() != Form::SYMM) throw std::runtime_error("ERROR: matrix not SYMMETRIC");
        ARrcSymStdEig<double> solver(matrix.rows(), nev_internal, config.get_ritz_string().data(), ncv_internal, config.tol.value(),
                                     static_cast<int>(config.maxIter.value()), residual, true);
        find_solution_rc(solver);
        copy_solution(solver);
    } else {
        eig::log->critical("Called eigs_sym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
        throw std::runtime_error("Called eigs_sym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}

template<typename MatrixType>
void eig::solver_arpack<MatrixType>::eigs_nsym() {
    if constexpr(std::is_same<Scalar, double>::value) {
        if(config.form != Form::NSYM) throw std::runtime_error("ERROR: config not NSYM");
        if(matrix.get_form() != Form::NSYM) throw std::runtime_error("ERROR: matrix not NSYM");
        if(nev_internal == 1) { nev_internal++; }

        ARNonSymStdEig<double, MatrixType> solver(matrix.rows(), nev_internal, &matrix, &MatrixType::MultAx, config.get_ritz_string().data(), ncv_internal,
                                                  config.tol.value(), config.maxIter.value(), residual);
        find_solution(solver, config.maxNev.value());
        copy_solution(solver);
    } else {
        throw std::runtime_error("Called eigs_nsym() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}

template<typename MatrixType>
void eig::solver_arpack<MatrixType>::eigs_nsym_rc() {
    if constexpr(std::is_same<Scalar, double>::value) {
        if(config.form != Form::NSYM) throw std::runtime_error("ERROR: config not NSYM");
        if(matrix.get_form() != Form::NSYM) throw std::runtime_error("ERROR: matrix not NSYM");
        //        if(nev_internal == 1) { nev_internal++; }
        ARrcNonSymStdEig<double> solver(matrix.rows(), nev_internal, config.get_ritz_string().data(), ncv_internal, config.tol.value(),
                                        static_cast<int>(config.maxIter.value()), residual, true);

        find_solution_rc(solver);
        copy_solution(solver);
    } else {
        throw std::runtime_error("Called eigs_nsym_rc() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}

template<typename MatrixType>
void eig::solver_arpack<MatrixType>::eigs_comp() {
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        ARCompStdEig<double, MatrixType> solver(matrix.rows(), nev_internal, &matrix, &MatrixType::MultAx, config.get_ritz_string().data(), ncv_internal,
                                                config.tol.value(), config.maxIter.value(), residual);
        find_solution(solver, config.maxNev.value());
        copy_solution(solver);
    } else {
        throw std::runtime_error("Called eigs_comp() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}

template<typename MatrixType>
void eig::solver_arpack<MatrixType>::eigs_comp_rc() {
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        ARrcCompStdEig<double> solver(matrix.rows(), nev_internal, config.get_ritz_string().data(), ncv_internal, config.tol.value(),
                                      static_cast<int>(config.maxIter.value()), residual, true);
        find_solution_rc(solver);
        copy_solution(solver);
    } else {
        throw std::runtime_error("Called eigs_comp_rc() with wrong type: " + std::string(tc::type_name<MatrixType>()));
    }
}

template<typename MatrixType>
template<typename Derived>
void eig::solver_arpack<MatrixType>::find_solution(Derived &solver, eig::size_type nev) {
    using ShiftType  = decltype(solver.GetShift());
    auto t_tot_token = t_tot->tic_token();

    // Start by preparing the matrix for solving. Apply shifts, do inversion/factorization and then compress
    auto t_pre_token = t_pre->tic_token();
    if constexpr(MatrixType::can_shift) {
        if(config.sigma) {
            eig::log->trace("Setting shift with sigma = {}", std::real(config.sigma.value()));
            matrix.set_shift(config.sigma.value());
            if constexpr(MatrixType::can_shift_invert) {
                if(config.shift_invert == Shinv::ON) {
                    eig::log->trace("Enabling shift-invert mode");
                    matrix.FactorOP();
                    if constexpr(std::is_same_v<ShiftType, cplx>) solver.SetShiftInvertMode(config.sigma.value(), &matrix, &MatrixType::MultOPv);
                    if constexpr(std::is_same_v<ShiftType, real>) solver.SetShiftInvertMode(std::real(config.sigma.value()), &matrix, &MatrixType::MultOPv);
                }
            } else if(config.shift_invert == Shinv::ON)
                throw std::runtime_error("Tried to shift-invert an incompatible matrix");
        }
    } else if(config.sigma) {
        throw std::runtime_error("Tried to apply shift on an incompatible matrix");
    }

    if constexpr(MatrixType::can_compress) {
        if(config.compress and config.compress.value()) matrix.compress();
    }
    t_pre_token.toc();

    auto t_mul_token = t_mul->tic_token();
    solver.FindArnoldiBasis();
    result.meta.arnoldi_found = solver.ArnoldiBasisFound();
    t_mul_token.toc();
    auto t_fnd_token = t_fnd->tic_token();
    if(solver.ArnoldiBasisFound()) {
        if(config.compute_eigvecs) {
            eig::log->trace("Finding eigenvectors");
            solver.FindEigenvectors();
            if(not solver.EigenvaluesFound()) eig::log->warn("Eigenvalues were not found");
            if(not solver.EigenvectorsFound()) eig::log->warn("Eigenvectors were not found");
        } else {
            eig::log->trace("Finding eigenvalues");
            solver.FindEigenvalues();
            if(not solver.EigenvaluesFound()) eig::log->warn("Eigenvalues were not found");
        }
    }

    t_fnd_token.toc();
    t_tot_token.toc();
    if(config.side == eig::Side::L)
        result.meta.eigvecsL_found = solver.EigenvectorsFound();
    else
        result.meta.eigvecsR_found = solver.EigenvectorsFound(); // BOOL!

    result.meta.eigvals_found = solver.EigenvaluesFound(); // BOOL!
    result.meta.iter          = solver.GetIter();
    result.meta.n             = solver.GetN();
    result.meta.nev           = std::min(nev, static_cast<eig::size_type>(solver.GetNev()));
    result.meta.nev_converged = solver.ConvergedEigenvalues();
    result.meta.ncv           = solver.GetNcv();
    result.meta.rows          = solver.GetN();
    result.meta.cols          = result.meta.nev;
    result.meta.counter       = matrix.counter;
    result.meta.form          = config.form.value();
    result.meta.type          = config.type.value();
    result.meta.ritz          = solver.GetWhich();
    result.meta.sigma         = solver.GetShift();
    result.meta.time_total    = t_tot->get_time();
    result.meta.time_matprod  = t_mul->get_time() + t_fnd->get_time();
    result.meta.time_prep     = t_pre->get_time();

    /* clang-format off */
    eig::log->trace("- {:<30} = {}"     ,"arnoldi_found",  result.meta.arnoldi_found);
    eig::log->trace("- {:<30} = {}"     ,"eigvals_found",  result.meta.eigvals_found);
    eig::log->trace("- {:<30} = {}"     ,"eigvecsR_found", result.meta.eigvecsR_found);
    eig::log->trace("- {:<30} = {}"     ,"eigvecsL_found", result.meta.eigvecsL_found);
    eig::log->trace("- {:<30} = {}"     ,"iter",           result.meta.iter);
    eig::log->trace("- {:<30} = {}"     ,"counter",        result.meta.counter);
    eig::log->trace("- {:<30} = {}"     ,"rows",           result.meta.rows);
    eig::log->trace("- {:<30} = {}"     ,"cols",           result.meta.cols);
    eig::log->trace("- {:<30} = {}"     ,"n",              result.meta.n);
    eig::log->trace("- {:<30} = {}"     ,"nev",            result.meta.nev);
    eig::log->trace("- {:<30} = {}"     ,"nev_converged",  result.meta.nev_converged);
    eig::log->trace("- {:<30} = {}"     ,"ncv",            result.meta.ncv);
    eig::log->trace("- {:<30} = {}"     ,"ritz",           result.meta.ritz);
    eig::log->trace("- {:<30} = {}"     ,"sigma",          result.meta.sigma);
    eig::log->trace("- {:<30} = {:.3f}s","time total",     result.meta.time_total);
    eig::log->trace("- {:<30} = {:.3f}s","time matprod",   result.meta.time_matprod);
    eig::log->trace("- {:<30} = {:.3f}s","time prep",      result.meta.time_prep);
    /* clang-format on */

    if(result.meta.nev_converged < result.meta.nev) eig::log->warn("Not enough eigenvalues converged");
}

template<typename MatrixType>
template<typename Derived>
void eig::solver_arpack<MatrixType>::find_solution_rc(Derived &solver) {
    using ShiftType  = decltype(solver.GetShift());
    auto t_tot_token = t_tot->tic_token();

    // Start by prepare the matrix for solving. Apply shifts, do inversion/factorization and then compress
    auto t_pre_token = t_pre->tic_token();
    if constexpr(MatrixType::can_shift) {
        if(config.sigma) {
            eig::log->trace("Setting shift with sigma = {}", std::real(config.sigma.value()));
            matrix.set_shift(config.sigma.value());
            if constexpr(MatrixType::can_shift_invert) {
                if(config.shift_invert == Shinv::ON) {
                    eig::log->trace("Enabling shift-invert mode");
                    matrix.FactorOP();
                    if constexpr(std::is_same_v<ShiftType, cplx>) solver.SetShiftInvertMode(config.sigma.value());
                    if constexpr(std::is_same_v<ShiftType, real>) solver.SetShiftInvertMode(std::real(config.sigma.value()));
                }
            } else if(config.shift_invert == Shinv::ON)
                throw std::runtime_error("Tried to shift-invert an incompatible matrix");
        }
    } else if(config.sigma) {
        throw std::runtime_error("Tried to apply shift on an incompatible matrix");
    }

    if constexpr(MatrixType::can_compress) {
        if(config.compress and config.compress.value()) matrix.compress();
    }
    t_pre_token.toc();

    // Generate an Arnoldi basis
    int step = 0;
    int nops = 0;
    int iter = 0;
    while(not solver.ArnoldiBasisFound()) {
        if(step == 0) {
            bool        iter_lim = not config.iter_ncv_x.empty() and iter >= config.iter_ncv_x.back();
            bool        time_lim = not config.time_tol_x10.empty() and t_tot->get_time() >= config.time_tol_x10.back();
            std::string msg;
            if(time_lim) {
                // Do this when the iterations are taking a really long time. Remember, the next site may be easier to optimize,
                // and when we eventually return to this site with updated environments, it may be not be as hard to optimize anymore.
                auto tol = solver.GetTol();
                msg.append(fmt::format("| tol {:<.4e} -> {:<.4e}", tol, 1e1 * tol));
                solver.ChangeTol(1e1 * tol); // Does not restart iterations
                config.time_tol_x10.pop_back();
            }
            if(iter_lim) {
                // The convergence rate/iteration is seems poor, so we increase ncv and restart the arnoldi iterations
                auto ncv  = solver.GetNcv();
                auto ncvx = std::clamp(config.ncv_x_factor * ncv, 2 * solver.GetNev() + 1, std::max(2 * solver.GetNev() + 1, solver.GetN() / 8));
                msg.append(fmt::format("| ncv {:<4} -> {:<4}", ncv, ncvx));
                if(ncv != ncvx) solver.ChangeNcv(ncvx); // This will restart the arnoldi iterations
                config.iter_ncv_x.pop_back();
            }
            static double last_log_time = 0;
            if(last_log_time > t_tot->get_time()) last_log_time = 0; // If this function has been run before
            if(std::abs(t_tot->get_time() - last_log_time) > 5.0) {
                eig::log->info("iter {:<4} | ops {:<5} | time {:8.2f} s | dt {:8.2f} ms/op {}", iter, nops, t_tot->get_time(),
                               t_mul->get_last_interval() / nops * 1000, msg);
                last_log_time = t_tot->get_time();
            }
        }
        if(config.maxTime and config.maxTime.value() <= t_tot->get_time()) break;

        solver.TakeStep();
        if(std::abs(solver.GetIdo()) == 1) {
            t_mul->tic();
            if constexpr(MatrixType::can_shift_invert) {
                if(matrix.isReadyFactorOp())
                    matrix.MultOPv(solver.GetVector(), solver.PutVector());
                else
                    matrix.MultAx(solver.GetVector(), solver.PutVector());
            } else
                matrix.MultAx(solver.GetVector(), solver.PutVector());
            t_mul->toc();
        }
        step++;
        nops++;
        if((iter == 0 and step > solver.GetNev()) or (step > solver.GetNcv() - solver.GetNev())) {
            step = 0;
            iter++;
        }
    }

    eig::log->info("iter {:<4} | nops {:<5} | time {:8.2f} s | dt {:8.2f} ms | actual iters {}", iter, nops, t_tot->get_time(),
                   t_tot->get_time() / solver.GetIter() * 1000, solver.GetIter());
    result.meta.arnoldi_found = solver.ArnoldiBasisFound(); // Copy the value here because solver.FindEigenv...() will set BasisOk=false later
    if(not solver.ArnoldiBasisFound()) eig::log->warn("Arnoldi basis was not found");

    // Find the eigenvectors/eigenvalues
    auto t_fnd_token = t_fnd->tic_token();
    if(solver.ArnoldiBasisFound()) {
        if(config.compute_eigvecs) {
            solver.FindEigenvectors();
            if(not solver.EigenvaluesFound()) eig::log->warn("Eigenvalues were not found");
            if(not solver.EigenvectorsFound()) eig::log->warn("Eigenvectors were not found");
        } else {
            eig::log->trace("Finding eigenvalues");
            solver.FindEigenvalues();
            if(not solver.EigenvaluesFound()) eig::log->warn("Eigenvalues were not found");
        }
    }

    t_fnd_token.toc();
    t_tot_token.toc();

    if(config.side == eig::Side::L)
        result.meta.eigvecsL_found = solver.EigenvectorsFound();
    else
        result.meta.eigvecsR_found = solver.EigenvectorsFound(); // BOOL!
    result.meta.eigvals_found = solver.EigenvaluesFound();       // BOOL!
    result.meta.iter          = solver.GetIter();
    result.meta.n             = solver.GetN();
    result.meta.nev           = solver.GetNev();
    result.meta.nev_converged = solver.ConvergedEigenvalues();
    result.meta.ncv           = solver.GetNcv();
    result.meta.tol           = solver.GetTol();
    result.meta.rows          = solver.GetN();
    result.meta.cols          = result.meta.nev;
    result.meta.counter       = matrix.counter;
    result.meta.form          = config.form.value();
    result.meta.type          = config.type.value();
    result.meta.ritz          = solver.GetWhich();
    result.meta.sigma         = (config.sigma ? config.sigma.value() : result.meta.sigma); // solver.GetShift() is only for shift-invert
    result.meta.time_total    = t_tot->get_time();
    result.meta.time_matprod  = t_mul->get_time() + t_fnd->get_time();
    result.meta.time_prep     = t_pre->get_time();

    /* clang-format off */
    eig::log->trace("Arpack finished");
    eig::log->trace("- {:<30} = {}"     ,"arnoldi_found",  result.meta.arnoldi_found);
    eig::log->trace("- {:<30} = {}"     ,"eigvals_found",  result.meta.eigvals_found);
    eig::log->trace("- {:<30} = {}"     ,"eigvecsR_found", result.meta.eigvecsR_found);
    eig::log->trace("- {:<30} = {}"     ,"eigvecsL_found", result.meta.eigvecsL_found);
    eig::log->trace("- {:<30} = {}"     ,"iter",           result.meta.iter);
    eig::log->trace("- {:<30} = {}"     ,"counter",        result.meta.counter);
    eig::log->trace("- {:<30} = {}"     ,"rows",           result.meta.rows);
    eig::log->trace("- {:<30} = {}"     ,"cols",           result.meta.cols);
    eig::log->trace("- {:<30} = {}"     ,"n",              result.meta.n);
    eig::log->trace("- {:<30} = {}"     ,"nev",            result.meta.nev);
    eig::log->trace("- {:<30} = {}"     ,"nev_converged",  result.meta.nev_converged);
    eig::log->trace("- {:<30} = {}"     ,"ncv",            result.meta.ncv);
    eig::log->trace("- {:<30} = {}"     ,"ritz",           result.meta.ritz);
    eig::log->trace("- {:<30} = {}"     ,"sigma",          result.meta.sigma);
    eig::log->trace("- {:<30} = {:.3f}s","time total",     result.meta.time_total);
    eig::log->trace("- {:<30} = {:.3f}s","time matprod",   result.meta.time_matprod);
    eig::log->trace("- {:<30} = {:.3f}s","time prep",      result.meta.time_prep);

    /* clang-format on */

    if(result.meta.nev_converged < result.meta.nev) eig::log->warn("Not enough eigenvalues converged");
}

template<typename MatrixType>
template<typename Derived>
void eig::solver_arpack<MatrixType>::copy_solution(Derived &solver) {
    // In this function we copy the solution "naively"
    eig::log->trace("Copying solution");
    using eigval_type                         = std::remove_pointer_t<decltype(solver.RawEigenvalues())>;
    using eigvec_type                         = std::remove_pointer_t<decltype(solver.RawEigenvectors())>;
    auto           eigvecsize                 = result.meta.rows * result.meta.cols;
    auto           eigvalsize                 = result.meta.cols;
    auto           eigvalsize_t               = static_cast<size_t>(eigvalsize);
    auto           eigvecsize_t               = static_cast<size_t>(eigvecsize);
    constexpr auto eigval_has_imag_separately = eig::sfinae::has_RawEigenvaluesImag_v<Derived>;
    constexpr auto eigval_is_cplx             = std::is_same_v<cplx, eigval_type>;
    constexpr auto eigval_is_real             = std::is_same_v<real, eigval_type>;
    constexpr auto eigvec_is_cplx             = std::is_same_v<cplx, eigvec_type>;
    constexpr auto eigvec_is_real             = std::is_same_v<real, eigvec_type>;

    if(result.meta.eigvals_found) {
        if constexpr(eigval_has_imag_separately) {
            eig::log->trace("Copying eigenvalues from separate real and imaginary buffers");
            result.eigvals_imag.resize(eigvalsize_t);
            result.eigvals_real.resize(eigvalsize_t);
            std::copy(solver.RawEigenvaluesImag(), solver.RawEigenvaluesImag() + eigvalsize, result.eigvals_imag.begin());
            std::copy(solver.RawEigenvaluesReal(), solver.RawEigenvaluesReal() + eigvalsize, result.eigvals_real.begin());
        }
        if constexpr(eigval_is_real) {
            if(not solver.EigenvaluesFound()) throw std::runtime_error("Eigenvalues were not found");
            eig::log->trace("Copying real eigenvalues");
            result.eigvals_real.resize(eigvalsize_t);
            std::copy(solver.RawEigenvalues(), solver.RawEigenvalues() + eigvalsize, result.eigvals_real.begin());
        }
        if constexpr(eigval_is_cplx) {
            if(not solver.EigenvaluesFound()) throw std::runtime_error("Eigenvalues were not found");
            eig::log->trace("Copying complex eigenvalues");
            result.eigvals_cplx.resize(eigvalsize_t);
            std::copy(solver.RawEigenvalues(), solver.RawEigenvalues() + eigvalsize, result.eigvals_cplx.begin());
        }
    }
    if(result.meta.eigvecsR_found or result.meta.eigvecsL_found) {
        if constexpr(eigvec_is_real) {
            if(not solver.EigenvectorsFound()) throw std::runtime_error("Eigenvectors were not found");
            eig::log->trace("Copying real eigenvectors");
            result.eigvecsR_real.resize(eigvecsize_t);
            std::copy(solver.RawEigenvectors(), solver.RawEigenvectors() + eigvecsize, result.eigvecsR_real.begin());
        }

        if constexpr(eigvec_is_cplx) {
            if(not solver.EigenvectorsFound()) throw std::runtime_error("Eigenvectors were not found");
            eig::log->trace("Copying complex eigenvectors");
            if(config.side == Side::L) {
                result.eigvecsL_cplx.resize(eigvecsize_t);
                std::copy(solver.RawEigenvectors(), solver.RawEigenvectors() + eigvecsize, result.eigvecsL_cplx.begin());
            } else {
                result.eigvecsR_cplx.resize(eigvecsize_t);
                std::copy(solver.RawEigenvectors(), solver.RawEigenvectors() + eigvecsize, result.eigvecsR_cplx.begin());
            }
        }
    }

    tid::get("arpack") += t_tot->get_time();
    tid::get("arpack.matvec") += t_mul->get_time();
    tid::get("arpack.findeig") += t_fnd->get_time();
    tid::get("arpack.prep") += t_pre->get_time();
}

template class eig::solver_arpack<MatVecDense<real>>;
template class eig::solver_arpack<MatVecDense<cplx>>;
template class eig::solver_arpack<MatVecSparse<real>>;
template class eig::solver_arpack<MatVecSparse<cplx>>;
template class eig::solver_arpack<MatVecMPO<real>>;
template class eig::solver_arpack<MatVecMPO<cplx>>;
