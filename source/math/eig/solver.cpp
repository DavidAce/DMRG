#include "solver.h"
#include "debug/exceptions.h"
#include "log.h"
#include "matvec/matvec_dense.h"
#include "matvec/matvec_mpo.h"
#include "matvec/matvec_mpos.h"
#include "matvec/matvec_sparse.h"
#include "matvec/matvec_zero.h"
#include "solver_arpack/solver_arpack.h"
#include "tid/tid.h"

int eig::getBasisSize(long L, int nev, std::optional<int> basisSize) {
    if(not basisSize.has_value() or basisSize.value() <= 0) { basisSize = nev * safe_cast<int>(std::ceil(std::log2(L))); }
    return std::clamp(basisSize.value(), 2 * nev, safe_cast<int>(L));
}

eig::solver::solver() {
    if(config.loglevel)
        eig::setLevel(config.loglevel.value());
    else
        eig::setLevel(2);
    eig::setTimeStamp();
    log = eig::log;
}

eig::solver::solver(const eig::settings &config_) : eig::solver::solver() { config = config_; }

void eig::solver::setLogLevel(size_t loglevel) {
    config.loglevel = loglevel;
    eig::setLevel(loglevel);
}

template<typename Scalar>
void eig::solver::subtract_phase(std::vector<Scalar> &eigvecs, size_type L, int nev)
// The solution to  the eigenvalue equation Av = l*v is determined up to a constant phase factor, i.e., if v
// is a solution, so is v*exp(i*theta). By computing the complex angle of the first element in v, one can then
// remove it from all other elements of v.
{
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        if(eigvecs.empty()) return;
        if(nev > 0) {
            for(int i = 0; i < nev; i++) {
                if(eigvecs[safe_cast<size_t>(i * L)].imag() == 0.0) { continue; }
                Scalar inv_phase     = Scalar(0.0, -1.0) * std::arg(eigvecs[safe_cast<size_t>(i * L)]);
                auto   begin         = eigvecs.begin() + i * L;
                auto   end           = begin + L;
                Scalar exp_inv_phase = std::exp(inv_phase);
                std::transform(begin, end, begin, [exp_inv_phase](std::complex<double> num) -> std::complex<double> { return (num * exp_inv_phase); });
                std::transform(begin, end, begin,
                               [](std::complex<double> num) -> std::complex<double> { return std::abs(num.imag()) > 1e-15 ? num : std::real(num); });
            }
        } else {
            throw std::logic_error("Subtract phase requires nev > 0");
        }
    }
}

template void eig::solver::subtract_phase(std::vector<cplx> &eigvecs, size_type L, int nev);
template void eig::solver::subtract_phase(std::vector<real> &eigvecs, size_type L, int nev);

void eig::solver::eig_init(Form form, Type type, Vecs compute_eigvecs, Dephase remove_phase) {
    eig::log->trace("eig init");
    result.reset();
    config.compute_eigvecs = config.compute_eigvecs.value_or(compute_eigvecs);
    config.remove_phase    = config.remove_phase.value_or(remove_phase);
    config.type            = config.type.value_or(type);
    config.form            = config.form.value_or(form);
    config.side            = config.side.value_or(Side::LR);
    config.storage         = config.storage.value_or(Storage::DENSE);
}

template<eig::Form form, typename Scalar>
void eig::solver::eig(Scalar *matrix, size_type L, Vecs compute_eigvecs_, Dephase remove_phase_) {
    static_assert(!std::is_const_v<Scalar>);
    auto t_eig = tid::tic_scope("eig");
    int  info  = 0;
    try {
        if constexpr(std::is_same_v<Scalar, real>) {
            eig_init(form, Type::REAL, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) {
                if(config.tag.empty()) config.tag = "dsyevd";
                info = dsyevd(matrix, L);
            } else if constexpr(form == Form::NSYM) {
                if(config.tag.empty()) config.tag = "dgeev";
                info = dgeev(matrix, L);
            }
        } else if constexpr(std::is_same_v<Scalar, cplx>) {
            eig_init(form, Type::CPLX, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) {
                if(config.tag.empty()) config.tag = "zheevd";
                info = zheevd(matrix, L);
                //                if(config.tag.empty()) config.tag = "zheev";
                //                info = zheev(matrix, L);
            } else if constexpr(form == Form::NSYM) {
                if(config.tag.empty()) config.tag = "zgeev";
                info = zgeev(matrix, L);
            }
        } else {
            throw except::runtime_error("Unknown type");
        }

    } catch(std::exception &ex) {
        eig::log->error("Eigenvalue solver failed: {}", ex.what());
        throw except::runtime_error("Eigenvalue solver Failed: {}", ex.what());
    }

    result.build_eigvals_cplx();
    result.build_eigvecs_cplx();

    if(info == 0 and config.remove_phase and config.remove_phase.value() == Dephase::OFF) {
        // The solution to  the eigenvalue equation Av = l*v is determined up to a constant phase factor, i.e., if v
        // is a solution, so is v*exp(i*theta). By computing the complex angle of the first element in v, one can then
        // remove it from all other elements of v.
        subtract_phase(result.eigvecsL_cplx, L, result.meta.nev);
        subtract_phase(result.eigvecsR_cplx, L, result.meta.nev);
    }
}

template void eig::solver::eig<eig::Form::SYMM>(real *matrix, size_type, Vecs, Dephase);
template void eig::solver::eig<eig::Form::NSYM>(real *matrix, size_type, Vecs, Dephase);
template void eig::solver::eig<eig::Form::SYMM>(cplx *matrix, size_type, Vecs, Dephase);
template void eig::solver::eig<eig::Form::NSYM>(cplx *matrix, size_type, Vecs, Dephase);

template<eig::Form form, typename Scalar>
void eig::solver::eig(Scalar *matrixA, Scalar *matrixB, size_type L, Vecs compute_eigvecs_, Dephase remove_phase_) {
    static_assert(!std::is_const_v<Scalar>);
    auto t_eig = tid::tic_scope("eig");
    int  info  = 0;
    try {
        if(matrixB == nullptr) {
            if constexpr(std::is_same_v<Scalar, real>) {
                eig_init(form, Type::REAL, compute_eigvecs_, remove_phase_);
                if constexpr(form == Form::SYMM) {
                    if(config.tag.empty()) config.tag = "dsyevd";
                    info = dsyevd(matrixA, L);

                } else if constexpr(form == Form::NSYM) {
                    if(config.tag.empty()) config.tag = "dgeev";
                    info = dgeev(matrixA, L);
                }
            } else if constexpr(std::is_same_v<Scalar, cplx>) {
                eig_init(form, Type::CPLX, compute_eigvecs_, remove_phase_);
                if constexpr(form == Form::SYMM) {
                    if(config.tag.empty()) config.tag = "zheevd";
                    info = zheevd(matrixA, L);
                    //                if(config.tag.empty()) config.tag = "zheev";
                    //                info = zheev(matrix, L);
                } else if constexpr(form == Form::NSYM) {
                    if(config.tag.empty()) config.tag = "zgeev";
                    info = zgeev(matrixA, L);
                }
            } else {
                throw except::runtime_error("Unknown type");
            }
        } else {
            throw except::logic_error("The generalized solvers have not been implemented");
            if constexpr(std::is_same_v<Scalar, real>) {
                eig_init(form, Type::REAL, compute_eigvecs_, remove_phase_);
                if constexpr(form == Form::SYMM) {
                    if(config.tag.empty()) config.tag = "dsygvd";
                    // info = dsygvd(matrixA, matrixB, L);

                } else if constexpr(form == Form::NSYM) {
                    if(config.tag.empty()) config.tag = "dggev";
                    // info = dggev(matrixA, matrixB, L);
                }
            } else if constexpr(std::is_same_v<Scalar, cplx>) {
                eig_init(form, Type::CPLX, compute_eigvecs_, remove_phase_);
                if constexpr(form == Form::SYMM) {
                    if(config.tag.empty()) config.tag = "zhegvd";
                    // info = zhegvd(matrixA, matrixB, L);
                    //                if(config.tag.empty()) config.tag = "zheev";
                    //                info = zheev(matrix, L);
                } else if constexpr(form == Form::NSYM) {
                    if(config.tag.empty()) config.tag = "zggev";
                    // info = zggev(matrixA, matrixB, L);
                }
            } else {
                throw except::runtime_error("Unknown type");
            }
        }

    } catch(std::exception &ex) {
        eig::log->error("Eigenvalue solver failed: {}", ex.what());
        throw except::runtime_error("Eigenvalue solver Failed: {}", ex.what());
    }

    result.build_eigvals_cplx();
    result.build_eigvecs_cplx();

    if(info == 0 and config.remove_phase and config.remove_phase.value() == Dephase::OFF) {
        // The solution to  the eigenvalue equation Av = l*v is determined up to a constant phase factor, i.e., if v
        // is a solution, so is v*exp(i*theta). By computing the complex angle of the first element in v, one can then
        // remove it from all other elements of v.
        subtract_phase(result.eigvecsL_cplx, L, result.meta.nev);
        subtract_phase(result.eigvecsR_cplx, L, result.meta.nev);
    }
}
template void eig::solver::eig<eig::Form::SYMM>(real *matrixA, real *matrixB, size_type, Vecs, Dephase);
template void eig::solver::eig<eig::Form::NSYM>(real *matrixA, real *matrixB, size_type, Vecs, Dephase);
template void eig::solver::eig<eig::Form::SYMM>(cplx *matrixA, cplx *matrixB, size_type, Vecs, Dephase);
template void eig::solver::eig<eig::Form::NSYM>(cplx *matrixA, cplx *matrixB, size_type, Vecs, Dephase);

template<eig::Form form, typename Scalar>
void eig::solver::eig(Scalar *matrix, size_type L, char range, int il, int iu, double vl, double vu, Vecs compute_eigvecs_, Dephase remove_phase_) {
    static_assert(form == Form::SYMM);
    static_assert(!std::is_const_v<Scalar>);
    auto t_eig = tid::tic_scope("eig");
    int  info  = 0;
    try {
        if constexpr(std::is_same_v<Scalar, real>) {
            eig_init(form, Type::REAL, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) {
                if(config.tag.empty()) config.tag = "dsyevr";
                info = dsyevr(matrix, L, range, il, iu, vl, vu);
            }
            if constexpr(form == Form::NSYM) throw std::logic_error("?sygvr not implemented");
        } else if constexpr(std::is_same_v<Scalar, cplx>) {
            eig_init(form, Type::CPLX, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) {
                if(config.tag.empty()) config.tag = "zheevr";
                info = zheevr(matrix, L, range, il, iu, vl, vu);
            }
            if constexpr(form == Form::NSYM) throw std::logic_error("?sygvx not implemented");
        } else {
            throw except::runtime_error("Unknown type");
        }

    } catch(std::exception &ex) {
        eig::log->error("Eigenvalue solver failed: {}", ex.what());
        throw except::runtime_error("Eigenvalue solver Failed: {}", ex.what());
    }

    result.build_eigvals_cplx();
    result.build_eigvecs_cplx();

    if(info == 0 and config.remove_phase and config.remove_phase.value() == Dephase::OFF) {
        // The solution to  the eigenvalue equation Av = l*v is determined up to a constant phase factor, i.e., if v
        // is a solution, so is v*exp(i*theta). By computing the complex angle of the first element in v, one can then
        // remove it from all other elements of v.
        subtract_phase(result.eigvecsL_cplx, L, result.meta.nev);
        subtract_phase(result.eigvecsR_cplx, L, result.meta.nev);
    }
}
template void eig::solver::eig<eig::Form::SYMM>(real *matrix, size_type, char, int, int, double, double, Vecs, Dephase);
template void eig::solver::eig<eig::Form::SYMM>(cplx *matrix, size_type, char, int, int, double, double, Vecs, Dephase);

template<eig::Form form, typename Scalar>
void eig::solver::eig(Scalar *matrixA, Scalar *matrixB, size_type L, char range, int il, int iu, double vl, double vu, Vecs compute_eigvecs_,
                      Dephase remove_phase_) {
    static_assert(form == Form::SYMM);
    static_assert(!std::is_const_v<Scalar>);
    auto t_eig = tid::tic_scope("eig");
    int  info  = 0;
    try {
        if(matrixB == nullptr) {
            if constexpr(std::is_same_v<Scalar, real>) {
                eig_init(form, Type::REAL, compute_eigvecs_, remove_phase_);
                if constexpr(form == Form::SYMM) {
                    if(config.tag.empty()) config.tag = "dsyevr";
                    info = dsyevr(matrixA, L, range, il, iu, vl, vu);
                }
                if constexpr(form == Form::NSYM) throw std::logic_error("?sygvr not implemented");
            } else if constexpr(std::is_same_v<Scalar, cplx>) {
                eig_init(form, Type::CPLX, compute_eigvecs_, remove_phase_);
                if constexpr(form == Form::SYMM) {
                    if(config.tag.empty()) config.tag = "zheevr";
                    info = zheevr(matrixA, L, range, il, iu, vl, vu);
                }
                if constexpr(form == Form::NSYM) throw std::logic_error("?sygvx not implemented");
            } else {
                throw except::runtime_error("Unknown type");
            }
        }else {
            if constexpr(std::is_same_v<Scalar, real>) {
                eig_init(form, Type::REAL, compute_eigvecs_, remove_phase_);
                if constexpr(form == Form::SYMM) {
                    if(config.tag.empty()) config.tag = "dsygvx";
                    info = dsygvx(matrixA, matrixB, L, range, il, iu, vl, vu);
                }
                if constexpr(form == Form::NSYM) throw std::logic_error("?sygvr not implemented");
            } else if constexpr(std::is_same_v<Scalar, cplx>) {
                eig_init(form, Type::CPLX, compute_eigvecs_, remove_phase_);
                if constexpr(form == Form::SYMM) {
                    if(config.tag.empty()) config.tag = "zhegvx";
                    throw except::logic_error("The generalized solvers have not been implemented (zhegvx)");
                    // info = zheevx(matrix, L, range, il, iu, vl, vu);
                }
                if constexpr(form == Form::NSYM) throw std::logic_error("?sygvx not implemented");
            } else {
                throw except::runtime_error("Unknown type");
            }

        }

    } catch(std::exception &ex) {
        eig::log->error("Eigenvalue solver failed: {}", ex.what());
        throw except::runtime_error("Eigenvalue solver Failed: {}", ex.what());
    }

    result.build_eigvals_cplx();
    result.build_eigvecs_cplx();

    if(info == 0 and config.remove_phase and config.remove_phase.value() == Dephase::OFF) {
        // The solution to  the eigenvalue equation Av = l*v is determined up to a constant phase factor, i.e., if v
        // is a solution, so is v*exp(i*theta). By computing the complex angle of the first element in v, one can then
        // remove it from all other elements of v.
        subtract_phase(result.eigvecsL_cplx, L, result.meta.nev);
        subtract_phase(result.eigvecsR_cplx, L, result.meta.nev);
    }
}
template void eig::solver::eig<eig::Form::SYMM>(real *matrixA, real *matrixB, size_type, char, int, int, double, double, Vecs, Dephase);
template void eig::solver::eig<eig::Form::SYMM>(cplx *matrixA, cplx *matrixB, size_type, char, int, int, double, double, Vecs, Dephase);

template<typename Scalar>
void eig::solver::eigs_init(size_type L, int nev, int ncv, Ritz ritz, Form form, Type type, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                            Storage storage, Vecs compute_eigvecs, Dephase remove_phase, Scalar *residual, Lib lib) {
    if(config.loglevel) eig::setLevel(config.loglevel.value());
    eig::log->trace("eigs init");
    result.reset();
    config.lib = config.lib.value_or(lib);

    // Precision settings which are overridden manually
    config.tol     = config.tol.value_or(1e-12);
    config.maxIter = config.maxIter.value_or(1000);

    // Other settings that we pass on each invocation of eigs
    config.compute_eigvecs = config.compute_eigvecs.value_or(compute_eigvecs);
    config.remove_phase    = config.remove_phase.value_or(remove_phase);
    config.maxNev          = config.maxNev.value_or(nev);
    config.maxNcv          = config.maxNcv.value_or(ncv);
    config.sigma           = config.sigma.value_or(sigma.value_or(cplx(0.0, 0.0)));
    config.shift_invert    = config.shift_invert.value_or(shift_invert);
    config.type            = config.type.value_or(type);
    config.form            = config.form.value_or(form);
    config.ritz            = config.ritz.value_or(ritz);
    config.side            = config.side.value_or(side);
    config.storage         = config.storage.value_or(storage);

    if(config.initial_guess.empty()) config.initial_guess.push_back({residual, 0});
    config.maxNev.value() = std::clamp(config.maxNev.value(), 1, safe_cast<int>(L));
    if(config.form == Form::NSYM) {
        if(config.maxNev.value() == 1) { config.maxNev = 2; }
    }

    config.maxNcv = getBasisSize(L, config.maxNev.value(), config.maxNcv); // Adjust ncv if <= 0 (autoselect) <= 2nev(clamp) or > L (clamp)
    assert(1 <= config.maxNev.value() and config.maxNev.value() <= L);
    assert(config.maxNev <= config.maxNcv.value() and config.maxNcv.value() <= L);

    if(config.shift_invert == Shinv::ON and not config.sigma) throw std::runtime_error("Sigma must be set to use shift-invert mode");
    config.checkRitz();

    if(not config.logTime) {
        if(eig::log->level() == spdlog::level::trace) config.logTime = 10.0;
        if(eig::log->level() == spdlog::level::debug) config.logTime = 60.0;
        if(eig::log->level() >= spdlog::level::info) config.logTime = 60.0 * 10;
    }
    if(not config.logIter) {
        if(eig::log->level() == spdlog::level::trace) config.logIter = 100;
        if(eig::log->level() == spdlog::level::debug) config.logIter = 1000;
        if(eig::log->level() >= spdlog::level::info) config.logIter = 5000;
    }
}
template void eig::solver::eigs_init(size_type L, int nev, int ncv, Ritz ritz, Form form, Type type, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                                     Storage storage, Vecs compute_eigvecs, Dephase remove_phase, real *residual, Lib lib);
template void eig::solver::eigs_init(size_type L, int nev, int ncv, Ritz ritz, Form form, Type type, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                                     Storage storage, Vecs compute_eigvecs, Dephase remove_phase, cplx *residual, Lib lib);

template<typename MatrixProductType>
void eig::solver::set_default_config(const MatrixProductType &matrix) {
    if(config.loglevel) eig::setLevel(config.loglevel.value());
    eig::log->trace("eigs init");

    if(not config.lib) config.lib = eig::Lib::PRIMME;

    // Precision settings which are overridden manually
    config.tol     = config.tol.value_or(1e-12);
    config.maxIter = config.maxIter.value_or(2000);

    // Other settings that we pass on each invocation of eigs
    config.compute_eigvecs = config.compute_eigvecs.value_or(eig::Vecs::OFF);
    config.remove_phase    = config.remove_phase.value_or(eig::Dephase::OFF);
    config.maxNev          = config.maxNev.value_or(1);
    config.maxNcv          = config.maxNcv.value_or(-1);
    config.shift_invert    = config.shift_invert.value_or(Shinv::OFF);
    config.type            = config.type.value_or(matrix.get_type());
    config.form            = config.form.value_or(matrix.get_form());
    config.side            = config.side.value_or(matrix.get_side());
    config.ritz            = config.ritz.value_or(eig::Ritz::SA);
    config.storage         = config.storage.value_or(MatrixProductType::storage);

    config.maxNev.value() = std::clamp(config.maxNev.value(), 1, safe_cast<int>(matrix.rows()));
    if(config.form == Form::NSYM) {
        if(config.maxNev.value() == 1) { config.maxNev = 2; }
    }

    config.maxNcv = getBasisSize(matrix.rows(), config.maxNev.value(), config.maxNcv); // Adjust ncv if <= 0 (autoselect) <= 2nev(clamp) or > L (clamp)
    assert(1 <= config.maxNev.value() and config.maxNev.value() <= matrix.rows());
    assert(config.maxNev <= config.maxNcv.value() and config.maxNcv.value() <= matrix.rows());

    if(config.shift_invert == Shinv::ON and not config.sigma) throw std::runtime_error("Sigma must be set to use shift-invert mode");
    config.checkRitz();

    if(not config.logTime) {
        if(eig::log->level() == spdlog::level::trace) config.logTime = 10.0;
        if(eig::log->level() == spdlog::level::debug) config.logTime = 60.0;
        if(eig::log->level() >= spdlog::level::info) config.logTime = 60.0 * 10;
    }
    if(not config.logIter) {
        if(eig::log->level() == spdlog::level::trace) config.logIter = 100;
        if(eig::log->level() == spdlog::level::debug) config.logIter = 1000;
        if(eig::log->level() == spdlog::level::info) config.logIter = 5000;
    }
}

template void eig::solver::set_default_config(const MatVecDense<real> &matrix);
template void eig::solver::set_default_config(const MatVecDense<cplx> &matrix);
template void eig::solver::set_default_config(const MatVecSparse<real> &matrix);
template void eig::solver::set_default_config(const MatVecSparse<cplx> &matrix);
template void eig::solver::set_default_config(const MatVecMPO<real> &matrix);
template void eig::solver::set_default_config(const MatVecMPO<cplx> &matrix);
template void eig::solver::set_default_config(const MatVecMPOS<real> &matrix);
template void eig::solver::set_default_config(const MatVecMPOS<cplx> &matrix);
template void eig::solver::set_default_config(const MatVecZero<real> &matrix);
template void eig::solver::set_default_config(const MatVecZero<cplx> &matrix);
template<typename Scalar, eig::Storage storage>
void eig::solver::eigs(const Scalar *matrix, size_type L, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                       Vecs compute_eigvecs, Dephase remove_phase, Scalar *residual) {
    auto t_eig = tid::tic_scope("eig");
    static_assert(storage != Storage::MPS and "Can't build MatVecMPO from a pointer");
    if constexpr(storage == Storage::DENSE) {
        auto matrix_dense = MatVecDense<Scalar>(matrix, L, true, form, side);
        eigs(matrix_dense, nev, ncv, ritz, form, side, sigma, shift_invert, compute_eigvecs, remove_phase, residual);
    } else if constexpr(storage == Storage::SPARSE) {
        auto matrix_sparse = MatVecSparse<Scalar, false>(matrix, L, true);
        eigs(matrix_sparse, nev, ncv, ritz, form, side, sigma, shift_invert, compute_eigvecs, remove_phase, residual);
    }
}

template void eig::solver::eigs(const real *matrix, size_type L, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma,
                                Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, real *residual);

template void eig::solver::eigs(const cplx *matrix, size_type L, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma,
                                Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, cplx *residual);

template<typename MatrixProductType>
void eig::solver::eigs(MatrixProductType &matrix, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                       Vecs compute_eigvecs, Dephase remove_phase, typename MatrixProductType::Scalar *residual) {
    auto t_eig   = tid::tic_scope("eig");
    using Scalar = typename MatrixProductType::Scalar;
    Type type;
    if constexpr(std::is_same_v<Scalar, real>)
        type = Type::REAL;
    else if constexpr(std::is_same_v<Scalar, cplx>)
        type = Type::CPLX;
    else
        throw std::runtime_error("Unsupported type");
    eigs_init(matrix.rows(), nev, ncv, ritz, form, type, side, sigma, shift_invert, matrix.storage, compute_eigvecs, remove_phase, residual);
    switch(config.lib.value()) {
        case Lib::ARPACK: {
            solver_arpack<MatrixProductType> solver(matrix, config, result);
            solver.eigs();
            break;
        }
        case Lib::PRIMME: {
            eigs_primme(matrix);
            break;
        }
    }
}

template void eig::solver::eigs(MatVecDense<real> &matrix, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                                Vecs compute_eigvecs, Dephase remove_phase, real *residual);

template void eig::solver::eigs(MatVecDense<cplx> &matrix, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                                Vecs compute_eigvecs, Dephase remove_phase, cplx *residual);

template void eig::solver::eigs(MatVecSparse<real> &matrix, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                                Vecs compute_eigvecs, Dephase remove_phase, real *residual);

template void eig::solver::eigs(MatVecSparse<cplx> &matrix, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                                Vecs compute_eigvecs, Dephase remove_phase, cplx *residual);

template void eig::solver::eigs(MatVecMPO<real> &matrix, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                                Vecs compute_eigvecs, Dephase remove_phase, real *residual);

template void eig::solver::eigs(MatVecMPO<cplx> &matrix, int nev, int ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                                Vecs compute_eigvecs, Dephase remove_phase, cplx *residual);

template<typename MatrixProductType>
void eig::solver::eigs(MatrixProductType &matrix) {
    auto t_eig = tid::tic_scope("eig");
    result.reset();
    set_default_config(matrix);
    switch(config.lib.value()) {
        case Lib::ARPACK: {
            config.tag = fmt::format("{}{}arpack", config.tag, config.tag.empty() ? "" : " ");
            solver_arpack<MatrixProductType> solver(matrix, config, result);
            solver.eigs();
            break;
        }
        case Lib::PRIMME: {
            config.tag = fmt::format("{}{}primme", config.tag, config.tag.empty() ? "" : " ");
            eigs_primme(matrix);
            break;
        }
    }
}

template void eig::solver::eigs(MatVecDense<real> &matrix);
template void eig::solver::eigs(MatVecDense<cplx> &matrix);
template void eig::solver::eigs(MatVecSparse<real> &matrix);
template void eig::solver::eigs(MatVecSparse<cplx> &matrix);
template void eig::solver::eigs(MatVecMPO<real> &matrix);
template void eig::solver::eigs(MatVecMPO<cplx> &matrix);
template void eig::solver::eigs(MatVecMPOS<real> &matrix);
template void eig::solver::eigs(MatVecMPOS<cplx> &matrix);
template void eig::solver::eigs(MatVecZero<real> &matrix);
template void eig::solver::eigs(MatVecZero<cplx> &matrix);