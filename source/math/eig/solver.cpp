#include "solver.h"
#include "debug/exceptions.h"
#include "log.h"
#include "matvec/matvec_dense.h"
#include "matvec/matvec_mpo.h"
#include "matvec/matvec_sparse.h"
#include "solver_arpack/solver_arpack.h"
#include <tid/tid.h>

eig::solver::solver() {
    if(config.loglevel)
        eig::setLevel(config.loglevel.value());
    else
        eig::setLevel(2);
    eig::setTimeStamp();
    log = eig::log;
}

void eig::solver::setLogLevel(size_t loglevel) {
    config.loglevel = loglevel;
    eig::setLevel(loglevel);
}

template<typename Scalar>
void eig::solver::subtract_phase(std::vector<Scalar> &eigvecs, size_type L, size_type nev)
// The solution to  the eigenvalue equation Av = l*v is determined up to a constant phase factor, i.e., if v
// is a solution, so is v*exp(i*theta). By computing the complex angle of the first element in v, one can then
// remove it from all other elements of v.
{
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        if(eigvecs.empty()) return;
        if(nev > 0) {
            for(size_type i = 0; i < nev; i++) {
                if(eigvecs[static_cast<size_t>(i * L)].imag() == 0.0) { continue; }
                Scalar inv_phase     = Scalar(0.0, -1.0) * std::arg(eigvecs[static_cast<size_t>(i * L)]);
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

template void eig::solver::subtract_phase(std::vector<eig::cplx> &eigvecs, size_type L, size_type nev);
template void eig::solver::subtract_phase(std::vector<eig::real> &eigvecs, size_type L, size_type nev);

void eig::solver::eig_init(Form form, Type type, Vecs compute_eigvecs, Dephase remove_phase) {
    eig::log->trace("eig init");
    result.reset();
    /* clang-format off */
    if(not config.compute_eigvecs) config.compute_eigvecs  = compute_eigvecs;
    if(not config.remove_phase)    config.remove_phase     = remove_phase;
    if(not config.type)            config.type             = type;
    if(not config.form)            config.form             = form;
    if(not config.side)            config.side             = Side::LR;
    if(not config.storage)         config.storage          = Storage::DENSE;
//    config.checkRitz(); // Not checked on eig
    /* clang-format on */
}

template<eig::Form form, typename Scalar>
void eig::solver::eig(const Scalar *matrix, size_type L, Vecs compute_eigvecs_, Dephase remove_phase_) {
    auto t_eig = tid::tic_scope("eig");
    int  info  = 0;
    try {
        if constexpr(std::is_same_v<Scalar, real>) {
            eig_init(form, Type::REAL, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) {
                if(config.tag.empty()) config.tag = "dsyevd";
                info = dsyevd(matrix, L);
            }
            if constexpr(form == Form::NSYM) {
                if(config.tag.empty()) config.tag = "dgeev";
                info = dgeev(matrix, L);
            }
        } else if constexpr(std::is_same_v<Scalar, cplx>) {
            eig_init(form, Type::CPLX, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) {
                if(config.tag.empty()) config.tag = "zheevd";
                info = zheevd(matrix, L);
            }
            if constexpr(form == Form::NSYM) {
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

template void eig::solver::eig<eig::Form::SYMM>(const real *matrix, size_type, Vecs, Dephase);
template void eig::solver::eig<eig::Form::NSYM>(const real *matrix, size_type, Vecs, Dephase);
template void eig::solver::eig<eig::Form::SYMM>(const cplx *matrix, size_type, Vecs, Dephase);
template void eig::solver::eig<eig::Form::NSYM>(const cplx *matrix, size_type, Vecs, Dephase);

template<eig::Form form, typename Scalar>
void eig::solver::eig(const Scalar *matrix, size_type L, char range, int il, int iu, double vl, double vu, int m, Vecs compute_eigvecs_,
                      Dephase remove_phase_) {
    auto t_eig = tid::tic_scope("eig");
    int  info  = 0;
    try {
        if constexpr(std::is_same_v<Scalar, real>) {
            eig_init(form, Type::REAL, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) {
                if(config.tag.empty()) config.tag = "dsyevr";
                info = dsyevr(matrix, L, range, il, iu, vl, vu, m);
            }
            if constexpr(form == Form::NSYM) throw std::logic_error("?sygvx not implemented");
        } else if constexpr(std::is_same_v<Scalar, cplx>) {
            eig_init(form, Type::CPLX, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) throw std::logic_error("zheevx not implemented");
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
template void eig::solver::eig<eig::Form::SYMM>(const real *matrix, size_type, char, int, int, double, double, int, Vecs, Dephase);
template void eig::solver::eig<eig::Form::SYMM>(const cplx *matrix, size_type, char, int, int, double, double, int, Vecs, Dephase);
// template void eig::solver::eig<eig::Form::NSYM>(const real *matrix, size_type, Vecs, Dephase);
// template void eig::solver::eig<eig::Form::SYMM>(const cplx *matrix, size_type, Vecs, Dephase);
// template void eig::solver::eig<eig::Form::NSYM>(const cplx *matrix, size_type, Vecs, Dephase);

template<typename Scalar>
void eig::solver::eigs_init(size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Type type, Side side, std::optional<cplx> sigma,
                            Shinv shift_invert, Storage storage, Vecs compute_eigvecs, Dephase remove_phase, Scalar *residual, Lib lib) {
    if(config.loglevel) eig::setLevel(config.loglevel.value());
    eig::log->trace("eigs init");
    result.reset();
    /* clang-format off */
    if(not config.lib) config.lib = lib;

    // Precision settings which are overridden manually
    if(not config.tol) config.tol = 1e-12;
    if(not config.maxIter) config.maxIter = 2000;

    // Other settings that we pass on each invocation of eigs
    if(not config.compute_eigvecs) config.compute_eigvecs = compute_eigvecs;
    if(not config.remove_phase)    config.remove_phase    = remove_phase;
    if(not config.maxNev)          config.maxNev          = nev;
    if(not config.maxNcv)          config.maxNcv          = ncv;
    if(not config.sigma)           config.sigma           = sigma;
    if(not config.shift_invert)    config.shift_invert    = shift_invert;
    if(not config.type)            config.type            = type;
    if(not config.form)            config.form            = form;
    if(not config.ritz)            config.ritz            = ritz;
    if(not config.side)            config.side            = side;
    if(not config.storage)         config.storage         = storage;
    if(not config.compress)        config.compress        = false;

    if(config.initial_guess.empty()) config.initial_guess.push_back({residual,0});
    /* clang-format on */

    if(config.maxNev.value() < 1) config.maxNev = 1;
    if(config.form == Form::NSYM) {
        if(config.maxNev.value() == 1) { config.maxNev = 2; }
    }
    // If ncv < nev we should set a new one as recommended
    if(config.maxNcv.value() <= config.maxNev.value()) config.maxNcv = std::clamp<long>(2l * config.maxNev.value() + 1, config.maxNev.value() + 1, L);
    // Make sure ncv is valid
    config.maxNcv = std::clamp<long>(config.maxNcv.value(), config.maxNev.value() + 1, L);

    if(config.maxNcv.value() > L) throw std::logic_error("ncv > L");
    if(config.maxNev.value() > L) throw std::logic_error("nev > L");
    if(config.maxNcv.value() <= config.maxNev.value()) throw std::logic_error("ncv <= nev");

    if(config.shift_invert == Shinv::ON and not config.sigma) throw std::runtime_error("Sigma must be set to use shift-invert mode");
    config.checkRitz();

    if(not config.logTime) {
        if(eig::log->level() == spdlog::level::trace) config.logTime = 10.0;
        if(eig::log->level() == spdlog::level::debug) config.logTime = 60.0;
        if(eig::log->level() >= spdlog::level::info) config.logTime = 60.0 * 10;
        if(eig::log->level() == spdlog::level::trace) config.logIter = 100;
        if(eig::log->level() == spdlog::level::debug) config.logIter = 1000;
        if(eig::log->level() == spdlog::level::info) config.logIter = 10000;
    }
}
template void eig::solver::eigs_init(size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Type type, Side side, std::optional<cplx> sigma,
                                     Shinv shift_invert, Storage storage, Vecs compute_eigvecs, Dephase remove_phase, real *residual, Lib lib);
template void eig::solver::eigs_init(size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Type type, Side side, std::optional<cplx> sigma,
                                     Shinv shift_invert, Storage storage, Vecs compute_eigvecs, Dephase remove_phase, cplx *residual, Lib lib);

template<typename MatrixProductType>
void eig::solver::set_default_config(const MatrixProductType &matrix) {
    if(config.loglevel) eig::setLevel(config.loglevel.value());
    eig::log->trace("eigs init");
    /* clang-format off */

    if(not config.lib) config.lib = eig::Lib::ARPACK;

    // Precision settings which are overridden manually
    if(not config.tol) config.tol = 1e-12;
    if(not config.maxIter) config.maxIter = 2000;

    // Other settings that we pass on each invocation of eigs
    if(not config.compute_eigvecs) config.compute_eigvecs = eig::Vecs::OFF;
    if(not config.remove_phase)    config.remove_phase    = eig::Dephase::OFF;
    if(not config.maxNev)          config.maxNev          = 1;
    if(not config.maxNcv)          config.maxNcv          = -1;
    if(not config.sigma)           config.sigma           = std::nullopt;
    if(not config.shift_invert)    config.shift_invert    = Shinv::OFF;
    if(not config.type)            config.type            = matrix.get_type();
    if(not config.form)            config.form            = matrix.get_form();
    if(not config.side)            config.side            = matrix.get_side();
    if(not config.ritz)            config.ritz            = eig::Ritz::LM;
    if(not config.storage)         config.storage         = MatrixProductType::storage;
    if(not config.compress)        config.compress        = false;

    /* clang-format on */

    if(config.maxNev.value() < 1) config.maxNev = 1;
    if(config.form == Form::NSYM) {
        if(config.maxNev.value() == 1) { config.maxNev = 2; }
    }
    // If ncv < nev we should set a new one as recommended
    if(config.maxNcv.value() <= config.maxNev.value())
        config.maxNcv = std::clamp<long>(2l * config.maxNev.value() + 1, config.maxNev.value() + 1, matrix.rows());
    // Make sure ncv is valid
    config.maxNcv = std::clamp<long>(config.maxNcv.value(), config.maxNev.value() + 1, matrix.rows());

    if(config.maxNcv.value() > matrix.rows()) throw std::logic_error("ncv > L");
    if(config.maxNev.value() > matrix.rows()) throw std::logic_error("nev > L");
    if(config.maxNcv.value() <= config.maxNev.value()) throw std::logic_error("ncv <= nev");

    if(config.shift_invert == Shinv::ON and not config.sigma) throw std::runtime_error("Sigma must be set to use shift-invert mode");
    config.checkRitz();

    if(not config.logTime) {
        if(eig::log->level() == spdlog::level::trace) config.logTime = 10.0;
        if(eig::log->level() == spdlog::level::debug) config.logTime = 60.0;
        if(eig::log->level() >= spdlog::level::info) config.logTime = 60.0 * 10;
        if(eig::log->level() == spdlog::level::trace) config.logIter = 100;
        if(eig::log->level() == spdlog::level::debug) config.logIter = 1000;
        if(eig::log->level() == spdlog::level::info) config.logIter = 10000;
    }
}

template void eig::solver::set_default_config(const MatVecDense<eig::real> &matrix);
template void eig::solver::set_default_config(const MatVecDense<eig::cplx> &matrix);
template void eig::solver::set_default_config(const MatVecSparse<eig::real> &matrix);
template void eig::solver::set_default_config(const MatVecSparse<eig::cplx> &matrix);
template void eig::solver::set_default_config(const MatVecMPO<eig::real> &matrix);
template void eig::solver::set_default_config(const MatVecMPO<eig::cplx> &matrix);

template<typename Scalar, eig::Storage storage>
void eig::solver::eigs(const Scalar *matrix, size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma,
                       Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, Scalar *residual) {
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

template void eig::solver::eigs(const eig::real *matrix, size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::real *residual);

template void eig::solver::eigs(const eig::cplx *matrix, size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::cplx *residual);

template<typename MatrixProductType>
void eig::solver::eigs(MatrixProductType &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
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

template void eig::solver::eigs(MatVecDense<eig::real> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<eig::cplx> sigma,
                                Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::real *residual);

template void eig::solver::eigs(MatVecDense<eig::cplx> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<eig::cplx> sigma,
                                Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::cplx *residual);

template void eig::solver::eigs(MatVecSparse<eig::real> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<eig::cplx> sigma,
                                Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::real *residual);

template void eig::solver::eigs(MatVecSparse<eig::cplx> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<eig::cplx> sigma,
                                Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::cplx *residual);

template void eig::solver::eigs(MatVecMPO<eig::real> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<eig::cplx> sigma,
                                Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::real *residual);

template void eig::solver::eigs(MatVecMPO<eig::cplx> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<eig::cplx> sigma,
                                Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::cplx *residual);

template<typename MatrixProductType>
void eig::solver::eigs(MatrixProductType &matrix) {
    auto t_eig = tid::tic_scope("eig");
    result.reset();
    set_default_config(matrix);
    switch(config.lib.value()) {
        case Lib::ARPACK: {
            if(config.tag.empty()) config.tag = "arpack";
            solver_arpack<MatrixProductType> solver(matrix, config, result);
            solver.eigs();
            break;
        }
        case Lib::PRIMME: {
            if(config.tag.empty()) config.tag = "primme";
            eigs_primme(matrix);
            break;
        }
    }
}

template void eig::solver::eigs(MatVecDense<eig::real> &matrix);
template void eig::solver::eigs(MatVecDense<eig::cplx> &matrix);
template void eig::solver::eigs(MatVecSparse<eig::real> &matrix);
template void eig::solver::eigs(MatVecSparse<eig::cplx> &matrix);
template void eig::solver::eigs(MatVecMPO<eig::real> &matrix);
template void eig::solver::eigs(MatVecMPO<eig::cplx> &matrix);