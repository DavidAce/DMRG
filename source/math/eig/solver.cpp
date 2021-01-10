//
// Created by david on 2020-06-04.
//

#include "solver.h"
#include "arpack_solver/arpack_solver.h"
#include "arpack_solver/matrix_product_dense.h"
#include "arpack_solver/matrix_product_hamiltonian.h"
#include "arpack_solver/matrix_product_sparse.h"

eig::solver::solver() {
    if(not eig::log) eig::log = tools::Logger::setLogger("eig", 2, true);
}

eig::solver::solver(size_t logLevel) : solver() { tools::Logger::setLogLevel(eig::log, logLevel); }
eig::solver::solver(std::shared_ptr<spdlog::logger> logger) { eig::log = std::move(logger); }

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
    int info = 0;
    try {
        if constexpr(std::is_same_v<Scalar, real>) {
            eig_init(form, Type::REAL, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) info = dsyevd(matrix, L);
            if constexpr(form == Form::NSYM) info = dgeev(matrix, L);
        } else if constexpr(std::is_same_v<Scalar, cplx>) {
            eig_init(form, Type::CPLX, compute_eigvecs_, remove_phase_);
            if constexpr(form == Form::SYMM) info = zheevd(matrix, L);
            if constexpr(form == Form::NSYM) info = zgeev(matrix, L);
        } else {
            throw std::runtime_error("Unknown type");
        }

    } catch(std::exception &ex) {
        eig::log->error("Eigenvalue solver failed: {}", ex.what());
        throw std::runtime_error(fmt::format("Eigenvalue solver Failed: {}", ex.what()));
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

void eig::solver::eigs_init(size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Type type, Side side, std::optional<cplx> sigma,
                            Shinv shift_invert, Storage storage, Vecs compute_eigvecs, Dephase remove_phase) {
    eig::log->trace("eigs init");
    result.reset();
    /* clang-format off */
    // Precision settings which are overridden manually
    if(not config.eigThreshold) config.eigThreshold = 1e-12;
    if(not config.eigMaxIter) config.eigMaxIter = 2000;

    // Other settings that we pass on each invocation of eigs
    if(not config.compute_eigvecs) config.compute_eigvecs = compute_eigvecs;
    if(not config.remove_phase)    config.remove_phase    = remove_phase;
    if(not config.eigMaxNev)       config.eigMaxNev       = nev;
    if(not config.eigMaxNcv)       config.eigMaxNcv       = ncv;
    if(not config.sigma)           config.sigma           = sigma;
    if(not config.shift_invert)    config.shift_invert    = shift_invert;
    if(not config.type)            config.type            = type;
    if(not config.form)            config.form            = form;
    if(not config.ritz)            config.ritz            = ritz;
    if(not config.side)            config.side            = side;
    if(not config.storage)         config.storage         = storage;
    /* clang-format on */

    if(config.eigMaxNev.value() < 1) config.eigMaxNev = 1;
    if(config.eigMaxNcv.value() <= 1) config.eigMaxNcv = 16;

    if(config.eigMaxNcv.value() < config.eigMaxNev.value()) {
        if(config.eigMaxNev.value() >= 1 and config.eigMaxNev.value() <= 16) {
            config.eigMaxNcv = 8 + static_cast<int>(std::ceil(1.5 * static_cast<double>(nev)));
        } else if(nev > 16 and nev <= L) {
            config.eigMaxNcv = 2 * nev;
        }
    }

    if(config.form == Form::NSYM) {
        if(config.eigMaxNev.value() == 1) { config.eigMaxNev = 2; }
    }
    if(config.eigMaxNcv.value() >= L) { config.eigMaxNcv = (L + config.eigMaxNev.value()) / 2; }

    assert(config.eigMaxNcv.value() <= L and "Ncv > L");
    assert(config.eigMaxNcv.value() >= config.eigMaxNev.value() and "Ncv < Nev");
    assert(config.eigMaxNev.value() <= L and "Nev > L");

    if(config.shift_invert == Shinv::ON and not config.sigma) throw std::runtime_error("Sigma must be set to use shift-invert mode");
    config.checkRitz();
}

template<typename Scalar, eig::Storage storage>
void eig::solver::eigs(const Scalar *matrix, size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma,
                       Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, Scalar *residual) {
    bool is_cplx = std::is_same<std::complex<double>, Scalar>::value;
    Type type    = is_cplx ? Type::CPLX : Type::REAL;
    eigs_init(L, nev, ncv, ritz, form, type, side, sigma, shift_invert, storage, compute_eigvecs, remove_phase);
    if constexpr(storage == Storage::DENSE) {
        auto                                      matrix_dense = MatrixProductDense<Scalar>(matrix, L, true);
        arpack_solver<MatrixProductDense<Scalar>> solver(matrix_dense, config, result, residual);
        solver.eigs();
    } else if constexpr(storage == Storage::SPARSE) {
        auto                                              matrix_sparse = MatrixProductSparse<Scalar, false>(matrix, L, true);
        arpack_solver<MatrixProductSparse<Scalar, false>> solver(matrix_sparse, config, result, residual);
        solver.eigs();
    }
}

template void eig::solver::eigs(const eig::real *matrix, size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::real *residual);

template void eig::solver::eigs(const eig::cplx *matrix, size_type L, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::cplx *residual);

template<typename MatrixProductType>
void eig::solver::eigs(MatrixProductType &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side, std::optional<cplx> sigma, Shinv shift_invert,
                       Vecs compute_eigvecs, Dephase remove_phase, typename MatrixProductType::Scalar *residual) {
    using Scalar = typename MatrixProductType::Scalar;
    Type type;
    if constexpr(std::is_same_v<Scalar, real>)
        type = Type::REAL;
    else if constexpr(std::is_same_v<Scalar, cplx>)
        type = Type::CPLX;
    else
        throw std::runtime_error("Unsupported type");
    eigs_init(matrix.rows(), nev, ncv, ritz, form, type, side, sigma, shift_invert, matrix.storage, compute_eigvecs, remove_phase);
    arpack_solver<MatrixProductType> solver(matrix, config, result, residual);
    solver.eigs();
}

template void eig::solver::eigs(MatrixProductDense<eig::real> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::real *residual);

template void eig::solver::eigs(MatrixProductDense<eig::cplx> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::cplx *residual);

template void eig::solver::eigs(MatrixProductSparse<eig::real> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::real *residual);

template void eig::solver::eigs(MatrixProductSparse<eig::cplx> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::cplx *residual);

template void eig::solver::eigs(MatrixProductHamiltonian<eig::real> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::real *residual);

template void eig::solver::eigs(MatrixProductHamiltonian<eig::cplx> &matrix, size_type nev, size_type ncv, Ritz ritz, Form form, Side side,
                                std::optional<eig::cplx> sigma, Shinv shift_invert, Vecs compute_eigvecs, Dephase remove_phase, eig::cplx *residual);
