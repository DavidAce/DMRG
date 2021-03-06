#pragma once
#include "enums.h"
#include "settings.h"
#include "solution.h"
namespace eig {

    class solver {
        public:
        eig::settings config;
        eig::solution result;

        solver();
        explicit solver(size_t logLevel);
        explicit solver(std::shared_ptr<spdlog::logger> logger);
        void setLogLevel(spdlog::level::level_enum level) const;
        void setLogLevel(size_t level) const;
        template<typename Scalar>
        void subtract_phase(std::vector<Scalar> &eigvecs, size_type L, size_type nev);

        // Functions for full diagonalization of explicit matrix
        int dsyevd(const real *matrix, size_type L);
        int zheevd(const cplx *matrix, size_type L);
        int dgeev(const real *matrix, size_type L);
        int zgeev(const cplx *matrix, size_type L);

        void eig_init(Form form, Type type, Vecs compute_eigvecs, Dephase remove_phase_);
        template<Form form = Form::SYMM, typename Scalar>
        void eig(const Scalar *matrix, size_type L, Vecs compute_eigvecs = Vecs::ON, Dephase remove_phase_ = Dephase::OFF);

        // Functions for few eigensolutions

        void eigs_init(size_type L, size_type nev, size_type ncv, Ritz ritz = Ritz::LM, Form form = Form::SYMM, Type type = Type::REAL, Side side = Side::R,
                       std::optional<cplx> sigma = std::nullopt, Shinv shinv = Shinv::OFF, Storage storage = Storage::DENSE, Vecs compute_eigvecs_ = Vecs::OFF,
                       Dephase remove_phase_ = Dephase::OFF, Lib lib = Lib::ARPACK);

        template<typename Scalar, Storage storage = Storage::DENSE>
        void eigs(const Scalar *matrix, size_type L, size_type nev, size_type ncv, Ritz ritz = Ritz::SR, Form form = Form::SYMM, Side side = Side::R,
                  std::optional<cplx> sigma = std::nullopt, Shinv shinv = Shinv::OFF, Vecs vecs = Vecs::ON, Dephase remove_phase = Dephase::OFF,
                  Scalar *residual = nullptr);

        template<typename MatrixProductType>
        void eigs(MatrixProductType &matrix, size_type nev, size_type ncv, Ritz ritz = Ritz::SR, Form form = Form::SYMM, Side side = Side::R,
                  std::optional<cplx> sigma = std::nullopt, Shinv shinv = Shinv::OFF, Vecs vecs = Vecs::ON, Dephase remove_phase = Dephase::OFF,
                  typename MatrixProductType::Scalar *residual = nullptr);

        template<typename MatrixProductType>
        int eigs_primme(MatrixProductType &matrix, typename MatrixProductType::Scalar *residual = nullptr);
    };
}
