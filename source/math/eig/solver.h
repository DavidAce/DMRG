#pragma once
#include "enums.h"
#include "settings.h"
#include "solution.h"

struct primme_params;

namespace eig {

    class solver {
        public:
        eig::settings config;
        eig::solution result;

        solver();
        explicit solver(size_t loglevel);
        void setLogLevel(size_t loglevel);
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
        template<typename Scalar>
        void eigs_init(size_type L, size_type nev, size_type ncv, Ritz ritz = Ritz::LM, Form form = Form::SYMM, Type type = Type::REAL, Side side = Side::R,
                       std::optional<cplx> sigma = std::nullopt, Shinv shinv = Shinv::OFF, Storage storage = Storage::DENSE, Vecs compute_eigvecs_ = Vecs::OFF,
                       Dephase remove_phase_ = Dephase::OFF, Scalar* residual = nullptr, Lib lib = Lib::ARPACK);

        template<typename MatrixProductType>
        void set_default_config(const MatrixProductType & matrix);


        template<typename Scalar, Storage storage = Storage::DENSE>
        void eigs(const Scalar *matrix, size_type L, size_type nev, size_type ncv, Ritz ritz = Ritz::SR, Form form = Form::SYMM, Side side = Side::R,
                  std::optional<cplx> sigma = std::nullopt, Shinv shinv = Shinv::OFF, Vecs vecs = Vecs::ON, Dephase remove_phase = Dephase::OFF,
                  Scalar *residual = nullptr);

        template<typename MatrixProductType>
        void eigs(MatrixProductType &matrix, size_type nev, size_type ncv, Ritz ritz = Ritz::SR, Form form = Form::SYMM, Side side = Side::R,
                  std::optional<cplx> sigma = std::nullopt, Shinv shinv = Shinv::OFF, Vecs vecs = Vecs::ON, Dephase remove_phase = Dephase::OFF,
                  typename MatrixProductType::Scalar *residual = nullptr);

        template<typename MatrixProductType>
        void eigs(MatrixProductType &matrix);


        template<typename MatrixProductType>
        int eigs_primme(MatrixProductType &matrix);


        private:
        template<typename MatrixProductType>
        static void MultAx_wrapper(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr);
        template<typename MatrixProductType>
        static void MultOPv_wrapper(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr);

        template<typename MatrixProductType>
        static void GradientConvTest(double *eval, void *evec, double *rNorm, int *isconv,
            struct primme_params *primme, int *ierr);

    };
}
