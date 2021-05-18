#pragma once
#include "enums.h"
#include <typeindex>
#include <vector>

namespace eig {
    class solver;
    template<typename MatrixType>
    class arpack_solver;

    class solution {
        public:
        friend class solver;
        template<typename MatrixType>
        friend class arpack_solver;

        private:
        std::vector<real> eigvals_real;
        std::vector<real> eigvals_imag;
        std::vector<cplx> eigvals_cplx;
        std::vector<real> eigvecsR_real;
        std::vector<real> eigvecsR_imag;
        std::vector<real> eigvecsL_real;
        std::vector<real> eigvecsL_imag;
        std::vector<cplx> eigvecsR_cplx;
        std::vector<cplx> eigvecsL_cplx;

        void build_eigvecs_cplx();
        void build_eigvecs_real();
        void build_eigvals_cplx();
        void build_eigvals_real();

        struct Meta {
            eig::size_type rows           = 0;
            eig::size_type cols           = 0;
            eig::size_type iter           = 0;
            eig::size_type nev            = 0; // Requested eigenvectors. aka cols
            eig::size_type nev_converged  = 0; // Converged eigenvectors
            eig::size_type n              = 0; // Linear dimension of the input matrix to diagonalize, aka rows.
            eig::size_type ncv            = 0;
            double         tol            = 0;
            int            counter        = 0;
            double         time_total     = 0;
            double         time_matprod   = 0;
            double         time_prep      = 0;
            bool           eigvals_found  = false;
            bool           eigvecsR_found = false;
            bool           eigvecsL_found = false;
            bool           arnoldi_found  = false;
            Form           form;
            Type           type;
            //            Side  side;
            std::string          ritz;
            std::complex<double> sigma = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        };

        public:
        Meta meta;

        template<typename Scalar, Side side = Side::R>
        auto &get_eigvecs() {
            if constexpr(std::is_same_v<Scalar, real>) {
                build_eigvecs_real();
                if constexpr(side == Side::R) return eigvecsR_real;
                if constexpr(side == Side::L) return eigvecsL_real;
                if constexpr(side == Side::LR) return std::pair(eigvecsL_real, eigvecsR_real);
            }
            if constexpr(std::is_same_v<Scalar, cplx>) {
                build_eigvecs_cplx();
                if constexpr(side == Side::R) return eigvecsR_cplx;
                if constexpr(side == Side::L) return eigvecsL_cplx;
                if constexpr(side == Side::LR) return std::pair(eigvecsL_cplx, eigvecsR_cplx);
            }
        }

        template<typename Scalar>
        auto &get_eigvecs(Side side = Side::R) {
            if(side == Side::R) return get_eigvecs<Scalar, Side::R>();
            if(side == Side::L) return get_eigvecs<Scalar, Side::L>();
            throw std::runtime_error("Cannot return both L and R eigenvectors");
        }

        template<Form form, Type type = Type::CPLX, Side side = Side::R>
        auto &get_eigvecs() {
            if constexpr(type == Type::REAL) {
                if constexpr(form == Form::SYMM) return get_eigvecs<real, side>();
                if constexpr(form == Form::NSYM) return get_eigvecs<cplx, side>();
            }
            if constexpr(type == Type::CPLX) {
                if constexpr(form == Form::SYMM) return get_eigvecs<cplx, side>();
                if constexpr(form == Form::NSYM) return get_eigvecs<cplx, side>();
            }
        }

        template<Type type, Form form, Side side = Side::R>
        auto &get_eigvecs() {
            return get_eigvecs<form, type, side>();
        }

        template<typename Scalar, Form form, Side side = Side::R>
        auto &get_eigvecs() {
            if constexpr(std::is_same<real, Scalar>::value) return get_eigvecs<form, Type::REAL, side>();
            if constexpr(std::is_same<cplx, Scalar>::value) return get_eigvecs<form, Type::CPLX, side>();
        }

        template<typename Scalar>
        auto &get_eigvals() {
            if constexpr(std::is_same_v<Scalar, real>) {
                build_eigvals_real();
                return eigvals_real;
            }
            if constexpr(std::is_same_v<Scalar, cplx>) {
                build_eigvals_cplx();
                return eigvals_cplx;
            }
        }

        template<Form form = Form::SYMM>
        auto &get_eigvals() {
            if constexpr(form == Form::SYMM) { return get_eigvals<real>(); }
            if constexpr(form == Form::NSYM) { return get_eigvals<cplx>(); }
        }

        void reset() {
            eigvals_real.clear();
            eigvals_imag.clear();
            eigvals_cplx.clear();
            eigvecsR_real.clear();
            eigvecsR_imag.clear();
            eigvecsL_real.clear();
            eigvecsL_imag.clear();
            eigvecsR_cplx.clear();
            eigvecsL_cplx.clear();
            meta = Meta();
        }

        bool eigvecs_are_real() { return meta.form == Form::SYMM and meta.type == Type::REAL; }

        bool eigvals_are_real() { return meta.form == Form::SYMM; }

        std::type_index get_eigvecs_type() {
            if(eigvecs_are_real())
                return typeid(real);
            else
                return typeid(cplx);
        }
    };
}
