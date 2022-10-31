#pragma once
#include "enums.h"
#include <typeindex>
#include <vector>

namespace eig {
    class solver;
    template<typename MatrixType>
    class solver_arpack;

    class solution {
        public:
        friend class solver;
        template<typename MatrixType>
        friend class solver_arpack;

        private:
        mutable std::vector<real> eigvals_real;
        mutable std::vector<real> eigvals_imag;
        mutable std::vector<cplx> eigvals_cplx;
        mutable std::vector<real> eigvecsR_real;
        mutable std::vector<real> eigvecsR_imag;
        mutable std::vector<real> eigvecsL_real;
        mutable std::vector<real> eigvecsL_imag;
        mutable std::vector<cplx> eigvecsR_cplx;
        mutable std::vector<cplx> eigvecsL_cplx;

        void build_eigvecs_cplx() const;
        void build_eigvecs_real() const;
        void build_eigvals_cplx() const;
        void build_eigvals_real() const;

        struct Meta {
            eig::size_type rows           = 0;
            eig::size_type cols           = 0;
            eig::size_type iter           = 0; // Number of iterations reported by solver
            eig::size_type nev            = 0; // Requested eigenvectors. aka cols
            eig::size_type nev_converged  = 0; // Converged eigenvectors
            eig::size_type n              = 0; // Linear dimension of the input matrix to diagonalize, aka rows.
            eig::size_type ncv            = 0;
            double         tol            = 0;
            long           num_mv         = 0; // Number of matrix-vector multiplications counted
            long           num_pc         = 0; // Number of preconditioner operations y = A⁻¹x counted
            long           num_op         = 0; // Number of operations, for instance in primme if maxblocksize > 1
            double         time_total     = 0;
            double         time_mv        = 0;
            double         time_op        = 0;
            double         time_pc        = 0;
            double         time_prep      = 0;
            bool           eigvals_found  = false;
            bool           eigvecsR_found = false;
            bool           eigvecsL_found = false;
            bool           arnoldi_found  = false;

            std::vector<double> residual_norms = {}; /*!< For each eigenpair, stores the norm of its residual_norm vector (Ax - λx) */

            double last_res_norm  = std::numeric_limits<double>::quiet_NaN();
            double last_grad_max  = std::numeric_limits<double>::quiet_NaN();
            long   last_grad_iter = 0;
            double last_grad_time = 0;
            long   last_log_iter  = 0;
            double last_log_time  = 0;

            Form                 form;
            Type                 type;
            std::string          tag;
            std::string          ritz;
            std::complex<double> sigma = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        };

        public:
        Meta meta;

        template<typename Scalar, Side side = Side::R>
        std::vector<Scalar> &get_eigvecs() const;

        template<typename Scalar>
        std::vector<Scalar> &get_eigvecs(Side side = Side::R) const;

        template<Form form, Type type = Type::CPLX, Side side = Side::R>
        auto &get_eigvecs() const {
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
        auto &get_eigvecs() const {
            return get_eigvecs<form, type, side>();
        }

        template<typename Scalar, Form form, Side side = Side::R>
        std::vector<Scalar> &get_eigvecs() const;

        template<typename Scalar>
        std::vector<Scalar> &get_eigvals() const;

        template<Form form = Form::SYMM>
        auto &get_eigvals() const {
            if constexpr(form == Form::SYMM) { return get_eigvals<real>(); }
            if constexpr(form == Form::NSYM) { return get_eigvals<cplx>(); }
        }

        const std::vector<double> &get_resnorms() const;

        void reset();

        bool eigvecs_are_real() const;

        bool eigvals_are_real() const;

        std::type_index get_eigvecs_type() const;
    };
}
