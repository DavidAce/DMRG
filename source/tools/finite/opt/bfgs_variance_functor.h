#pragma once

#include "bfgs_base_functor.h"
#include <unsupported/Eigen/CXX11/Tensor>

namespace tools::finite::opt::internal {
    template<typename Scalar, LagrangeNorm lagrangeNorm>
    class bfgs_variance_functor : public bfgs_base_functor {
        private:
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        mutable VectorType               residual;
        mutable Eigen::Tensor<Scalar, 3> Hn_tensor, H2n_tensor, H2r_tensor, H2H2n_tensor;
        Eigen::Tensor<Scalar, 3>         envL, envR;
        Eigen::Tensor<Scalar, 3>         env2L, env2R;
        Eigen::Tensor<Scalar, 4>         mpo, mpo2;
        mutable Scalar                   shift;
        mutable bool                     print_path    = true;
        mutable bool                     readyCompress = false;
        mutable bool                     readyShift    = false;
        void                             get_Hn(const VectorType &n) const;
        void                             get_H2n(const VectorType &n) const;
        void                             get_H2r(const VectorType &r) const;
        void                             get_H2H2n(const VectorType &n) const;

        public:
        template<typename>
        friend class NormParametrization;
        explicit bfgs_variance_functor(const TensorsFinite &tensors, const AlgorithmStatus &status);
        bool Evaluate(const double *v_double_double, double *fx, double *grad_double_double) const final;
        void compress();
        void set_shift(double shift_);
    };

}
