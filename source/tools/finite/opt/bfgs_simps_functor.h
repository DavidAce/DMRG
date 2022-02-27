#pragma once

#include <ceres/first_order_function.h>
#include <unsupported/Eigen/CXX11/Tensor>
namespace tid {
    class ur;
}

namespace spdlog {
    class logger;
}
namespace tools::finite::opt::internal {
    template<typename Scalar>
    class bfgs_simps_functor : public ceres::FirstOrderFunction {
        private:
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        using TensorMap3 = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>;
        using TensorMap4 = Eigen::TensorMap<const Eigen::Tensor<Scalar, 4>>;

        int                    num_parameters; // Includes lagrangian multiplier(s)
        long                   size;
        Eigen::DSizes<long, 3> dims;

        mutable Eigen::Tensor<Scalar, 3> Hm_tensor, Hv_tensor;
        const TensorMap3                &mps, &envL, envR;
        const TensorMap4                &mpo;

        void compute_Hv(const VectorType &v) const; // Computes  H|φ>
        void compute_Hm(const VectorType &m) const; // Computes  H|φ>

        public:
        explicit bfgs_simps_functor(const TensorMap3 &mps_,  /*!< The current mps  */
                                    const TensorMap3 &envL_, /*!< The left block tensor for H.  */
                                    const TensorMap3 &envR_, /*!< The right block tensor for H.  */
                                    const TensorMap4 &mpo_   /*!< The H MPO's  */
        );
        int  NumParameters() const final;
        bool Evaluate(const double *v_double_double, double *fx, double *grad_double_double) const final;
    };

}
