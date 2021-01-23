//
// Created by david on 2019-07-15.
//

#pragma once
#include <general/eigen_tensor_fwd_decl.h>
#include <tools/finite/opt-internal/ceres_base.h>

namespace tools::finite::opt::internal {
    template<typename Scalar>
    class ceres_subspace_functor : public ceres_base_functor {
        private:
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        const MatrixType &     H2;
        const Eigen::VectorXd &eigvals;

        public:
        explicit ceres_subspace_functor(const class_tensors_finite &tensors, const class_algorithm_status &status, const MatrixType &H2_subspace_,
                                        const Eigen::VectorXd &eigvals_);
        bool Evaluate(const double *v_double_double, double *fx, double *grad_double_double) const final;
    };
}
