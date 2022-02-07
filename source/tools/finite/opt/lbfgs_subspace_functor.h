#pragma once
#include "lbfgs_base_functor.h"
#include <math/tenx/fwd_decl.h>

namespace tools::finite::opt::internal {
    template<typename Scalar>
    class lbfgs_subspace_functor : public lbfgs_base_functor {
        protected:
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
        const MatrixType      &H2;
        const Eigen::VectorXd &eigvals;

        public:
        template<typename>
        friend class NormParametrization;
        explicit lbfgs_subspace_functor(const TensorsFinite &tensors, const AlgorithmStatus &status, const MatrixType &H2_subspace_,
                                        const Eigen::VectorXd &eigvals_);
        bool Evaluate(const double *v_double_double, double *fx, double *grad_double_double) const final;
    };
}
