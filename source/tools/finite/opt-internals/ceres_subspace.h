//
// Created by david on 2019-07-15.
//

#ifndef DMRG_CERES_SUBSPACE_H
#define DMRG_CERES_SUBSPACE_H

#include <tools/finite/opt.h>

namespace tools::finite::opt{
    namespace internals{
        template<typename Scalar>
        class ceres_subspace_functor : public ceres_base_functor{
        private:
            using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
            using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
            const Eigen::MatrixXcd &eigvecs;
            const Eigen::VectorXd  &eigvals;
            MatrixType             H2;
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            explicit ceres_subspace_functor(const class_finite_state & state,
                                            const class_simulation_status & sim_status,
                                            const Eigen::MatrixXcd & eigvecs_,
                                            const Eigen::VectorXd  & eigvals_);
            bool Evaluate(const double* v_double_double,
                          double* fx,
                          double* grad_double_double) const override;
        };
    }
}

#endif //DMRG_CERES_SUBSPACE_H
