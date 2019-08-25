//
// Created by david on 2019-07-15.
//

#ifndef DMRG_CERES_DIRECT_FUNCTOR_H
#define DMRG_CERES_DIRECT_FUNCTOR_H

#include <tools/finite/opt.h>

namespace tools::finite::opt{
    namespace internals{


        template<typename Scalar>
        class ceres_direct_functor : public ceres_base_functor{
        private:
            using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
            using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
            Eigen::Tensor<Scalar,3> envL, envR;
            Eigen::Tensor<Scalar,4> env2L, env2R;
            Eigen::Tensor<Scalar,4> mpo;
            Eigen::DSizes<long,3>   dsizes;
            mutable bool print_path = true;
            mutable Eigen::Tensor<Scalar,3> Hv_tensor, H2v_tensor;
            void get_Hv   (const VectorType &v) const;
            void get_H2v  (const VectorType &v) const;
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            explicit ceres_direct_functor(const class_finite_state & state, const class_simulation_status &sim_status);
            bool Evaluate(const double* v_double_double,
                          double* fx,
                          double* grad_double_double) const override;
        };

    }

}




#endif //DMRG_CERES_DIRECT_FUNCTOR_H
