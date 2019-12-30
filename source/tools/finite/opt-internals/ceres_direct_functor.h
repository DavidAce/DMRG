//
// Created by david on 2019-07-15.
//

#pragma once

#include <tools/finite/opt.h>

namespace tools::finite::opt::internal{
        template<typename Scalar>
        class ceres_direct_functor : public ceres_base_functor{
        private:
            using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
            using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
            mutable Eigen::Tensor<Scalar,3> Hv_tensor, H2v_tensor;
            Eigen::DSizes<long,3>   dsizes;
            Eigen::Tensor<Scalar,3> envL, envR;
            Eigen::Tensor<Scalar,4> env2L, env2R;
            Eigen::Tensor<Scalar,4> mpo;
            mutable bool print_path = true;
            void get_Hv   (const VectorType &v) const;
            void get_H2v  (const VectorType &v) const;
        public:
            explicit ceres_direct_functor(const class_state_finite & state, const class_simulation_status &sim_status);
            bool Evaluate(const double* v_double_double,
                          double* fx,
                          double* grad_double_double) const final;
        };

    }



