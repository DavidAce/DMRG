#pragma once

// Testing cppoptlib
#include "cppoptlib/function.h"
#include "cppoptlib/solver/lbfgs.h"

class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
class AlgorithmStatus;


using FunctionXd = cppoptlib::function::Function<double>;

class Function : public FunctionXd {
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using FunctionXd::hessian_t;
    using FunctionXd::vector_t;

    using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorType = Eigen::Matrix<double, Eigen::Dynamic, 1>;
    mutable double max_grad_norm = 0;
    mutable double norm;
    mutable double var;
    mutable double nHn, nH2n;
    mutable double vv, log10var;
    mutable Eigen::Tensor<double, 3> Hn_tensor, H2n_tensor;
    long                             size; //Does not include lagrange multiplier (i.e. num_paramters - 1 if LagrangeNorm::ON)
    Eigen::DSizes<long, 3>           dims;
    Eigen::Tensor<double, 3>         envL, envR;
    Eigen::Tensor<double, 3>         env2L, env2R;
    Eigen::Tensor<double, 4>         mpo, mpo2;
    mutable bool                     print_path    = true;
    mutable bool                     readyCompress = false;
    void                             get_Hn(const VectorType &v) const;
    void                             get_H2n(const VectorType &v) const;

    Function(const TensorsFinite &tensors, const AlgorithmStatus &status);

    scalar_t operator()(const vector_t &v) const;

    void Gradient(const vector_t &x, vector_t *grad) const;
};