#include "cppoptlib_direct_functor.h"
#include "opt-internal.h"
#include <algorithms/AlgorithmStatus.h>
#include <math/num.h>
#include <math/svd.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/state/StateFinite.h>
#include <tensors/TensorsFinite.h>
#include <tools/common/contraction.h>

Function::Function(const TensorsFinite &tensors, const AlgorithmStatus &status) {
    mpo                 = tensors.get_multisite_mpo().real();
    mpo2                = tensors.get_multisite_mpo_squared().real();
    const auto &env_ene = tensors.get_multisite_env_ene_blk();
    const auto &env_var = tensors.get_multisite_env_var_blk();
    envL                = env_ene.L.real();
    envR                = env_ene.R.real();
    env2L               = env_var.L.real();
    env2R               = env_var.R.real();
    dims                = tensors.active_problem_dims();
    Hn_tensor.resize(dims);
    H2n_tensor.resize(dims);
    size = dims[0] * dims[1] * dims[2];
}

Function::scalar_t Function::operator()(const vector_t &v) const {
    norm = v.norm();
    get_H2n(v);
    get_Hn(v);
    auto Hn  = Eigen::Map<VectorType>(Hn_tensor.data(), Hn_tensor.size());
    auto H2n = Eigen::Map<VectorType>(H2n_tensor.data(), H2n_tensor.size());
    nHn      = v.dot(Hn);
    nH2n     = v.dot(H2n);
    var      = nH2n - nHn * nHn;
    log10var = std::log10(std::abs(var));
    return log10var;
}

void Function::Gradient(const vector_t &v, vector_t *grad) const {
    double pref = 2.0; // Prefactor
    auto Hn  = Eigen::Map<VectorType>(Hn_tensor.data(), Hn_tensor.size());
    auto H2n = Eigen::Map<VectorType>(H2n_tensor.data(), H2n_tensor.size());
    if(grad != nullptr) {
        auto one_over_norm = 1.0 / norm;
        auto var_1         = 1.0 / std::real(var) / std::log(10);
        *grad              = var_1 * pref * one_over_norm * (H2n - 2.0 * nHn * Hn - (nH2n - 2.0 * nHn * nHn) * v);
        max_grad_norm      = (*grad).template lpNorm<Eigen::Infinity>();
    }

    // Here we define the norm constraint by using the a lagrange multiplier trick:
    //      f(x)
    // is replaced with
    //      L(x,lambda) = f(x) + lambda * g(x)
    // where g(x) = | <x|x> - 1 |

    //    double lambda     = 1.0;                // a dummy
    //    double constraint = std::abs(vv - 1.0); // aka g(x)
    //    if(fx != nullptr) fx[0] += lambda * constraint;
    //    if(grad_double_double != nullptr) {
    //        Eigen::Map<VectorType> grad_w_multiplier(reinterpret_cast<Scalar *>(grad_double_double), size + 1);
    //        grad_w_multiplier.topRows(size) += lambda * pref * num::sign(vv - 1.0) * v; // aka  += lambda * dg(x)/dx = lambda * sign(x) * x
    //        grad_w_multiplier.bottomRows(1)[0] = constraint;
    //        //            max_grad_norm           = grad_w_multiplier.template lpNorm<Eigen::Infinity>();
    //    }
}

void Function::get_H2n(const VectorType &v) const {
    auto v_tensor = Eigen::TensorMap<const Eigen::Tensor<const double, 3>>(v.derived().data(), dims);
    tools::common::contraction::matrix_vector_product(H2n_tensor, v_tensor, mpo2, env2L, env2R);
}

void Function::get_Hn(const VectorType &v) const {
    auto v_tensor = Eigen::TensorMap<const Eigen::Tensor<const double, 3>>(v.derived().data(), dims);
    tools::common::contraction::matrix_vector_product(Hn_tensor, v_tensor, mpo, envL, envR);
}

