//
// Created by david on 2019-12-17.
//
#include "ceres_pthread.h"

//template<typename Scalar>
//opt::Rosenbrock<Scalar>::Rosenbrock(MatrixType &H_) : RosenbrockBase(H_) {}
template<typename Scalar>
opt::Rosenbrock<Scalar>::Rosenbrock(const MatrixType &H_, const MatrixType &H2_) : RosenbrockBase(H_,H2_) {}

template<typename Scalar>
bool opt::Rosenbrock<Scalar>::Evaluate
            (
              const double *v_ptr,
              double *fx,
              double *grad_ptr) const
      {
        Eigen::Map<const VectorType> v(v_ptr, NumParameters());
        double vv = v.squaredNorm();
        double norm = std::sqrt(vv);
        VectorType Hv = H * v;
        VectorType H2v = H2 * v;
        double ene = v.adjoint() * Hv;
        double ene2 = v.adjoint() * H2v;

        double var = std::abs(ene2 - ene);
        variance = var;
        energy = ene;
        double norm_offset = std::abs(1 - norm);
        double log10var = std::log10(var);
        if (fx != nullptr) {
        fx[0] = log10var + norm_offset;
        }

        if (grad_ptr != nullptr) {
        Eigen::Map<VectorType> grad(grad_ptr, NumParameters());
        auto vv_1 = std::pow(vv, -1);
        auto var_1 = 1.0 / var / std::log(10);
        grad = var_1 * vv_1 * 2.0 * (H2v - 2.0 * ene * Hv - (ene2 - 2.0 * ene * ene) * v);
        grad += 2.0 * norm_offset * v;
        }
        return true;
}

template class opt::Rosenbrock<double>;
//template class opt::Rosenbrock<std::complex<double>>;


