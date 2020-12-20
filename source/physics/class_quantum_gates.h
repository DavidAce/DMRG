#pragma once
#include <Eigen/Core>
#include <complex>
#include <general/eigen_tensor_fwd_decl.h>
#include <vector>

namespace qm {
    /* clang-format off */
    using Scalar = std::complex<double>;
    using cplx = std::complex<double>;
    using real = double;

    struct Gate {
        Eigen::Tensor<Scalar,2> op;
        std::vector<size_t> pos;
        Gate(const Eigen::Tensor<Scalar,2>  & op_, const std::vector<size_t> & pos_);
        Gate(const Eigen::Tensor<Scalar,2>  & op_, const std::vector<size_t> & pos_, Scalar delta);
        void exp_inplace(Scalar delta);
        [[nodiscard]] Gate exp(Scalar delta) const;
        [[nodiscard]] bool isUnitary(double prec = 1e-12) const;
        [[nodiscard]] Gate adjoint() const;
    };
}