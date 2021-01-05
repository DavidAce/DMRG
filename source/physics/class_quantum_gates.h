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

    class Gate {
        private:
        Eigen::Tensor<Scalar,2> exp_internal(const Eigen::Tensor<Scalar,2> & op_, Scalar alpha) const;
        mutable std::optional<Eigen::Tensor<Scalar,2>> adj = std::nullopt;
        public:
        Eigen::Tensor<Scalar,2> op;
        std::vector<size_t> pos;
        Gate(const Eigen::Tensor<Scalar,2>  & op_, const std::vector<size_t> & pos_);
        Gate(const Eigen::Tensor<Scalar,2>  & op_, const std::vector<size_t> & pos_, Scalar alpha);
        void exp_inplace(Scalar alpha);
        [[nodiscard]] Gate exp(Scalar alpha) const;
        [[nodiscard]] bool isUnitary(double prec = 1e-12) const;
        [[nodiscard]] Eigen::Tensor<Scalar,2>& adjoint() const;
    };
}