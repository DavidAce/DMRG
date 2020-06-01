#pragma once
#include <complex>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>

namespace tools::finite::opt::internal {
    class candidate_tensor {
        using cplx = std::complex<double>;
        using real = double;

        private:
        std::optional<Eigen::Tensor<cplx, 3>> tensor   = std::nullopt;
        std::optional<double>                 eigval   = std::nullopt;
        std::optional<double>                 energy_r = std::nullopt;
        std::optional<double>                 energy   = std::nullopt;
        std::optional<double>                 variance = std::nullopt;
        std::optional<double>                 overlap  = std::nullopt;

        public:
        bool is_basis_vector = false;

        candidate_tensor() = default;
        candidate_tensor(const Eigen::Tensor<cplx, 3> &tensor_, double eigval, double energy_reduced, std::optional<double> variance_, double overlap_);

        const Eigen::Tensor<cplx, 3> &     get_tensor() const;
        Eigen::Map<const Eigen::VectorXcd> get_vector() const;
        double                             get_energy() const;
        double                             get_eigval() const;
        double                             get_variance() const;
        double                             get_overlap() const;
        void                               clear();
        void                               set_variance(double var);

        bool operator<(const candidate_tensor &rhs) const;
        bool operator>(const candidate_tensor &rhs) const;
    };
}
