#pragma once
#include <complex>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>
#include <config/enums.h>
#include <tensors/state/class_mps_site.h>

namespace tools::finite::opt {
    class opt_mps {
        using cplx = std::complex<double>;
        using real = double;

        private:
        // All of these values are supposed to be for the full system size
        std::optional<std::string>            name      = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 3>> tensor    = std::nullopt;
        std::optional<std::vector<size_t>>    sites     = std::nullopt;
        std::optional<double>                 eigval    = std::nullopt;
        std::optional<double>                 energy_r  = std::nullopt;
        std::optional<double>                 energy    = std::nullopt;
        std::optional<double>                 variance  = std::nullopt;
        std::optional<double>                 overlap   = std::nullopt;
        std::optional<double>                 alpha     = std::nullopt;
        std::optional<double>                 norm      = std::nullopt;
        std::optional<size_t>                 length    = std::nullopt;
        std::optional<size_t>                 iter      = std::nullopt;
        std::optional<size_t>                 counter   = std::nullopt;
        std::optional<double>                 time      = std::nullopt;
        std::optional<OptMode>                optMode   = std::nullopt;
        std::optional<OptSpace>               optSpace  = std::nullopt;

        public:
        bool is_basis_vector = false;
        std::vector<class_mps_site> mps_backup; // Used during subspace expansion to keep track of compatible neighbor mps

        opt_mps() = default;
        // Constructor used for candidates
        opt_mps(const std::string &name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double eigval_, double energy_reduced_,
                  std::optional<double> variance_, double overlap_, size_t length);
        // Constructor used for results
        opt_mps(const std::string &name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double energy_, double variance_,
                  double overlap_, size_t length, size_t iter_, size_t counter_, size_t time_);

        [[nodiscard]] const std::string &                get_name() const;
        [[nodiscard]] const Eigen::Tensor<cplx, 3> &     get_tensor() const;
        [[nodiscard]] Eigen::Map<const Eigen::VectorXcd> get_vector() const;
        [[nodiscard]] Eigen::Map<Eigen::VectorXd>        get_vector_cplx_as_2xreal();
        [[nodiscard]] Eigen::VectorXd                    get_vector_cplx_as_1xreal() const;
        [[nodiscard]] const std::vector<size_t> &        get_sites() const;
        [[nodiscard]] double                             get_energy() const;
        [[nodiscard]] double                             get_energy_reduced() const;
        [[nodiscard]] double                             get_energy_per_site() const;
        [[nodiscard]] double                             get_eigval() const;
        [[nodiscard]] double                             get_variance() const;
        [[nodiscard]] double                             get_variance_per_site() const;
        [[nodiscard]] double                             get_overlap() const;
        [[nodiscard]] double                             get_alpha() const;
        [[nodiscard]] double                             get_norm() const;
        [[nodiscard]] size_t                             get_length() const;
        [[nodiscard]] size_t                             get_iter() const;
        [[nodiscard]] size_t                             get_counter() const;
        [[nodiscard]] double                             get_time() const;
        [[nodiscard]] OptSpace                           get_optspace() const;
        [[nodiscard]] OptMode                            get_optmode() const;
        void                                             clear();
        void                                             normalize();
        void                                             set_name(const std::string &name_);
        void                                             set_tensor(const Eigen::Tensor<cplx, 3> &tensor_);
        void                                             set_tensor(const Eigen::VectorXcd &vector, const Eigen::DSizes<long, 3> &dims);
        void                                             set_sites(const std::vector<size_t> &sites_);
        void                                             set_eigval(double eigval_);
        void                                             set_energy_reduced(double energy_reduced_);
        void                                             set_energy(double energy_);
        void                                             set_energy_per_site(double energy_per_site_);
        void                                             set_variance(double variance_);
        void                                             set_variance_per_site(double variance_per_site_);
        void                                             set_overlap(double overlap_);
        void                                             set_alpha(std::optional<double> alpha_);
        void                                             set_length(size_t length);
        void                                             set_iter(size_t iter_);
        void                                             set_counter(size_t counter_);
        void                                             set_time(double time_);
        void                                             set_tensor_cplx(const double *data, const Eigen::DSizes<long, 3> &dims);
        void                                             set_tensor_real(const double *data, const Eigen::DSizes<long, 3> &dims);
        void                                             set_optspace(OptSpace optspace_);
        void                                             set_optmode(OptMode optmode_);
        void                                             validate_candidate() const;
        void                                             validate_result() const;
        bool                                             operator<(const opt_mps &rhs) const;
        bool                                             operator>(const opt_mps &rhs) const;
        bool                                             has_nan() const;
    };
}
