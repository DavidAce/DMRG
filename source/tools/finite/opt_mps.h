#pragma once
#include <complex>
#include <config/enums.h>
#include <optional>
#include <tensors/state/class_mps_site.h>
#include <unsupported/Eigen/CXX11/Tensor>

namespace tools::finite::opt {
    class opt_mps {
        using cplx = std::complex<double>;
        using real = double;

        private:
        // All of these values are supposed to be for the full system size
        std::optional<std::string>            name          = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 3>> tensor        = std::nullopt;
        std::optional<std::vector<size_t>>    sites         = std::nullopt;
        std::optional<double>                 eigval        = std::nullopt;
        std::optional<double>                 energy_r      = std::nullopt;
        std::optional<double>                 energy        = std::nullopt;
        std::optional<double>                 variance      = std::nullopt;
        std::optional<double>                 overlap       = std::nullopt;
        std::optional<double>                 alpha         = std::nullopt;
        std::optional<double>                 norm          = std::nullopt;
        std::optional<size_t>                 length        = std::nullopt;
        std::optional<size_t>                 iter          = std::nullopt;
        std::optional<size_t>                 counter       = std::nullopt;
        std::optional<double>                 time          = std::nullopt;
        std::optional<double>                 delta_f       = std::nullopt;
        std::optional<double>                 grad_norm     = std::nullopt;
        std::optional<double>                 relchange     = std::nullopt;
        std::optional<long>                   krylov_nev    = std::nullopt;
        std::optional<long>                   krylov_ncv    = std::nullopt;
        std::optional<double>                 krylov_tol    = std::nullopt;
        std::optional<double>                 krylov_eigval = std::nullopt;
        std::optional<std::string>            krylov_ritz   = std::nullopt;
        std::optional<cplx>                   krylov_shift  = std::nullopt;
        std::optional<OptMode>                optMode       = std::nullopt;
        std::optional<OptSpace>               optSpace      = std::nullopt;
        std::optional<OptExit>                optExit       = std::nullopt;

        public:
        bool                        is_basis_vector = false;
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
        [[nodiscard]] double                             get_delta_f() const;
        [[nodiscard]] double                             get_grad_norm() const;
        [[nodiscard]] double                             get_relchange() const;
        [[nodiscard]] long                               get_krylov_nev() const;
        [[nodiscard]] long                               get_krylov_ncv() const;
        [[nodiscard]] double                             get_krylov_tol() const;
        [[nodiscard]] double                             get_krylov_eigval() const;
        [[nodiscard]] std::string                        get_krylov_ritz() const;
        [[nodiscard]] cplx                               get_krylov_shift() const;
        [[nodiscard]] OptSpace                           get_optspace() const;
        [[nodiscard]] OptMode                            get_optmode() const;
        [[nodiscard]] OptExit                            get_optexit() const;
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
        void                                             set_delta_f(double delta_f_);
        void                                             set_grad_norm(double grad_norm_);
        void                                             set_relchange(double relative_change_);
        void                                             set_krylov_nev(long nev_);
        void                                             set_krylov_ncv(long ncv_);
        void                                             set_krylov_tol(double tol_);
        void                                             set_krylov_eigval(double krylov_eigval_);
        void                                             set_krylov_ritz(const std::string &ritz_);
        void                                             set_krylov_shift(const cplx &ritz_);
        void                                             set_tensor_cplx(const double *data, const Eigen::DSizes<long, 3> &dims);
        void                                             set_tensor_real(const double *data, const Eigen::DSizes<long, 3> &dims);
        void                                             set_optspace(OptSpace optspace_);
        void                                             set_optmode(OptMode optmode_);
        void                                             set_optexit(OptExit optexit_);
        void                                             validate_candidate() const;
        void                                             validate_result() const;
        bool                                             operator<(const opt_mps &rhs) const;
        bool                                             operator>(const opt_mps &rhs) const;
        bool                                             has_nan() const;
    };
}
