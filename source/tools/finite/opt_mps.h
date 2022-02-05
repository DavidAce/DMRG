#pragma once
#include <complex>
#include <config/enums.h>
#include <optional>
#include <tensors/site/mps/MpsSite.h>
#include <unsupported/Eigen/CXX11/Tensor>

namespace tools::finite::opt {
    class opt_mps {
        using cplx = std::complex<double>;
        using real = double;

        private:
        // All of these values are supposed to be for the full system size
        std::optional<std::string>            name        = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 3>> tensor      = std::nullopt;
        std::optional<std::vector<size_t>>    sites       = std::nullopt;
        std::optional<double>                 eigval      = std::nullopt;
        std::optional<double>                 energy_r    = std::nullopt;
        std::optional<double>                 energy      = std::nullopt;
        std::optional<double>                 variance    = std::nullopt;
        std::optional<double>                 overlap     = std::nullopt;
        std::optional<double>                 alpha       = std::nullopt;
        std::optional<double>                 norm        = std::nullopt;
        std::optional<size_t>                 length      = std::nullopt;
        std::optional<size_t>                 iter        = std::nullopt;
        std::optional<size_t>                 num_op      = std::nullopt;
        std::optional<size_t>                 num_mv      = std::nullopt;
        std::optional<double>                 time        = std::nullopt;
        std::optional<double>                 delta_f     = std::nullopt;
        std::optional<double>                 max_grad    = std::nullopt;
        std::optional<double>                 relchange   = std::nullopt;
        std::optional<long>                   eigs_idx    = std::nullopt;
        std::optional<long>                   eigs_nev    = std::nullopt;
        std::optional<long>                   eigs_ncv    = std::nullopt;
        std::optional<double>                 eigs_tol    = std::nullopt;
        std::optional<double>                 eigs_eigval = std::nullopt;
        std::optional<std::string>            eigs_ritz   = std::nullopt;
        std::optional<cplx>                   eigs_shift  = std::nullopt;
        std::optional<double>                 eigs_resid  = std::nullopt;
        std::optional<OptMode>                optMode     = std::nullopt;
        std::optional<OptSolver>              optSolver   = std::nullopt;
        std::optional<OptRitz>                optRitz     = std::nullopt;
        std::optional<OptExit>                optExit     = std::nullopt;

        public:
        bool                 is_basis_vector = false;
        std::vector<MpsSite> mps_backup; // Used during subspace expansion to keep track of compatible neighbor mps

        opt_mps() = default;
        // Constructor used for candidates
        opt_mps(std::string_view name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double eigval_, double energy_shift_,
                std::optional<double> variance_, double overlap_, size_t length);
        // Constructor used for results
        opt_mps(std::string_view name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double energy_, double variance_,
                double overlap_, size_t length, size_t iter_, size_t counter_, size_t time_);

        [[nodiscard]] std::string_view                   get_name() const;
        [[nodiscard]] const Eigen::Tensor<cplx, 3>      &get_tensor() const;
        [[nodiscard]] Eigen::Map<const Eigen::VectorXcd> get_vector() const;
        [[nodiscard]] Eigen::Map<Eigen::VectorXd>        get_vector_cplx_as_2xreal();
        [[nodiscard]] Eigen::VectorXd                    get_vector_cplx_as_1xreal() const;

        template<OptType optType>
        [[nodiscard]] Eigen::VectorXd get_initial_state_with_lagrange_multiplier() const;

        [[nodiscard]] const std::vector<size_t> &get_sites() const;
        [[nodiscard]] double                     get_energy() const;
        [[nodiscard]] double                     get_energy_shift() const;
        [[nodiscard]] double                     get_energy_per_site() const;
        [[nodiscard]] double                     get_eigval() const;
        [[nodiscard]] double                     get_variance() const;
        [[nodiscard]] double                     get_variance_per_site() const;
        [[nodiscard]] double                     get_overlap() const;
        [[nodiscard]] double                     get_alpha() const;
        [[nodiscard]] double                     get_norm() const;
        [[nodiscard]] size_t                     get_length() const;
        [[nodiscard]] size_t                     get_iter() const;
        [[nodiscard]] size_t                     get_op() const;
        [[nodiscard]] size_t                     get_mv() const;
        [[nodiscard]] double                     get_time() const;
        [[nodiscard]] double                     get_delta_f() const;
        [[nodiscard]] double                     get_max_grad() const;
        [[nodiscard]] double                     get_relchange() const;
        [[nodiscard]] long                       get_eigs_idx() const;
        [[nodiscard]] long                       get_eigs_nev() const;
        [[nodiscard]] long                       get_eigs_ncv() const;
        [[nodiscard]] double                     get_eigs_tol() const;
        [[nodiscard]] double                     get_eigs_eigval() const;
        [[nodiscard]] std::string_view           get_eigs_ritz() const;
        [[nodiscard]] cplx                       get_eigs_shift() const;
        [[nodiscard]] double                     get_eigs_resid() const;
        [[nodiscard]] OptSolver                  get_optsolver() const;
        [[nodiscard]] OptMode                    get_optmode() const;
        [[nodiscard]] OptExit                    get_optexit() const;
        void                                     clear();
        void                                     normalize();
        void                                     set_name(std::string_view name_);
        void                                     set_tensor(const Eigen::Tensor<cplx, 3> &tensor_);
        void                                     set_tensor(const Eigen::VectorXcd &vector, const Eigen::DSizes<long, 3> &dims);
        void                                     set_sites(const std::vector<size_t> &sites_);
        void                                     set_eigval(double eigval_);
        void                                     set_energy_shift(double energy_shift_);
        void                                     set_energy(double energy_);
        void                                     set_energy_per_site(double energy_per_site_);
        void                                     set_variance(double variance_);
        void                                     set_variance_per_site(double variance_per_site_);
        void                                     set_overlap(double overlap_);
        void                                     set_alpha(std::optional<double> alpha_);
        void                                     set_length(size_t length);
        void                                     set_iter(size_t iter_);
        void                                     set_op(size_t op_);
        void                                     set_mv(size_t mv_);
        void                                     set_time(double time_);
        void                                     set_delta_f(double delta_f_);
        void                                     set_max_grad(double grad_norm_);
        void                                     set_relchange(double relative_change_);
        void                                     set_eigs_idx(long idx_);
        void                                     set_eigs_nev(long nev_);
        void                                     set_eigs_ncv(long ncv_);
        void                                     set_eigs_tol(double tol_);
        void                                     set_eigs_eigval(double eigval_);
        void                                     set_eigs_ritz(std::string_view ritz_);
        void                                     set_eigs_shift(const cplx &shift_);
        void                                     set_eigs_resid(const double &resid);
        void                                     set_tensor_cplx(const double *data, const Eigen::DSizes<long, 3> &dims);
        void                                     set_tensor_real(const double *data, const Eigen::DSizes<long, 3> &dims);
        void                                     set_optsolver(OptSolver optspace_);
        void                                     set_optmode(OptMode optmode_);
        void                                     set_optexit(OptExit optexit_);
        void                                     validate_basis_vector() const;
        void                                     validate_result() const;
        bool                                     operator<(const opt_mps &rhs) const;
        bool                                     operator>(const opt_mps &rhs) const;
        bool                                     has_nan() const;
    };
}
