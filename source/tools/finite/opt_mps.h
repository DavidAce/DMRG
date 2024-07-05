#pragma once
#include "config/enums.h"
#include "math/float.h"
#include "tensors/site/mps/MpsSite.h"
#include <complex>
#include <optional>
#include <unsupported/Eigen/CXX11/Tensor>

namespace tools::finite::opt {
    class opt_mps {
        private:
        struct NextEv {
            // If the eigenvalues are  x_0, x_1, x_2 ... x_n, this struct saves
            // information about the eigensolution x_1,x_2..., when nev > 1
            // We save information about the next neighboring eigensolution
            long   eigs_idx;
            double eigs_eigval;
            double eigs_rnorm;
            double energy;
            double variance;
            double overlap;
                   NextEv(const opt_mps &res);
        };
        // All of these values are supposed to be for the full system size
        std::optional<std::string>            name             = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 3>> tensor           = std::nullopt;
        std::optional<Eigen::Tensor<cplx, 2>> bond             = std::nullopt;
        std::optional<std::vector<size_t>>    sites            = std::nullopt;
        std::optional<double>                 eshift           = std::nullopt; /*!< current energy shift in the energy MPOs: (H-eshift)  */
        std::optional<double>                 eshift_eigval    = std::nullopt; /*!< eigenvalue of (H-eshift)  */
        std::optional<double>                 eshift_eigval_sq = std::nullopt; /*!< eigenvalue of (H-eshift)²  */
        std::optional<double>                 energy           = std::nullopt; /*!< Energy: eigenvalue of H: (E - eshift) + eshift   */
        std::optional<double>                 variance         = std::nullopt; /*!< Variance: H²-E² = <(H-eshift)²> - <H-eshift>² */
        std::optional<double>                 rnorm            = std::nullopt;
        std::optional<double>                 overlap          = std::nullopt;
        std::optional<double>                 alpha            = std::nullopt;
        std::optional<double>                 norm             = std::nullopt;
        std::optional<size_t>                 length           = std::nullopt;
        std::optional<size_t>                 iter             = std::nullopt;
        std::optional<size_t>                 num_op           = std::nullopt; /*!< Number of inverse-matrix-vector products */
        std::optional<size_t>                 num_mv           = std::nullopt; /*!< Number of matrix-vector products */
        std::optional<size_t>                 num_pc           = std::nullopt; /*!< Number of preconditioner calls */
        std::optional<double>                 time             = std::nullopt;
        std::optional<double>                 time_mv          = std::nullopt;
        std::optional<double>                 delta_f          = std::nullopt;
        std::optional<double>                 grad_tol         = std::nullopt;
        std::optional<double>                 grad_max         = std::nullopt;
        std::optional<double>                 relchange        = std::nullopt;
        std::optional<long>                   eigs_idx         = std::nullopt;
        std::optional<long>                   eigs_nev         = std::nullopt;
        std::optional<long>                   eigs_ncv         = std::nullopt;
        std::optional<double>                 eigs_tol         = std::nullopt;
        std::optional<double>                 eigs_eigval      = std::nullopt;
        std::optional<std::string>            eigs_ritz        = std::nullopt;
        std::optional<cplx>                   eigs_shift       = std::nullopt;
        std::optional<OptCost>                optCost          = std::nullopt;
        std::optional<OptAlgo>                optAlgo          = std::nullopt;
        std::optional<OptSolver>              optSolver        = std::nullopt;
        std::optional<OptRitz>                optRitz          = std::nullopt;
        std::optional<OptExit>                optExit          = std::nullopt;
        std::optional<long>                   bond_lim         = std::nullopt;
        std::optional<double>                 trnc_lim         = std::nullopt;

        public:
        bool                 is_basis_vector = false;
        std::vector<MpsSite> mps_backup; // Used during subspace expansion to keep track of compatible neighbor mps
        std::vector<NextEv>  next_evs;

        opt_mps() = default;
        // Constructor used for candidates
        opt_mps(std::string_view name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double eshift_eigval_, double eshift_,
                std::optional<double> variance_, double overlap_, size_t length);
        // Constructor used for results
        opt_mps(std::string_view name_, const Eigen::Tensor<cplx, 3> &tensor_, const std::vector<size_t> &sites_, double energy_, double variance_,
                double overlap_, size_t length, size_t iter_, size_t counter_, size_t time_);

        [[nodiscard]] bool                               is_initialized() const;
        [[nodiscard]] std::string_view                   get_name() const;
        [[nodiscard]] const Eigen::Tensor<cplx, 3>      &get_tensor() const;
        [[nodiscard]] const Eigen::Tensor<cplx, 2>      &get_bond() const;
        [[nodiscard]] Eigen::Map<const Eigen::VectorXcd> get_vector() const;
        [[nodiscard]] Eigen::Map<const Eigen::VectorXd>  get_vector_cplx_as_2xreal() const;
        [[nodiscard]] Eigen::VectorXd                    get_vector_cplx_as_1xreal() const;

        template<OptType optType>
        [[nodiscard]] Eigen::VectorXd get_initial_state_with_lagrange_multiplier() const;

        template<OptType optType>
        [[nodiscard]] std::vector<double> get_stl_initial_state_with_lagrange_multiplier() const;

        [[nodiscard]] const std::vector<size_t> &get_sites() const;
        [[nodiscard]] double                     get_energy() const;
        [[nodiscard]] double                     get_energy_shift() const;
        [[nodiscard]] double                     get_energy_per_site() const;
        [[nodiscard]] double                     get_eshift_eigval() const;
        [[nodiscard]] double                     get_variance() const;
        [[nodiscard]] double                     get_rnorm() const;
        [[nodiscard]] double                     get_overlap() const;
        [[nodiscard]] double                     get_alpha() const;
        [[nodiscard]] double                     get_norm() const;
        [[nodiscard]] size_t                     get_length() const;
        [[nodiscard]] size_t                     get_iter() const;
        [[nodiscard]] size_t                     get_op() const;
        [[nodiscard]] size_t                     get_mv() const;
        [[nodiscard]] size_t                     get_pc() const;
        [[nodiscard]] double                     get_time() const;
        [[nodiscard]] double                     get_time_mv() const;
        [[nodiscard]] double                     get_delta_f() const;
        [[nodiscard]] double                     get_grad_tol() const;
        [[nodiscard]] double                     get_grad_max() const;
        [[nodiscard]] double                     get_relchange() const;
        [[nodiscard]] long                       get_eigs_idx() const;
        [[nodiscard]] long                       get_eigs_nev() const;
        [[nodiscard]] long                       get_eigs_ncv() const;
        [[nodiscard]] double                     get_eigs_tol() const;
        [[nodiscard]] double                     get_eigs_eigval() const;
        [[nodiscard]] std::string_view           get_eigs_ritz() const;
        [[nodiscard]] cplx                       get_eigs_shift() const;
        [[nodiscard]] OptSolver                  get_optsolver() const;
        [[nodiscard]] OptCost                    get_optcost() const;
        [[nodiscard]] OptAlgo                    get_optalgo() const;
        [[nodiscard]] OptExit                    get_optexit() const;
        [[nodiscard]] long                       get_bond_lim() const;
        [[nodiscard]] double                     get_trnc_lim() const;
        [[nodiscard]] std::vector<NextEv>        get_next_evs() const;

        void clear();
        void normalize();
        void set_name(std::string_view name_);
        void set_tensor(const Eigen::Tensor<cplx, 3> &tensor_);
        void set_tensor(const Eigen::VectorXcd &vector, const Eigen::DSizes<long, 3> &dims);
        void set_bond(const Eigen::MatrixXcd &matrix);
        void set_bond(const Eigen::Tensor<cplx, 2> &bond_);
        void set_sites(const std::vector<size_t> &sites_);
        void set_eshift_eigval(double eshift_eigval_);
        void set_energy_shift(double energy_shift_);
        void set_energy(double energy_);
        void set_energy_per_site(double energy_per_site_);
        void set_variance(double variance_);
        void set_rnorm(const double rnorm_);
        void set_overlap(double overlap_);
        void set_alpha(std::optional<double> alpha_);
        void set_length(size_t length);
        void set_iter(size_t iter_);
        void set_op(size_t op_);
        void set_mv(size_t mv_);
        void set_pc(size_t pc_);
        void set_time(double time_);
        void set_time_mv(double time_mv_);
        void set_delta_f(double delta_f_);
        void set_grad_tol(double grad_tol_);
        void set_grad_max(double grad_max_);
        void set_relchange(double relative_change_);
        void set_eigs_idx(long idx_);
        void set_eigs_nev(long nev_);
        void set_eigs_ncv(long ncv_);
        void set_eigs_tol(double tol_);
        void set_eigs_eigval(double eigval_);
        void set_eigs_ritz(std::string_view ritz_);
        void set_eigs_shift(const cplx shift_);
        void set_tensor_cplx(const double *data, const Eigen::DSizes<long, 3> &dims);
        void set_tensor_real(const double *data, const Eigen::DSizes<long, 3> &dims);
        void set_optsolver(OptSolver optsolver_);
        void set_optcost(OptCost optcost_);
        void set_optalgo(OptAlgo optalgo_);
        void set_optexit(OptExit optexit_);
        void set_bond_limit(long bond_);
        void set_trnc_limit(double trnc_);
        void validate_initial_mps() const;
        void validate_basis_vector() const;
        void validate_result() const;
        bool operator<(const opt_mps &rhs) const;
        bool operator>(const opt_mps &rhs) const;
        bool has_nan() const;
    };
}
