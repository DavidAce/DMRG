#pragma once
#include "math/tenx/fwd_decl.h"
#include "tools/common/log.h"
#include "tools/finite/opt.h"
#include "tools/finite/opt_mps.h"
#include <ceres/gradient_problem_solver.h>

/* clang-format off */
namespace tools::finite::opt::internal{

    using cplx = std::complex<double>;
    using real = double;

    template<typename T>
    using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    template<typename T, auto rank = 3>
    using TensorType = Eigen::Tensor<T, rank>;

    extern opt_mps lbfgs_optimize_variance    (const TensorsFinite & tensors, const opt_mps & initial_mps, const AlgorithmStatus & status, OptMeta & meta);
    extern opt_mps eigs_optimize_energy       (const TensorsFinite & tensors, const opt_mps & initial_mps, const AlgorithmStatus & status, OptMeta & meta);
    extern opt_mps eigs_optimize_subspace     (const TensorsFinite & tensors, const opt_mps & initial_mps, const AlgorithmStatus & status, OptMeta & meta);
    extern opt_mps eigs_optimize_variance     (const TensorsFinite & tensors, const opt_mps & initial_mps, const AlgorithmStatus & status, OptMeta & meta);
    extern opt_mps eigs_optimize_overlap      (const TensorsFinite & tensors, const opt_mps & initial_mps, const AlgorithmStatus & status, OptMeta & meta);
    extern void eigs_extract_results          (const TensorsFinite & tensors, const opt_mps & initial_mps, const OptMeta & meta, const eig::solver &solver,
                                               std::vector<opt_mps> &results, bool converged_only = true, double max_overlap_sq_sum = 0.7);


    namespace subspace{
        extern std::vector<int> generate_nev_list(int shape);

        template<typename Scalar>
        extern std::pair<Eigen::MatrixXcd, Eigen::VectorXd>
        find_subspace_part(const TensorsFinite & tensors, double energy_target, double target_subspace_error, const OptMeta & meta);

        template<typename Scalar>
        extern std::pair<Eigen::MatrixXcd, Eigen::VectorXd>
        find_subspace_full(const TensorsFinite & tensors);

        template<typename Scalar>
        extern std::vector<opt_mps>
        find_subspace(const TensorsFinite &tensors, double target_subspace_error, const OptMeta & meta);

        extern void
        filter_subspace(std::vector<opt_mps> & subspace, size_t max_accept);


        extern std::optional<size_t> get_idx_to_eigvec_with_highest_overlap(const std::vector<opt_mps> & eigvecs, double energy_llim_per_site, double energy_ulim_per_site);
        extern std::optional<size_t> get_idx_to_eigvec_with_lowest_variance(const std::vector<opt_mps> & eigvecs, double energy_llim_per_site, double energy_ulim_per_site);
        extern std::vector<size_t> get_idx_to_eigvec_with_highest_overlap(const std::vector<opt_mps> & eigvecs, size_t max_eigvecs, double energy_llim_per_site, double energy_ulim_per_site);

        extern Eigen::MatrixXcd get_eigvecs(const std::vector<opt_mps> & eigvecs);
        extern Eigen::VectorXd get_eigvals(const std::vector<opt_mps> & eigvecs);
        extern Eigen::VectorXd get_energies(const std::vector<opt_mps> & eigvecs);
        extern Eigen::VectorXd get_energies_per_site(const std::vector<opt_mps> & eigvecs);
        extern Eigen::VectorXd get_overlaps(const std::vector<opt_mps> & eigvecs);
        extern std::vector<double> get_subspace_errors(const std::vector<opt_mps> & eigvecs);
        extern double get_subspace_error(const std::vector<opt_mps> & eigvecs, std::optional<size_t> max_eigvecs = std::nullopt);
        extern double get_subspace_error(const std::vector<double> &overlaps);
        extern Eigen::VectorXcd get_vector_in_subspace(const std::vector<opt_mps> & eigvecs, size_t idx);
        extern Eigen::VectorXcd get_vector_in_subspace(const std::vector<opt_mps> & eigvecs, const Eigen::VectorXcd & subspace_vector);
        extern Eigen::VectorXcd get_vector_in_fullspace(const std::vector<opt_mps> & eigvecs, const Eigen::VectorXcd & subspace_vector);
        extern TensorType<cplx,3> get_tensor_in_fullspace(const std::vector<opt_mps> & eigvecs, const Eigen::VectorXcd & subspace_vector, const std::array<Eigen::Index,3> & dims);
        template <typename T>
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_hamiltonian_squared_in_subspace(const ModelFinite & model,
                                                                                                         const EdgesFinite & edges,
                                                                                                         const std::vector<opt_mps> & eigvecs);
    }

    inline ceres::GradientProblemSolver::Options lbfgs_default_options;

    inline bool no_state_in_window = false;



    extern double windowed_func_abs(double x,double window);
    extern double windowed_grad_abs(double x,double window);
    extern double windowed_func_pow(double x,double window);
    extern double windowed_grad_pow(double x,double window);
    extern std::pair<double,double> windowed_func_grad(double x,double window);
    extern long get_ops(long d, long chiL, long chiR, long m);
    extern long get_ops_R(long d, long chiL, long chiR, long m);



}
