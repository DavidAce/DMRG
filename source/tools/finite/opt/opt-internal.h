#pragma once
// This trick avoids collision between the preprocessor
// symbol I and a template type size_t I in Eigen getting
// overridden and causing trouble at compile time.
//#include <complex.h>
//#undef I
#include "../opt.h"
#include <ceres/gradient_problem_solver.h>
#include <general/eigen_tensor_fwd_decl.h>
#include <tools/common/log.h>

/* clang-format off */
namespace tools::finite::opt::internal{

    using cplx = std::complex<double>;
    using real = double;

    template<typename T>
    using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    template<typename T, auto rank = 3>
    using TensorType = Eigen::Tensor<T, rank>;

    extern opt_mps
    ceres_direct_optimization(const TensorsFinite & tensors, const AlgorithmStatus &status, OptType optType, OptMode optMode,OptSpace optSpace);

    extern opt_mps
    ceres_direct_optimization(const TensorsFinite & tensors, const opt_mps &initial_mps,
                              const AlgorithmStatus &status, OptType optType, OptMode optMode,OptSpace optSpace);

    extern opt_mps
    ceres_subspace_optimization (const TensorsFinite & tensors, const AlgorithmStatus & status, OptType optType, OptMode optMode,OptSpace optSpace);

    extern opt_mps
    ceres_subspace_optimization(const TensorsFinite & tensors, const opt_mps &initial_mps,
                                const AlgorithmStatus &status, OptType optType, OptMode optMode,OptSpace optSpace);

    extern opt_mps
    ceres_optimize_subspace(const TensorsFinite & tensors, const opt_mps &initial_mps, const std::vector<opt_mps> & candidate_list,
                            const Eigen::MatrixXcd &H2_subspace, const AlgorithmStatus &status, OptType optType, OptMode optMode,OptSpace optSpace);


    extern Eigen::Tensor<std::complex<double>,3> cppoptlib_optimization      (const TensorsFinite & tensors, const AlgorithmStatus & status);
    extern opt_mps ground_state_optimization  (const TensorsFinite & tensors, const AlgorithmStatus & status, StateRitz ritz);
    extern opt_mps ground_state_optimization  (const opt_mps & initial_mps,const TensorsFinite & tensors, const AlgorithmStatus & status, std::string_view ritz);
    extern opt_mps krylov_energy_optimization (const TensorsFinite & tensors, const opt_mps &initial_mps,
                                                 const AlgorithmStatus &status, OptType optType, OptMode optMode,
                                                 OptSpace optSpace);
    extern opt_mps krylov_variance_optimization(const TensorsFinite &tensors, const opt_mps &initial_mps,
                                                  const AlgorithmStatus &status, OptType optType, OptMode optMode,
                                                  OptSpace optSpace);
    extern void krylov_extract_solutions(const opt_mps &initial_mps,
                                         const TensorsFinite &tensors,
                                         eig::solver &solver,
                                         std::vector<tools::finite::opt::opt_mps> &eigvecs_mps,
                                         const std::string & tag = "",
                                         bool converged_only = true);
    extern opt_mps primme_variance_optimization(const TensorsFinite &tensors, const opt_mps &initial_mps,
                                                const AlgorithmStatus &status, OptType optType, OptMode optMode,
                                                OptSpace optSpace);
    extern Eigen::Tensor<std::complex<double>,3> ham_sq_optimization         (const TensorsFinite & tensors, OptType optType, OptMode optMode, OptSpace optSpace);
    extern Eigen::Tensor<std::complex<double>,3> ceres_rosenbrock_optimization (const StateFinite & state);


    namespace subspace{
        template<typename Scalar>
        extern std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
        find_subspace_part(const MatrixType<Scalar> & H_local, const TensorType<cplx,3> &multisite_mps, double energy_target, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

        template<typename Scalar>
        extern std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
        find_subspace_full(const MatrixType<Scalar> & H_local, const TensorType<cplx,3> &multisite_tensor);

        template<typename Scalar>
        extern std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
        find_subspace(const TensorsFinite &tensors, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

        template<typename Scalar>
        extern std::vector<opt_mps>
        find_candidates(const TensorsFinite &tensors, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

        extern void
        filter_candidates(std::vector<opt_mps> & candidate_list, double maximum_subspace_error, size_t max_accept);


        extern std::optional<size_t> get_idx_to_candidate_with_highest_overlap(const std::vector<opt_mps> & candidate_list, double energy_llim_per_site, double energy_ulim_per_site);
        extern std::optional<size_t> get_idx_to_candidate_with_lowest_variance(const std::vector<opt_mps> & candidate_list, double energy_llim_per_site, double energy_ulim_per_site);
        extern std::vector<size_t> get_idx_to_candidates_with_highest_overlap(const std::vector<opt_mps> & candidate_list, size_t max_candidates, double energy_llim_per_site, double energy_ulim_per_site);

        extern Eigen::MatrixXcd get_eigvecs(const std::vector<opt_mps> & candidate_list);
        extern Eigen::VectorXd get_eigvals(const std::vector<opt_mps> & candidate_list);
        extern Eigen::VectorXd get_energies(const std::vector<opt_mps> & candidate_list);
        extern Eigen::VectorXd get_energies_per_site(const std::vector<opt_mps> & candidate_list);
        extern Eigen::VectorXd get_overlaps(const std::vector<opt_mps> & candidate_list);
        extern std::vector<double> get_subspace_errors(const std::vector<opt_mps> & candidate_list);
        extern double get_subspace_error(const std::vector<opt_mps> & candidate_list, std::optional<size_t> max_candidates = std::nullopt);
        extern double get_subspace_error(const std::vector<double> &overlaps);
        extern Eigen::VectorXcd get_vector_in_subspace(const std::vector<opt_mps> & candidate_list, size_t idx);
        extern Eigen::VectorXcd get_vector_in_subspace(const std::vector<opt_mps> & candidate_list, const Eigen::VectorXcd & subspace_vector);
        extern Eigen::VectorXcd get_vector_in_fullspace(const std::vector<opt_mps> & candidate_list, const Eigen::VectorXcd & subspace_vector);
        extern TensorType<cplx,3> get_tensor_in_fullspace(const std::vector<opt_mps> & candidate_list, const Eigen::VectorXcd & subspace_vector, const std::array<Eigen::Index,3> & dims);
    }

    template <typename T>
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multisite_hamiltonian_matrix(const ModelFinite & model, const EdgesFinite & edges, double energy_shift = 0);
    template <typename T>
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multisite_hamiltonian_squared_matrix(const ModelFinite & model, const EdgesFinite & edges);
    template <typename T>
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multisite_hamiltonian_squared_subspace_matrix(const ModelFinite & model,
                                                                                                     const EdgesFinite & edges,
                                                                                                     const std::vector<opt_mps> & candidate_list,
                                                                                                     double energy_shift = 0.0);


    inline ceres::GradientProblemSolver::Options ceres_default_options;

    inline bool no_state_in_window = false;

    extern std::vector<int> generate_size_list(int shape);


    template <typename T>  int sgn(const T val) {return (T(0) < val) - (val < T(0)); }
    extern double windowed_func_abs(double x,double window);
    extern double windowed_grad_abs(double x,double window);
    extern double windowed_func_pow(double x,double window);
    extern double windowed_grad_pow(double x,double window);
    extern std::pair<double,double> windowed_func_grad(double x,double window);
    extern long get_ops(long d, long chiL, long chiR, long m);
    extern long get_ops_R(long d, long chiL, long chiR, long m);



}
