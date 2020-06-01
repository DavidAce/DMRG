//
// Created by david on 2019-03-18.
//

#pragma once
// This trick avoids collision between the preprocessor
// symbol I and a template type size_t I in Eigen getting
// overridden and causing trouble at compile time.
#include <complex.h>
#undef I
#include <ceres/ceres.h>
#include <config/enums.h>
#include <general/nmspc_tensor_omp.h>
#include <tools/finite/opt-internals/candidate_tensor.h>
#include <tools/finite/opt-internals/enum_classes.h>
class class_state_finite;
class class_model_finite;
class class_edges_finite;
class class_tensors_finite;
class class_algorithm_status;

namespace tools::finite::opt {
    using Scalar = std::complex<double>;
    extern Eigen::Tensor<Scalar, 3> find_excited_state(const class_tensors_finite &tensors, const class_algorithm_status &status, OptMode optMode,
                                                       OptSpace optSpace, OptType optType);
    extern Eigen::Tensor<Scalar, 3> find_ground_state(const class_tensors_finite &tensors, StateRitz ritz);
}

/* clang-format off */
namespace tools::finite::opt::internal{
        extern Eigen::Tensor<std::complex<double>,3>
               ceres_direct_optimization(const class_tensors_finite & tensors, const class_algorithm_status &status, OptType optType, OptMode optMode,OptSpace optSpace);

        extern Eigen::Tensor<std::complex<double>,3>
               ceres_direct_optimization(const class_tensors_finite & tensors, const Eigen::Tensor<std::complex<double>,3> &multisite_tensor,
                                         const class_algorithm_status &status, OptType optType, OptMode optMode,OptSpace optSpace);

        extern Eigen::Tensor<std::complex<double>,3>
               ceres_subspace_optimization (const class_tensors_finite & tensors, const class_algorithm_status & status, OptType optType, OptMode optMode,OptSpace optSpace);

        extern Eigen::Tensor<std::complex<double>,3> cppoptlib_optimization      (const class_tensors_finite & tensors, const class_algorithm_status & status);
        extern Eigen::Tensor<std::complex<double>,3> ground_state_optimization   (const class_tensors_finite & tensors, StateRitz ritz);
        extern Eigen::Tensor<std::complex<double>,3> ground_state_optimization   (const class_tensors_finite & tensors, std::string_view ritz);
        extern Eigen::Tensor<std::complex<double>,3> ham_sq_optimization         (const class_tensors_finite & tensors, OptType optType, OptMode optMode, OptSpace optSpace);
        extern Eigen::Tensor<std::complex<double>,3> ceres_rosenbrock_optimization (const class_state_finite & state);

        namespace subspace{
            using cplx = std::complex<double>;
            using real = double;

            template<typename T>
            using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
            template<typename T, auto rank = 3>
            using TensorType = Eigen::Tensor<T, rank>;



            template<typename Scalar>
            extern std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
            find_subspace_part(const MatrixType<Scalar> & H_local, const TensorType<cplx,3> &multisite_tensor, double energy_target, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

            template<typename Scalar>
            extern std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
            find_subspace_full(const MatrixType<Scalar> & H_local, const TensorType<cplx,3> &multisite_tensor);

            template<typename Scalar>
            extern std::tuple<Eigen::MatrixXcd, Eigen::VectorXd>
            find_subspace(const class_tensors_finite &tensors, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

            template<typename Scalar>
            extern std::vector<candidate_tensor>
            find_candidates(const class_tensors_finite &tensors, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

            extern void
            filter_candidates(std::vector<candidate_tensor> & candidate_list, double maximum_subspace_error, size_t max_accept);

            extern std::tuple<Eigen::MatrixXcd, Eigen::VectorXd, Eigen::VectorXd, double>
            filter_states(const Eigen::MatrixXcd &eigvecs, const Eigen::VectorXd &eigvals,
                          Eigen::VectorXd &overlaps, double maximum_subspace_error,
                          size_t max_accept);

            extern std::pair<double, int> get_best_variance_in_window(const class_tensors_finite &tensors, const Eigen::MatrixXcd &eigvecs,
                                                                      const Eigen::VectorXd &energies_per_site, double lbound, double ubound);

            extern std::pair<double, int> get_best_overlap_in_window(const Eigen::VectorXd &overlaps, const Eigen::VectorXd &energies_per_site, double lbound, double ubound);
            extern std::optional<size_t> get_idx_to_candidate_with_highest_overlap(const std::vector<candidate_tensor> & candidate_list, double energy_lbound, double energy_ubound);
            extern std::vector<size_t> get_idx_to_candidates_with_highest_overlap(const std::vector<candidate_tensor> & candidate_list, size_t max_candidates, double energy_lbound, double energy_ubound);

            extern std::vector<std::pair<double, int>> get_best_candidates_in_window(const Eigen::VectorXd &overlaps, const Eigen::VectorXd &energies_per_site, double lbound,double ubound);

            extern Eigen::MatrixXcd get_eigvecs(const std::vector<candidate_tensor> & candidate_list);
            extern Eigen::VectorXd get_eigvals(const std::vector<candidate_tensor> & candidate_list);
            extern Eigen::VectorXd get_energies(const std::vector<candidate_tensor> & candidate_list);
            extern Eigen::VectorXd get_overlaps(const std::vector<candidate_tensor> & candidate_list);
            extern double get_subspace_error(const std::vector<candidate_tensor> & candidate_list);
            extern double get_subspace_error(const std::vector<double> &overlaps);

            extern Eigen::VectorXcd get_vector_in_subspace(const std::vector<candidate_tensor> & candidate_list, size_t idx);
            extern Eigen::VectorXcd get_vector_in_fullspace(const std::vector<candidate_tensor> & candidate_list, const Eigen::VectorXcd & subspace_vector);
        }


        namespace local_hamiltonians{
            extern Eigen::Tensor<std::complex<double>,6>   get_multi_hamiltonian_tensor(const class_model_finite & model, const class_edges_finite & edges);
            extern Eigen::Tensor<std::complex<double>,6>   get_multi_hamiltonian_squared_tensor(const class_model_finite & model, const class_edges_finite & edges);
            extern Eigen::MatrixXcd                        get_multi_hamiltonian_matrix(const class_model_finite & model, const class_edges_finite & edges);
            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_matrix(const class_model_finite & model, const class_edges_finite & edges);
            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_subspace_matrix(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs );
            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_subspace_matrix_new(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs );
            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_subspace_matrix_new(const class_model_finite & model, const class_edges_finite & edges, const std::vector<candidate_tensor> & candidate_list);
        }

        template <typename T>
        extern Eigen::Tensor<T,6> get_multi_hamiltonian_tensor(const class_model_finite & model, const class_edges_finite & edges);
        template <typename T>
        extern Eigen::Tensor<T,6> get_multi_hamiltonian_squared_tensor(const class_model_finite & model, const class_edges_finite & edges);
        template <typename T>
        extern Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_matrix(const class_model_finite & model, const class_edges_finite & edges);
        template <typename T>
        extern Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_squared_matrix(const class_model_finite & model, const class_edges_finite & edges);
        template <typename T>
        extern Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_squared_subspace_matrix(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs);
        template <typename T>
        extern Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_squared_subspace_matrix_new(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs);
        template <typename T>
        extern Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_squared_subspace_matrix_new(const class_model_finite & model, const class_edges_finite & edges, const std::vector<candidate_tensor> & candidate_list);


        inline ceres::GradientProblemSolver::Options ceres_default_options;

        inline bool no_state_in_window = false;

        extern std::vector<int> generate_size_list(int shape);


        template <typename T>  int sgn(const T val) {return (T(0) < val) - (val < T(0)); }
        extern double windowed_func_abs(double x,double window);
        extern double windowed_grad_abs(double x,double window);
        extern double windowed_func_pow(double x,double window);
        extern double windowed_grad_pow(double x,double window);
        extern std::pair<double,double> windowed_func_grad(double x,double window);



        namespace reports{
            using bfgs_tuple       = std::tuple<std::string,long,int,double,std::complex<double>,double,double,int,int,double>;
            using time_tuple       = std::tuple<double,double,double,double,double>;
            using eigs_tuple       = std::tuple<int,double,double,double,double,double,double>;
            inline std::vector<bfgs_tuple> bfgs_log;
            inline std::vector<time_tuple> time_log;
            inline std::vector<eigs_tuple> eigs_log;

            void print_bfgs_report();
            void print_time_report();
            void print_eigs_report();

            void bfgs_add_entry(const std::string & algorithm, long size, int rank, double energy, std::complex<double> variance,
                                double overlap, double norm, int iter, int counter, double time);
            void time_add_entry(double vH2v, double vHv, double vH2, double vH, double time);
            void eigs_add_entry(double nev, double max_olap, double min_olap, double eps, double eig_time,double ham_time, double lu_time);

        }

        void reset_timers();

        class ceres_base_functor : public ceres::FirstOrderFunction{
        public:
        protected:
            mutable double variance;
            mutable double energy  ;
            mutable double energy_reduced;
            mutable double energy_lower_bound;
            mutable double energy_upper_bound;
            mutable double energy_target;
            mutable double energy_min;
            mutable double energy_max;
            mutable double energy_dens;
            mutable double energy_target_dens;
            mutable double energy_window;
            mutable double energy_offset;
            mutable double norm_offset;
            mutable double norm;
            mutable size_t counter = 0;
            size_t length;
            size_t iteration;
            int    num_parameters;
            bool   have_bounds_on_energy = false;
            OMP omp;

        public:
            explicit ceres_base_functor(const class_tensors_finite & tensors, const class_algorithm_status &status);

            double get_variance   () const;
            double get_energy     () const;
            size_t get_count      () const;
            double get_norm       () const;
            int    NumParameters  () const final;
        };


    }

