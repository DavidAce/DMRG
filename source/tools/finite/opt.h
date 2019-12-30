//
// Created by david on 2019-03-18.
//

#pragma once
// This trick avoids collision between the preprocessor
// symbol I and a template type size_t I in Eigen getting
// overridden and causing trouble at compile time.
#include <complex.h>
#undef I

// It is VERY IMPORTANT to set this flag when using CERES
// For some reason the Eigen destructor starts calling std::free instead of
// the handmade allocator free in calls inside CERES. I do not know the reason.
// Setting this flag here or globally helps.
#define EIGEN_MALLOC_ALREADY_ALIGNED 0
#include <tools/nmspc_tools.h>
#include <general/class_tic_toc.h>
#include <general/nmspc_omp.h>
#include <ceres/ceres.h>

class class_tic_toc;


namespace tools::finite::opt::internal{
        extern Eigen::Tensor<std::complex<double>,3> old_subspace_optimization(const class_state_finite &state,
                                                                               const class_simulation_status &sim_status,
                                                                               OptType optType, OptMode optMode);
        extern Eigen::Tensor<std::complex<double>,3> old_direct_optimization(const class_state_finite &state,
                                                                             const class_simulation_status &sim_status,
                                                                             OptType optType);
        extern Eigen::Tensor<std::complex<double>,3> ceres_direct_optimization(const class_state_finite &state,
                                                                               const class_simulation_status &sim_status,
                                                                               OptType optType);
        extern Eigen::Tensor<std::complex<double>,3> ceres_direct_optimization(const class_state_finite &state,
                                                                               const Eigen::Tensor<std::complex<double>,3> &theta,
                                                                               const class_simulation_status &sim_status,
                                                                               OptType optType);
        extern Eigen::Tensor<std::complex<double>,3> ceres_subspace_optimization   (const class_state_finite & state,
                                                                                    const class_simulation_status & sim_status,
                                                                                    OptType optType, OptMode optMode);
        extern Eigen::Tensor<std::complex<double>,3> cppoptlib_optimization      (const class_state_finite & state, const class_simulation_status & sim_status);
        extern Eigen::Tensor<std::complex<double>,4> ground_state_optimization   (const class_state_finite & state, std::string ritzstring = "SR");
        extern Eigen::Tensor<std::complex<double>,3> ceres_rosenbrock_optimization (const class_state_finite & state);

        namespace local_hamiltonians{
            extern Eigen::Tensor<std::complex<double>,6>   get_multi_hamiltonian_tensor(const class_state_finite & state);
            extern Eigen::Tensor<std::complex<double>,6>   get_multi_hamiltonian_squared_tensor(const class_state_finite & state);
            extern Eigen::MatrixXcd                        get_multi_hamiltonian_matrix(const class_state_finite & state);
            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_matrix(const class_state_finite & state);
            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_subspace_matrix(const class_state_finite & state, const Eigen::MatrixXcd & eigvecs );
            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_subspace_matrix_new(const class_state_finite & state, const Eigen::MatrixXcd & eigvecs );
        }

        template <typename T>
        Eigen::Tensor<T,6> get_multi_hamiltonian_tensor(const class_state_finite & state);
        template <typename T>
        Eigen::Tensor<T,6> get_multi_hamiltonian_squared_tensor(const class_state_finite & state);
        template <typename T>
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_matrix(const class_state_finite & state);
        template <typename T>
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_squared_matrix(const class_state_finite & state);
        template <typename T>
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_squared_subspace_matrix(const class_state_finite & state, const Eigen::MatrixXcd & eigvecs);
        template <typename T>
        Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_squared_subspace_matrix_new(const class_state_finite & state, const Eigen::MatrixXcd & eigvecs);


        inline ceres::GradientProblemSolver::Options ceres_default_options;

        inline bool no_state_in_window = false;
        inline double subspace_error_threshold = 1e-8;

        extern std::vector<int> generate_size_list(size_t shape);


        template <typename T>  int sgn(const T val) {return (T(0) < val) - (val < T(0)); }
        extern double windowed_func_abs(double x,double window);
        extern double windowed_grad_abs(double x,double window);
        extern double windowed_func_pow(double x,double window);
        extern double windowed_grad_pow(double x,double window);
        extern std::pair<double,double> windowed_func_grad(double x,double window);



        namespace reports{
            using direct_opt_tuple = std::tuple<std::string,int,double,std::complex<double>,double,double,int,int,double>;
            using subspc_opt_tuple = std::tuple<std::string,int,double,double,double,double,int,int,double>;
            using lbfgs_tuple      = std::tuple<double,double,double,double,double>;
            using eig_tuple        = std::tuple<int,double,double,double,double,double,double>;
//            std::vector<log_tuple> opt_log;
            void print_report(const std::vector<direct_opt_tuple> &opt_log);
            void print_report(const std::vector<subspc_opt_tuple> &opt_log);
            void print_report(const std::vector<eig_tuple> &eig_log);
            void print_report(const lbfgs_tuple lbfgs_log);
        }

        void reset_timers();
        inline std::unique_ptr<class_tic_toc> t_opt  =  std::make_unique<class_tic_toc>(true,5,"t_opt ");
        inline std::unique_ptr<class_tic_toc> t_eig  =  std::make_unique<class_tic_toc>(true,5,"t_eig ");
        inline std::unique_ptr<class_tic_toc> t_ham  =  std::make_unique<class_tic_toc>(true,5,"t_ham ");
        inline std::unique_ptr<class_tic_toc> t_tot  =  std::make_unique<class_tic_toc>(true,5,"t_tot ");
        inline std::unique_ptr<class_tic_toc> t_vH2v =  std::make_unique<class_tic_toc>(true,5,"t_vH2v");
        inline std::unique_ptr<class_tic_toc> t_vHv  =  std::make_unique<class_tic_toc>(true,5,"t_vHv ");
        inline std::unique_ptr<class_tic_toc> t_vH2  =  std::make_unique<class_tic_toc>(true,5,"t_vH2 ");
        inline std::unique_ptr<class_tic_toc> t_vH   =  std::make_unique<class_tic_toc>(true,5,"t_vH  ");
        inline std::unique_ptr<class_tic_toc> t_op   =  std::make_unique<class_tic_toc>(true,5,"t_op  ");

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
            mutable int    counter = 0;
            size_t length;
            int    iteration;
            int    num_parameters;
            bool   have_bounds_on_energy = false;
            OMP omp;

        public:
            explicit ceres_base_functor(const class_state_finite & state, const class_simulation_status &sim_status);

            double get_variance   () const;
            double get_energy     () const;
            size_t get_count      () const;
            double get_norm       () const;
            int    NumParameters  () const final;
        };


    }

