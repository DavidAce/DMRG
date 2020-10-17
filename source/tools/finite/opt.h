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
#include <tools/common/log.h>

class class_state_finite;
class class_model_finite;
class class_edges_finite;
class class_tensors_finite;
class class_algorithm_status;
class class_tic_toc;

namespace tools::finite::opt {
    class opt_state;

    using cplx = std::complex<double>;
    using real = double;

    template<typename T>
    using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    template<typename T, auto rank = 3>
    using TensorType = Eigen::Tensor<T, rank>;

    using Scalar = std::complex<double>;
    extern opt_state find_excited_state(const class_tensors_finite &tensors, const opt_state &initial_tensor, const class_algorithm_status &status, OptMode optMode, OptSpace optSpace,
                                         OptType optType);
    extern opt_state find_excited_state(const class_tensors_finite &tensors, const class_algorithm_status &status, OptMode optMode, OptSpace optSpace,
                                         OptType optType);
    extern Eigen::Tensor<Scalar, 3> find_ground_state(const class_tensors_finite &tensors, StateRitz ritz);
}

/* clang-format off */
namespace tools::finite::opt::internal{

    extern opt_state
           ceres_direct_optimization(const class_tensors_finite & tensors, const class_algorithm_status &status, OptType optType, OptMode optMode,OptSpace optSpace);

    extern opt_state
           ceres_direct_optimization(const class_tensors_finite & tensors, const opt_state &initial_tensor,
                                     const class_algorithm_status &status, OptType optType, OptMode optMode,OptSpace optSpace);

    extern opt_state
           ceres_subspace_optimization (const class_tensors_finite & tensors, const class_algorithm_status & status, OptType optType, OptMode optMode,OptSpace optSpace);

    extern opt_state
            ceres_subspace_optimization(const class_tensors_finite & tensors, const opt_state &initial_tensor,
                                      const class_algorithm_status &status, OptType optType, OptMode optMode,OptSpace optSpace);


    extern Eigen::Tensor<std::complex<double>,3> cppoptlib_optimization      (const class_tensors_finite & tensors, const class_algorithm_status & status);
    extern Eigen::Tensor<std::complex<double>,3> ground_state_optimization   (const class_tensors_finite & tensors, StateRitz ritz);
    extern Eigen::Tensor<std::complex<double>,3> ground_state_optimization   (const class_tensors_finite & tensors, std::string_view ritz);
    extern Eigen::Tensor<std::complex<double>,3> ham_sq_optimization         (const class_tensors_finite & tensors, OptType optType, OptMode optMode, OptSpace optSpace);
    extern Eigen::Tensor<std::complex<double>,3> ceres_rosenbrock_optimization (const class_state_finite & state);


    namespace subspace{
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
        extern std::vector<opt_state>
        find_candidates(const class_tensors_finite &tensors, double subspace_error_threshold, OptMode optMode, OptSpace optSpace);

        extern void
        filter_candidates(std::vector<opt_state> & candidate_list, double maximum_subspace_error, size_t max_accept);


        extern std::optional<size_t> get_idx_to_candidate_with_highest_overlap(const std::vector<opt_state> & candidate_list, double energy_llim_per_site, double energy_ulim_per_site);
        extern std::vector<size_t> get_idx_to_candidates_with_highest_overlap(const std::vector<opt_state> & candidate_list, size_t max_candidates, double energy_llim_per_site, double energy_ulim_per_site);

        extern Eigen::MatrixXcd get_eigvecs(const std::vector<opt_state> & candidate_list);
        extern Eigen::VectorXd get_eigvals(const std::vector<opt_state> & candidate_list);
        extern Eigen::VectorXd get_energies(const std::vector<opt_state> & candidate_list);
        extern Eigen::VectorXd get_energies_per_site(const std::vector<opt_state> & candidate_list);
        extern Eigen::VectorXd get_overlaps(const std::vector<opt_state> & candidate_list);
        extern double get_subspace_error(const std::vector<opt_state> & candidate_list);
        extern double get_subspace_error(const std::vector<double> &overlaps);
        extern Eigen::VectorXcd get_vector_in_subspace(const std::vector<opt_state> & candidate_list, size_t idx);
        extern Eigen::VectorXcd get_vector_in_subspace(const std::vector<opt_state> & candidate_list, const Eigen::VectorXcd & subspace_vector);
        extern Eigen::VectorXcd get_vector_in_fullspace(const std::vector<opt_state> & candidate_list, const Eigen::VectorXcd & subspace_vector);
    }


//        namespace local_hamiltonians{
//            extern Eigen::Tensor<std::complex<double>,6>   get_multisite_hamiltonian_tensor(const class_model_finite & model, const class_edges_finite & edges);
//            extern Eigen::Tensor<std::complex<double>,6>   get_multi_hamiltonian_squared_tensor(const class_model_finite & model, const class_edges_finite & edges);
//            extern Eigen::MatrixXcd                        get_multisite_hamiltonian_matrix(const class_model_finite & model, const class_edges_finite & edges);
//            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_matrix(const class_model_finite & model, const class_edges_finite & edges);
//            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_subspace_matrix(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs );
//            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_subspace_matrix_new(const class_model_finite & model, const class_edges_finite & edges, const Eigen::MatrixXcd & eigvecs );
//            extern Eigen::MatrixXcd                        get_multi_hamiltonian_squared_subspace_matrix_new(const class_model_finite & model, const class_edges_finite & edges, const std::vector<opt_tensor> & candidate_list);
//        }

//        template <typename T>
//        extern Eigen::Tensor<T,6> get_multi_hamiltonian_tensor(const class_model_finite & model, const class_edges_finite & edges);
//        template <typename T>
//        extern Eigen::Tensor<T,6> get_multi_hamiltonian_squared_tensor(const class_model_finite & model, const class_edges_finite & edges);
    template <typename T>
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multisite_hamiltonian_matrix(const class_model_finite & model, const class_edges_finite & edges);
    template <typename T>
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> get_multisite_hamiltonian_squared_subspace_matrix(const class_model_finite & model, const class_edges_finite & edges, const std::vector<opt_state> & candidate_list);


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

    namespace reports{
        struct bfgs_entry {
            std::string description;
            long size,space;
            double energy,variance,overlap,norm;
            size_t iter, counter;
            double time;
        };
        struct time_entry{
            double vH2v,vHv,vH2,vH,bfgs;
        };
        struct eigs_entry{
            long nev;
            double max_olap,min_olap,eps,eig_time,ham_time,lu_time;
        };

        inline std::vector<bfgs_entry> bfgs_log;
        inline std::vector<time_entry> time_log;
        inline std::vector<eigs_entry> eigs_log;

        void print_bfgs_report();
        void print_time_report();
        void print_eigs_report();

        void bfgs_add_entry(const std::string & description, long size,long space, double energy, double variance,
                            double overlap, double norm, size_t iter, size_t counter, double time);
        void bfgs_add_entry(const std::string & mode,const std::string & tag, const opt_state & tensor, std::optional<long> space = std::nullopt);
        void time_add_dir_entry();
        void time_add_sub_entry();
        void eigs_add_entry(long nev, double max_olap, double min_olap, double eps, double eig_time,double ham_time, double lu_time);

    }


    class ceres_base_functor : public ceres::FirstOrderFunction{
    protected:
        mutable double energy;
        mutable double energy_per_site;
        mutable double variance;
        mutable double variance_per_site;
        mutable double energy_reduced;
        mutable double energy_llim_per_site;
        mutable double energy_ulim_per_site;
        mutable double energy_tgt_per_site;
        mutable double energy_min_per_site;
        mutable double energy_max_per_site;
        mutable double energy_dens;
        mutable double energy_dens_target;
        mutable double energy_dens_window;
        mutable double energy_offset;
        mutable double norm_offset;
        mutable double norm;
        mutable size_t counter = 0;
        mutable long  ops = 0;
        size_t length;
        size_t iteration;
        int    num_parameters;
        bool   have_bounds_on_energy = false;

    public:
        explicit ceres_base_functor(const class_tensors_finite & tensors, const class_algorithm_status &status);

        double get_energy              () const;
        double get_energy_per_site     () const;
        double get_variance            () const;
        double get_variance_per_site   () const;
        size_t get_count               () const;
        double get_norm                () const;
        long   get_ops                 () const;
        int    NumParameters           () const final;

        std::unique_ptr<class_tic_toc> t_bfgs;
        std::unique_ptr<class_tic_toc> t_vH2;
        std::unique_ptr<class_tic_toc> t_vH2v;
        std::unique_ptr<class_tic_toc> t_vH;
        std::unique_ptr<class_tic_toc> t_vHv;
    };
    /* clang-format on */

    template<typename FunctorType>
    class CustomLogCallback : public ceres::IterationCallback {
        public:
        std::shared_ptr<spdlog::logger> log;
        const FunctorType &             functor;
        int    freq_log_iter = 5; //Wait at least this many iterations between logs
        int    last_log_iter = 0;
        int    init_log_iter = 0;
        double last_log_time = 0;
        double freq_log_time = 5; //Wait at least this many seconds between logs
        double init_log_time = 5;
        size_t last_count    = 0;

        explicit CustomLogCallback(const FunctorType &functor_);
//            : functor(functor_)
//        {
//            if(not log) log = tools::Logger::setLogger("xDMRG");
//            log->set_level(tools::log->level());
//            if(log->level() == spdlog::level::debug){
//                freq_log_iter = 10;
//                freq_log_time = 10;
//                init_log_time = 1;
//            }
//        }
        ceres::CallbackReturnType operator()(const ceres::IterationSummary &summary);
//        {
//            if(not log) return ceres::SOLVER_CONTINUE;
//            if(log->level() > spdlog::level::debug) return ceres::SOLVER_CONTINUE;
//            if(summary.iteration - last_log_iter  < freq_log_iter and summary.iteration > init_log_iter  ) return ceres::SOLVER_CONTINUE;
//            if(summary.cumulative_time_in_seconds - last_log_time < freq_log_time and summary.cumulative_time_in_seconds > init_log_time) return ceres::SOLVER_CONTINUE;
//            last_log_time = summary.cumulative_time_in_seconds;
//            last_log_iter = summary.iteration;
//            /* clang-format off */
//            log->debug("LBFGS: iter {:>5} f {:>8.5f} |Δf| {:>3.2e} "
//                      "|∇f|∞ {:>3.2e} |ΔΨ| {:3.2e} ls {:3.2e} evals {:>4}/{:<4} "
//                      "t_step {:<6} t_iter {:<6} t_tot {:<5} GOp/s {:<4.2f} | energy {:<18.15f} log₁₀var {:<6.6f}",
//                      summary.iteration,
//                      summary.cost,
//                      summary.cost_change,
//                      summary.gradient_max_norm,
//                      summary.step_norm, // By lbfgs
//                      summary.step_size, // By line search
//                      functor.get_count() - last_count,
//                      functor.get_count(),
//
//                      fmt::format("{:>.0f} ms",summary.step_solver_time_in_seconds * 1000),
//                      fmt::format("{:>.0f} ms",summary.iteration_time_in_seconds * 1000),
//                      fmt::format("{:>.3f} s",summary.cumulative_time_in_seconds),
//                      static_cast<double>(functor.get_ops()) / functor.t_vH2->get_last_interval()/1e9,
//                      functor.get_energy_per_site(),
//                      std::log10(functor.get_variance_per_site())
//            );
//            last_count = functor.get_count();
//            /* clang-format on */
//            return ceres::SOLVER_CONTINUE;
//        }
    };

}
