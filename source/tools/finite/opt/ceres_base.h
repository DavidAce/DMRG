
#pragma once
#include <ceres/first_order_function.h>
#include <ceres/iteration_callback.h>
#include <memory>
class StateFinite;
class ModelFinite;
class EdgesFinite;
class TensorsFinite;
class AlgorithmStatus;

namespace tid {
    class ur;
}

namespace spdlog {
    class logger;
}
namespace tools::finite::opt::internal {

    class ceres_base_functor : public ceres::FirstOrderFunction {
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
        mutable double delta_f;
        mutable double max_grad_norm = 0;
        mutable double norm_offset;
        mutable double norm;
        mutable size_t counter = 0;
        mutable long   ops     = 0;
        size_t         length;
        size_t         iteration;
        int            num_parameters;
        bool           have_bounds_on_energy = false;

        public:
        explicit ceres_base_functor(const TensorsFinite &tensors, const AlgorithmStatus &status);

        double get_energy() const;
        double get_energy_per_site() const;
        double get_variance() const;
        double get_variance_per_site() const;
        size_t get_count() const;
        double get_norm() const;
        double get_norm_offset() const;
        double get_delta_f() const;
        double get_max_grad_norm() const;
        long   get_ops() const;
        int    NumParameters() const final;

        void set_delta_f(double delta_f_) const;
        void set_max_grad_norm(double max_grad_norm_) const;

        std::unique_ptr<tid::ur> t_step;
        std::unique_ptr<tid::ur> t_H2n;
        std::unique_ptr<tid::ur> t_nH2n;
        std::unique_ptr<tid::ur> t_Hn;
        std::unique_ptr<tid::ur> t_nHn;
    };
    /* clang-format on */

    template<typename FunctorType>
    class CustomLogCallback : public ceres::IterationCallback {
        public:
        std::shared_ptr<spdlog::logger> log;
        const FunctorType              &functor;
        // The default values below are modified during construction when settings::debug == true
        int    freq_log_iter = 10; // Permit logs after this many iterations
        int    last_log_iter = 0;
        int    init_log_iter = 0; // Wait at least this long before first log
        double freq_log_time = 5; // Permit logs after this many seconds
        double last_log_time = 0;
        double init_log_time = 0; // Wait at least this long before first log
        size_t last_count    = 0;

        explicit CustomLogCallback(const FunctorType &functor_);
        ceres::CallbackReturnType operator()(const ceres::IterationSummary &summary);
    };

}