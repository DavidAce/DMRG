
#pragma once
#include <ceres/first_order_function.h>
#include <ceres/local_parameterization.h>
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
    enum class LagrangeNorm { ON, OFF };

    class lbfgs_base_functor : public ceres::FirstOrderFunction {
        protected:
        int                    num_parameters; // Includes lagrangian multiplier(s)
        long                   size;           // Does not include lagrange multiplier (i.e. size + 1 = num_parameters if LagrangeNorm::ON)
        Eigen::DSizes<long, 3> dims;
        mutable double         energy;
        mutable double         energy_per_site;
        mutable double         variance;
        mutable double         variance_per_site;
        mutable double         energy_shift;
        mutable double         energy_llim_per_site;
        mutable double         energy_ulim_per_site;
        mutable double         energy_tgt_per_site;
        mutable double         energy_min_per_site;
        mutable double         energy_max_per_site;
        mutable double         energy_dens;
        mutable double         energy_dens_target;
        mutable double         energy_dens_window;
        mutable double         energy_offset;
        mutable double         delta_f;
        mutable double         max_grad_norm = 0;
        mutable double         norm_offset;
        mutable double         norm;
        mutable size_t         counter = 0;
        mutable long           ops     = 0;
        size_t                 length;
        size_t                 iteration;
        bool                   have_bounds_on_energy = false;

        public:
        //        template<typename> friend class NormParametrization;
        explicit lbfgs_base_functor(const TensorsFinite &tensors, const AlgorithmStatus &status);

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
    class NormParametrization : public ceres::LocalParameterization {
        public:
        const FunctorType &functor;
        explicit NormParametrization(const FunctorType &functor_);

        // Generalization of the addition operation,
        //
        //   x_plus_delta = Plus(x, delta)
        //
        // with the condition that Plus(x, 0) = x.
        bool Plus(const double *x, const double *delta, double *x_plus_delta) const final;

        // The jacobian of Plus(x, delta) w.r.t delta at delta = 0.
        // jacobian is a row-major GlobalSize() x LocalSize() matrix.
        bool ComputeJacobian(const double *x, double *jacobian) const final;

        int GlobalSize() const final; // Size of x.
        int LocalSize() const final;  // Size of delta.
    };

}