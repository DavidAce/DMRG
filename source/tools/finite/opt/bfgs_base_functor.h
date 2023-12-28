
#pragma once
#include <array>
#include <ceres/first_order_function.h>
#include <ceres/manifold.h>
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

    class bfgs_base_functor : public ceres::FirstOrderFunction {
        protected:
        int                 num_parameters; // Includes lagrangian multiplier(s)
        long                size;           // Does not include lagrange multiplier (i.e. size + 1 = num_parameters if LagrangeNorm::ON)
        std::array<long, 3> dims;
        mutable double      fval{};
        mutable double      energy{};
        mutable double      variance{};
        mutable double      energy_shift;
        mutable double      energy_offset{};
        mutable double      delta_f{};
        mutable double      max_grad_norm = 0;
        mutable double      norm_offset{};
        mutable double      norm{};
        mutable double      resnorm = 0;
        mutable size_t      counter = 0;
        mutable long        ops     = 0;
        size_t              length;
        size_t              iteration;
        bool                have_bounds_on_energy = false;

        public:
        //        template<typename> friend class NormParametrization;
        explicit bfgs_base_functor(const TensorsFinite &tensors, const AlgorithmStatus &status);

        double              get_fval() const;
        double              get_energy() const;
        double              get_variance() const;
        size_t              get_count() const;
        double              get_norm() const;
        double              get_norm_offset() const;
        double              get_resnorm() const;
        double              get_delta_f() const;
        double              get_max_grad_norm() const;
        long                get_ops() const;
        long                get_size() const;
        std::array<long, 3> get_dims() const;
        int                 NumParameters() const final;

        void set_delta_f(double delta_f_) const;
        void set_max_grad_norm(double max_grad_norm_) const;

        std::unique_ptr<tid::ur> t_step;
        std::unique_ptr<tid::ur> t_H2n;
        std::unique_ptr<tid::ur> t_H2r;
        std::unique_ptr<tid::ur> t_nH2n;
        std::unique_ptr<tid::ur> t_Hn;
        std::unique_ptr<tid::ur> t_nHn;
    };
    /* clang-format on */

    template<typename FunctorType>
    class NormParametrization : public ceres::Manifold {
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
        bool PlusJacobian(const double *x, double *jacobian) const final;

        int AmbientSize() const final; // Size of x.
        int TangentSize() const final; // Size of delta.
    };

}