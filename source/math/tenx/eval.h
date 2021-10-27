#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

namespace tenx {
    template<typename Derived, typename Device = Eigen::DefaultDevice>
    class selfCleaningEvaluator {
        private:
        using Evaluator = Eigen::TensorEvaluator<const Eigen::TensorForcedEvalOp<const Derived>, Device>;
        Evaluator m_eval;

        public:
        using Scalar                        = typename Eigen::internal::remove_const<typename decltype(m_eval)::Scalar>::type;
        using DimType                       = typename decltype(m_eval)::Dimensions::Base;
        static constexpr auto NumDimensions = DimType{}.size();

        selfCleaningEvaluator(const Evaluator &eval) : m_eval(eval) {}
        template<int AccessLevel>
        selfCleaningEvaluator(const Eigen::TensorBase<Derived, AccessLevel> &expr, const Device &device = Device()) : m_eval(Evaluator(expr.eval(), device)) {
            m_eval.evalSubExprsIfNeeded(nullptr);
        }

        ~selfCleaningEvaluator() {
            // This whole point of this object is to call cleanup automatically on destruct.
            // If there are pending operations to evaluate, m_eval will allocate a buffer to hold a result,
            // which needs to be deallocated.
            m_eval.cleanup();
        }

        constexpr auto rank() { return DimType{}.size(); }
        constexpr auto dimensions() { return m_eval.dimensions(); }
        constexpr auto size() { return m_eval.size(); }
        constexpr auto data() { return m_eval.data(); }

        const Evaluator *operator->() const { return &m_eval; }
        Evaluator       *operator->() { return &m_eval; }
        constexpr auto   map() { return Eigen::TensorMap<Eigen::Tensor<Scalar, NumDimensions>>(m_eval.data(), m_eval.dimensions()); }
    };

    // Evaluates expressions if needed
    template<typename T, int AccessLevel, typename Device = Eigen::DefaultDevice>
    auto asEval(const Eigen::TensorBase<T, AccessLevel> &expr, const Device &device = Device()) {
        return selfCleaningEvaluator(expr, device);
    }
}