#include "bfgs_simps_functor.h"
#include "config/debug.h"
#include "opt-internal.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"

namespace debug {
    template<typename Derived>
    bool hasNaN(const Eigen::EigenBase<Derived> &obj, [[maybe_unused]] std::string_view name = "") {
        return obj.derived().hasNaN();
    }

    template<typename Scalar, auto rank>
    bool hasNaN(const Eigen::Tensor<Scalar, rank> &tensor, std::string_view name = "") {
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> vector(tensor.data(), tensor.size());
        return hasNaN(vector, name);
    }
    template<typename Scalar, auto rank>
    bool hasNaN(const Eigen::TensorMap<const Eigen::Tensor<Scalar, rank>> &tensor, std::string_view name = "") {
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> vector(tensor.data(), tensor.size());
        return hasNaN(vector, name);
    }
}

using namespace tools::finite::opt::internal;

template<typename Scalar>
bfgs_simps_functor<Scalar>::bfgs_simps_functor(const TensorMap3 &mps_,  /*!< The current mps  */
                                               const TensorMap3 &envL_, /*!< The left block tensor for H.  */
                                               const TensorMap3 &envR_, /*!< The right block tensor for H.  */
                                               const TensorMap4 &mpo_   /*!< The H MPO's  */
                                               )
    : mps(mps_), envL(envL_), envR(envR_), mpo(mpo_)

{
    tools::log->trace("Constructing direct functor");
    dims           = mps.dimensions();
    size           = mps.size();
    num_parameters = static_cast<int>(size); // number of doubles to optimize (i.e. 2x if complex)
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) { num_parameters *= 2; }
    if constexpr(settings::debug) {
        std::string msg;
        if(debug::hasNaN(mps)) msg.append("\t mps\n");
        if(debug::hasNaN(mpo)) msg.append("\t mpo\n");
        if(debug::hasNaN(envL)) msg.append("\t envL\n");
        if(debug::hasNaN(envR)) msg.append("\t envR\n");
        if(not msg.empty()) throw std::runtime_error(fmt::format("The following objects have nan's:\n{}", msg));
    }

    tools::log->trace("- Allocating memory for matrix-vector products");
    Hm_tensor.resize(dims);
    Hv_tensor.resize(dims);
    auto m = Eigen::Map<const VectorType>(mps.data(), mps.size()); // |Ψ>
    compute_Hm(m);
}

template<typename Scalar>
int bfgs_simps_functor<Scalar>::NumParameters() const {
    return num_parameters;
}

template<typename Scalar>
bool bfgs_simps_functor<Scalar>::Evaluate(const double *v_double_double, double *fx, double *grad_double_double) const {
    Eigen::Map<const VectorType> v(reinterpret_cast<const Scalar *>(v_double_double), size); // |φ>
    Eigen::Map<const VectorType> m(mps.data(), mps.size());                                  // |Ψ>
    if constexpr(settings::debug) {
        if(v.hasNaN()) throw std::runtime_error(fmt::format("bfgs_simps_functor::Evaluate: v has nan's\n{}", v));
    }

    compute_Hv(v);
    auto   Hv = Eigen::Map<VectorType>(Hv_tensor.data(), Hv_tensor.size());
    double f  = (Hv - m).squaredNorm();
    //    if(fx != nullptr) fx[0] = f;
    if(fx != nullptr) fx[0] = std::log(f);

    Eigen::Map<VectorType> grad(reinterpret_cast<Scalar *>(grad_double_double), size);
    if(grad_double_double != nullptr) {
        compute_Hv(Hv);                                                        // TODO: YEEAH! WE USE H * H,  Not H²!! So just apply H on n again.
        auto HHv = Eigen::Map<VectorType>(Hv_tensor.data(), Hv_tensor.size()); // This is (H-E) * (H-E) |n>, not (H-E)²|n> as in variance.
        auto Hm  = Eigen::Map<VectorType>(Hm_tensor.data(), Hm_tensor.size()); // Hm does not change between iterations
        grad     = (HHv - Hm) / f;
        //        grad     = (HHv - Hm);
    }
    return true;
}

template<typename Scalar>
void bfgs_simps_functor<Scalar>::compute_Hv(const VectorType &v) const {
    // Computes  Hv_tensor = H|v>
    auto v_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(v.derived().data(), dims);
    tools::common::contraction::matrix_vector_product(Hv_tensor, v_tensor, mpo, envL, envR);
}
template<typename Scalar>
void bfgs_simps_functor<Scalar>::compute_Hm(const VectorType &m) const {
    // Computes  Hm_tensor = H|m>
    auto m_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(m.derived().data(), dims);
    tools::common::contraction::matrix_vector_product(Hm_tensor, m_tensor, mpo, envL, envR);
}

template class tools::finite::opt::internal::bfgs_simps_functor<real>;
template class tools::finite::opt::internal::bfgs_simps_functor<cplx>;
