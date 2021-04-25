//
// Created by david on 2019-07-15.
//
#include "ceres_subspace_functor.h"
#include <tensors/class_tensors_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt-internal/opt-internal.h>
using namespace tools::finite::opt::internal;

template<typename Scalar>
tools::finite::opt::internal::ceres_subspace_functor<Scalar>::ceres_subspace_functor(const class_tensors_finite &tensors, const class_algorithm_status &status,
                                                                                     const MatrixType &H2_subspace, const Eigen::VectorXd &eigvals_)
    : ceres_base_functor(tensors, status), H2(H2_subspace), eigvals(eigvals_) {
    double nonhermiticity = (H2 - H2.adjoint()).cwiseAbs().sum() / static_cast<double>(H2.size());
    if(nonhermiticity > 1e-12)
        tools::log->error("ceres_subspace_functor(): H2 is not Hermitian: {:.16f}", nonhermiticity);

    num_parameters = static_cast<int>(eigvals.size());

    // If the scalar is complex we optimize a problem of twice the linear size
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) num_parameters *= 2;
}

template<typename Scalar>
bool tools::finite::opt::internal::ceres_subspace_functor<Scalar>::Evaluate(const double *v_double_double, double *fx, double *grad_double_double) const {
    auto t_step_token  = t_step->tic_token();
    Scalar     nH2n, nHn, var;
    double     vv, log10var;
    VectorType Hn, H2n;
    int        vecSize = NumParameters();
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) { vecSize = NumParameters() / 2; }
    Eigen::Map<const VectorType> v(reinterpret_cast<const Scalar *>(v_double_double), vecSize);

    VectorType n = v.normalized();
    vv   = v.squaredNorm();
    norm = std::sqrt(vv);

    t_Hn->tic();
    Hn = eigvals.asDiagonal() * n;
    t_Hn->toc();

    t_nHn->tic();
    nHn = n.dot(Hn);
    t_nHn->toc();

    t_H2n->tic();
    H2n = H2.template selfadjointView<Eigen::Lower>() * n;
    t_H2n->toc();

    t_nH2n->tic();
    nH2n = n.dot(H2n);
    t_nH2n->toc();

    // Do this next bit carefully to avoid negative variance when numbers are very small
    var = nH2n - nHn * nHn;
    double eps = std::numeric_limits<double>::epsilon();
    if(std::real(var) < -eps or std::real(nH2n) < -eps)
        tools::log->debug("Counter = {} | SUBSPACE | "
                          "negative: "
                          "var  {:.16f} + {:.16f}i | "
                          "nHn  {:.16f} + {:.16f}i | "
                          "nH2n {:.16f} + {:.16f}i",
                          counter,
                          std::real(var), std::imag(var),
                          std::real(nHn), std::imag(nHn),
                          std::real(nH2n), std::imag(nH2n));

//    var = std::abs(var);
//    var = std::real(var) == 0.0 ? eps : var;

    energy            = std::real(nHn + energy_reduced);
    energy_per_site   = energy / static_cast<double>(length);
    variance          = std::abs(var);
    variance_per_site = variance / static_cast<double>(length);
    norm_offset       = std::abs(vv) - 1.0;
//    log10var          = std::log10(variance);
    // Here we work with a small offset on the variance:
    //      When the variance is very low, numerical noise will sometimes cause nH2n < 0 or var < 0.
    //      I think this is caused by H2 being slightly non-symmetric
    double var_offset = std::clamp(std::real(var) + 1e-12, eps, std::abs(var) + 1e-12 )  ;
    log10var          = std::log10(var_offset);
    if(fx != nullptr) { fx[0] = log10var; }

    Eigen::Map<VectorType> grad(reinterpret_cast<Scalar *>(grad_double_double), vecSize);
    if(grad_double_double != nullptr) {
        auto one_over_norm  = 1.0/norm;
        auto var_1 = 1.0 / var_offset / std::log(10);
        grad       = var_1 * one_over_norm * (H2n - 2.0 * nHn * Hn - (nH2n - 2.0 * nHn * nHn) * n);
        if constexpr(std::is_same<Scalar, double>::value) { grad *= 2.0; }
        grad_max_norm = grad.template lpNorm<Eigen::Infinity>();
    }

    if(std::isnan(log10var) or std::isinf(log10var)) {
        tools::log->warn("log₁₀ variance is invalid");
        tools::log->warn("log₁₀(var)      = {:.16f}", std::log10(variance));
        tools::log->warn("counter         = {}", counter);
        tools::log->warn("vecsize         = {}", vecSize);
        tools::log->warn("vv              = {:.16f} + i{:.16f}", std::real(vv), std::imag(vv));
        tools::log->warn("nH2n            = {:.16f} + i{:.16f}", std::real(nH2n), std::imag(nH2n));
        tools::log->warn("nHn             = {:.16f} + i{:.16f}", std::real(nHn), std::imag(nHn));
        tools::log->warn("var             = {:.16f} + i{:.16f}", std::real(var), std::imag(var));
        tools::log->warn("energy offset   = {:.16f}", energy_offset);
        tools::log->warn("energy reduced  = {:.16f}", energy_reduced);
        tools::log->warn("norm            = {:.16f}", norm);
        tools::log->warn("norm   offset   = {:.16f}", norm_offset);
        throw std::runtime_error("Direct functor failed at counter = " + std::to_string(counter));
    }

    counter++;
    return true;
}

template class tools::finite::opt::internal::ceres_subspace_functor<double>;
template class tools::finite::opt::internal::ceres_subspace_functor<std::complex<double>>;
