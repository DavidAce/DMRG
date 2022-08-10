#include "bfgs_subspace_functor.h"
#include "debug/exceptions.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/opt/opt-internal.h"
using namespace tools::finite::opt::internal;

template<typename Scalar>
tools::finite::opt::internal::bfgs_subspace_functor<Scalar>::bfgs_subspace_functor(const TensorsFinite &tensors, const AlgorithmStatus &status,
                                                                                   const MatrixType &H2_subspace, const Eigen::VectorXd &eigvals_)
    : bfgs_base_functor(tensors, status), H2(H2_subspace), eigvals(eigvals_) {
    double nonhermiticity = (H2 - H2.adjoint()).cwiseAbs().sum() / static_cast<double>(H2.size());
    if(nonhermiticity > 1e-12) tools::log->error("bfgs_subspace_functor(): H2 is not Hermitian: {:.16f}", nonhermiticity);

    num_parameters = static_cast<int>(eigvals.size());

    // If the scalar is complex we optimize a problem of twice the linear size
    if constexpr(std::is_same_v<Scalar, cplx>) num_parameters *= 2;
}

template<typename Scalar>
bool tools::finite::opt::internal::bfgs_subspace_functor<Scalar>::Evaluate(const double *v_double_double, double *fx, double *grad_double_double) const {
    auto       t_step_token = t_step->tic_token();
    Scalar     nH2n, nHn, var;
    double     vv, log10var;
    VectorType Hn, H2n;
    int        vecSize = NumParameters();
    if constexpr(std::is_same_v<Scalar, cplx>) vecSize = NumParameters() / 2;
    Eigen::Map<const VectorType> v(reinterpret_cast<const Scalar *>(v_double_double), vecSize);

    VectorType n = v.normalized();
    vv           = v.squaredNorm();
    norm         = std::sqrt(vv);

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
    var        = nH2n - nHn * nHn;
    double eps = std::numeric_limits<double>::epsilon();
    if(std::real(var) < -eps or std::real(nH2n) < -eps)
        tools::log->trace("Counter = {} | SHINV | "
                          "negative: "
                          "var  {:.16f} + {:.16f}i | "
                          "nHn  {:.16f} + {:.16f}i | "
                          "nH2n {:.16f} + {:.16f}i",
                          counter, std::real(var), std::imag(var), std::real(nHn), std::imag(nHn), std::real(nH2n), std::imag(nH2n));

    //    var = std::abs(var);
    //    var = std::real(var) == 0.0 ? eps : var;

    energy      = std::real(nHn + energy_shift);
    variance    = std::abs(var);
    norm_offset = std::abs(vv) - 1.0;
    //    log10var          = std::log10(variance);
    // Here we work with a small offset on the variance:
    //      When the variance is very low, numerical noise will sometimes cause nH2n < 0 or var < 0.
    //      I think this is caused by H2 being slightly non-symmetric
    double var_offset = std::clamp(std::real(var) + 1e-12, eps, std::abs(var) + 1e-12);
    log10var          = std::log10(var_offset);
    if(fx != nullptr) { fx[0] = log10var; }
    auto pref = std::is_same_v<Scalar, real> ? 2.0 : 1.0; // Factor 2 for real

    Eigen::Map<VectorType> grad(reinterpret_cast<Scalar *>(grad_double_double), vecSize);
    if(grad_double_double != nullptr) {
        auto one_over_norm = 1.0 / norm;
        auto var_1         = 1.0 / var_offset / std::log(10);
        grad               = pref * one_over_norm * (H2n - 2.0 * nHn * Hn - (nH2n - 2.0 * nHn * nHn) * n);
        max_grad_norm      = grad.template lpNorm<Eigen::Infinity>(); // To monitor the actual gradient norm of the optimization (not its logarithm)
        grad *= var_1;                                                // Because we are optimizing the logarithm.
    }
    fval = fx[0];

    if(std::isnan(log10var) or std::isinf(log10var)) {
        tools::log->warn("σ²H is invalid");
        tools::log->warn("σ²H             = {:8.2e}", variance);
        tools::log->warn("mv          = {}", counter);
        tools::log->warn("vecsize         = {}", vecSize);
        tools::log->warn("vv              = {:.16f} + i{:.16f}", std::real(vv), std::imag(vv));
        tools::log->warn("nH2n            = {:.16f} + i{:.16f}", std::real(nH2n), std::imag(nH2n));
        tools::log->warn("nHn             = {:.16f} + i{:.16f}", std::real(nHn), std::imag(nHn));
        tools::log->warn("var             = {:.16f} + i{:.16f}", std::real(var), std::imag(var));
        tools::log->warn("energy offset   = {:.16f}", energy_offset);
        tools::log->warn("energy shift    = {:.16f}", energy_shift);
        tools::log->warn("norm            = {:.16f}", norm);
        tools::log->warn("norm   offset   = {:.16f}", norm_offset);
        throw except::runtime_error("Direct functor failed at mv = {}", counter);
    }
    counter++;
    return true;
}

template class tools::finite::opt::internal::bfgs_subspace_functor<real>;
template class tools::finite::opt::internal::bfgs_subspace_functor<cplx>;
