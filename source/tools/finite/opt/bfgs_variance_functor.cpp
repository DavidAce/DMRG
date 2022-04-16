#include "bfgs_variance_functor.h"
#include "config/debug.h"
#include "math/num.h"
#include "math/svd.h"
#include "opt-internal.h"
#include "tensors/TensorsFinite.h"
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
}

using namespace tools::finite::opt::internal;

template<typename Scalar, LagrangeNorm lagrangeNorm>
bfgs_variance_functor<Scalar, lagrangeNorm>::bfgs_variance_functor(const TensorsFinite &tensors, const AlgorithmStatus &status)
    : bfgs_base_functor(tensors, status) {
    tools::log->trace("Constructing direct functor");

    if constexpr(std::is_same_v<Scalar, real>) {
        tools::log->trace("- Generating real-valued multisite components");
        mpo                 = tensors.get_multisite_mpo().real();
        mpo2                = tensors.get_multisite_mpo_squared().real();
        const auto &env_ene = tensors.get_multisite_env_ene_blk();
        const auto &env_var = tensors.get_multisite_env_var_blk();
        envL                = env_ene.L.real();
        envR                = env_ene.R.real();
        env2L               = env_var.L.real();
        env2R               = env_var.R.real();
    } else if constexpr(std::is_same_v<Scalar, cplx>) {
        tools::log->trace("- Generating complex-valued multisite components");
        mpo                 = tensors.get_multisite_mpo();
        mpo2                = tensors.get_multisite_mpo_squared();
        const auto &env_ene = tensors.get_multisite_env_ene_blk();
        const auto &env_var = tensors.get_multisite_env_var_blk();
        envL                = env_ene.L;
        envR                = env_ene.R;
        env2L               = env_var.L;
        env2R               = env_var.R;
    }
    if constexpr(settings::debug) {
        std::string msg;
        if(debug::hasNaN(mpo)) msg.append("\t mpo\n");
        if(debug::hasNaN(mpo2)) msg.append("\t mpo2\n");
        if(debug::hasNaN(envL)) msg.append("\t envL\n");
        if(debug::hasNaN(envR)) msg.append("\t envR\n");
        if(debug::hasNaN(env2L)) msg.append("\t env2L\n");
        if(debug::hasNaN(env2R)) msg.append("\t env2R\n");
        if(not msg.empty()) throw std::runtime_error(fmt::format("The following objects have nan's:\n{}", msg));
    }

    tools::log->trace("- Allocating memory for matrix-vector products");
    Hn_tensor.resize(dims);
    H2n_tensor.resize(dims);
    if constexpr(lagrangeNorm == LagrangeNorm::ON) {
        H2r_tensor.resize(dims);
        residual.resize(size);
    }
    // Reminder
    // size below is the number of scalars to optimize excluding lagrange multipliers
    // num_parameters is the number of doubles to optimize (i.e. 2x if complex) including lagrange multipliers
    num_parameters = static_cast<int>(size);
    if constexpr(lagrangeNorm == LagrangeNorm::ON) num_parameters += 2; // Include lagrange multiplier(s) here
    if constexpr(std::is_same_v<Scalar, cplx>) { num_parameters *= 2; }
}

template<typename Scalar, LagrangeNorm lagrangeNorm>
bool bfgs_variance_functor<Scalar, lagrangeNorm>::Evaluate(const double *v_double_double, double *fx, double *grad_double_double) const {
    t_step->tic();
    Scalar                       var;
    Scalar                       nHn, nH2n;
    double                       vv, log10var;
    Eigen::Map<const VectorType> v(reinterpret_cast<const Scalar *>(v_double_double), size); // v does not include the lagrange multiplier

    if constexpr(settings::debug) {
        if(v.hasNaN()) throw std::runtime_error(fmt::format("bfgs_variance_functor::Evaluate: v has nan's at mv {}\n{}", counter, v));
    }

    VectorType n_temp;
    auto       n_data = v.data();
    // When not using LagrangeNorm we enforce the normalization here
    if constexpr(lagrangeNorm == LagrangeNorm::OFF) {
        n_temp = v.normalized();
        n_data = n_temp.data();
    }
    Eigen::Map<const VectorType> n(n_data, v.size());

    vv   = v.squaredNorm();
    norm = std::sqrt(vv);

    get_H2n(n);
    get_Hn(n);

    auto Hn  = Eigen::Map<VectorType>(Hn_tensor.data(), Hn_tensor.size());
    auto H2n = Eigen::Map<VectorType>(H2n_tensor.data(), H2n_tensor.size());

    print_path = false;

    t_nHn->tic();
    nHn = n.dot(Hn);
    t_nHn->toc();
    t_nH2n->tic();
    nH2n = n.dot(H2n);
    t_nH2n->toc();

    // Do this next bit carefully to avoid negative variance when numbers are very small
    var = nH2n - nHn * nHn;
    //    var        = nH2n;
    double eps = std::numeric_limits<double>::epsilon();
    if((std::real(var) < -eps or std::real(nH2n) < -eps))
        tools::log->trace("Counter = {} | BFGS | "
                          "negative: "
                          "var  {:.16f} + {:.16f}i | "
                          "nHn  {:.16f} + {:.16f}i | "
                          "nH2n {:.16f} + {:.16f}i | "
                          "grad {:8.4e}",
                          counter, std::real(var), std::imag(var), std::real(nHn), std::imag(nHn), std::real(nH2n), std::imag(nH2n), max_grad_norm);

    //    var = std::abs(var);
    //    var = std::real(var) == 0.0 ? eps : var;

    energy      = std::real(nHn + energy_shift);
    variance    = std::abs(var);
    norm_offset = std::abs(vv - 1.0);
    //    log10var          = std::log10(variance);
    // Here we work with a small offset on the variance:
    //      When the variance is very low, numerical noise will sometimes cause nH2n < 0 or var < 0.
    //      I think this is caused by H2 being slightly non-symmetric
    //    double var_offset = std::clamp(std::real(var) + 1e-12, eps, std::abs(var) + 1e-12);
    log10var  = std::log10(std::abs(var));
    auto pref = std::is_same_v<Scalar, real> ? 2.0 : 1.0; // Factor 2 for real

    Eigen::Map<VectorType> grad(reinterpret_cast<Scalar *>(grad_double_double), size);
    if(grad_double_double != nullptr) {
        auto one_over_norm = 1.0 / norm;
        auto var_1         = 1.0 / std::real(var) / std::log(10);
        grad               = pref * one_over_norm * (H2n - 2.0 * nHn * Hn - (nH2n - 2.0 * nHn * nHn) * n);
        max_grad_norm      = grad.template lpNorm<Eigen::Infinity>(); // To monitor the actual gradient norm of the optimization (not its logarithm)
        grad *= var_1;                                                // Because we are optimizing the logarithm.
    }
    if(fx != nullptr) { fx[0] = log10var; }

    if constexpr(lagrangeNorm == LagrangeNorm::ON) {
        // Here we define the norm constraint by using the lagrange multiplier trick:
        //      f(x)
        // is replaced with
        //      L(x,lambda) = f(x) + lambda1 * g(x) + lambda2 * h(x)
        // where g(x) = | <x|x> - 1 |
        // where h(x) = |r| = | H²x - E²x | (the 2-norm of the residual on the eigenvalue eq H²x = E²x )
        // Note that dh(x)/dx = 1/2sqrt(f) * (H²r-E²r)
        residual = H2n - nH2n * n; // aka r
        resnorm  = residual.norm();

        double g = std::abs(vv - 1.0); // aka g(x)
        double h = resnorm;            // aka h(x)

        if(fx != nullptr) fx[0] += g + std::log10(h);
        if(grad_double_double != nullptr) {
            get_H2r(residual);
            auto H2r               = Eigen::Map<VectorType>(H2r_tensor.data(), H2r_tensor.size());
            auto grad_w_multiplier = Eigen::Map<VectorType>(reinterpret_cast<Scalar *>(grad_double_double), size + 2);
            grad_w_multiplier.topRows(size) += pref * num::sign(vv - 1.0) * v; // aka  += lambda * dg(x)/dx = lambda * sign(x) * x
            grad_w_multiplier.topRows(size) += (H2r - nH2n * residual) / (2 * std::pow(resnorm, 1.5) * std::log(10)); // aka dh(x)/dx = (H²r-E²r) / 2sqrt(f)
            grad_w_multiplier.bottomRows(2)[0] = g;
            grad_w_multiplier.bottomRows(2)[1] = std::log10(h);
        }
    }

    if(std::isnan(log10var) or std::isinf(log10var)) {
        tools::log->warn("σ²H is invalid");
        tools::log->warn("σ²H             = {:8.2e}", variance);
        tools::log->warn("mv         = {}", counter);
        tools::log->warn("size            = {}", size);
        tools::log->warn("vv              = {:.16f} + i{:.16f}", std::real(vv), std::imag(vv));
        tools::log->warn("nH2n            = {:.16f} + i{:.16f}", std::real(nH2n), std::imag(nH2n));
        tools::log->warn("nHn             = {:.16f} + i{:.16f}", std::real(nHn), std::imag(nHn));
        tools::log->warn("var             = {:.16f} + i{:.16f}", std::real(var), std::imag(var));
        tools::log->warn("energy offset   = {:.16f}", energy_offset);
        tools::log->warn("energy shift    = {:.16f}", energy_shift);
        tools::log->warn("norm            = {:.16f}", norm);
        tools::log->warn("norm   offset   = {:.16f}", norm_offset);
        throw std::runtime_error("Direct functor failed at mv = " + std::to_string(counter));
    }

    counter++;
    t_step->toc();
    return true;
}

template<typename Scalar, LagrangeNorm lagrangeNorm>
void bfgs_variance_functor<Scalar, lagrangeNorm>::get_H2n(const VectorType &v) const {
    t_H2n->tic();
    auto v_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(v.derived().data(), dims);
    tools::common::contraction::matrix_vector_product(H2n_tensor, v_tensor, mpo2, env2L, env2R);
    t_H2n->toc();

    ops = tools::finite::opt::internal::get_ops(dims[0], dims[1], dims[2], mpo.dimension(0));

    /*
     * NOTE 2020-10-05
     * I ran some benchmarks to figure out if the order of indices matter, for instance
     * if we should list "spin" dimensions before mpo dimension or vice versa.
     * It turns out there is a large difference in favor of taking spin dimensions before
     * mpo dimensions. The number of operations increases with the number of sites.
     * The percentages below are 100*(ops_m - ops_d) / ops_d, where _d and _m denote
     * taking spin dim "d" first or mpo dim "m" first.
     *
     *          chi=256          chi=512
     * l = 2:    -0.65%           -0.32%
     * l = 3:     3.68%            1.84%
     * l = 4:    25.17%           12.97%
     * l = 5:   110.31%           59.80%
     * l = 6:   399.50% (*)      234.72%
     * l = 7:  1243.46%          815.45%
     * l = 8:  3370.13%         2498.57%
     *
     * However, it turns out that Eigen already switches the order around to take the fastest route
     * So it does not matter for us in which order we do it... i.e. a line like
     *      .contract(mpo, tenx::idx({0, 3}, {2, 0}))
     * is equivalent to
     *      .contract(mpo, tenx::idx({3, 0}, {0, 2}))
     *
     * The following benchmark result revealed this fact:
     *
     * [2020-10-05 15:21:26][tensorbench][  info  ] H²|Ψ> version cpu3 m | psi dimensions {4, 256, 256} | iter 1/3 |  time   0.4296 s | GOp/s 104.1811
     * [2020-10-05 15:21:27][tensorbench][  info  ] H²|Ψ> version cpu3 m | psi dimensions {4, 256, 256} | iter 2/3 |  time   0.4258 s | GOp/s 105.0995
     * [2020-10-05 15:21:27][tensorbench][  info  ] H²|Ψ> version cpu3 m | psi dimensions {4, 256, 256} | iter 3/3 |  time   0.4179 s | GOp/s 107.0952
     * [2020-10-05 15:21:27][tensorbench][  info  ] H²|Ψ> version cpu3 m | total time 1.2733 s
     * [2020-10-05 15:21:28][tensorbench][  info  ] H²|Ψ> version cpu3 d | psi dimensions {4, 256, 256} | iter 1/3 |  time   0.4572 s | GOp/s 98.4875
     * [2020-10-05 15:21:28][tensorbench][  info  ] H²|Ψ> version cpu3 d | psi dimensions {4, 256, 256} | iter 2/3 |  time   0.4225 s | GOp/s 106.5780
     * [2020-10-05 15:21:29][tensorbench][  info  ] H²|Ψ> version cpu3 d | psi dimensions {4, 256, 256} | iter 3/3 |  time   0.4460 s | GOp/s 100.9666
     * [2020-10-05 15:21:29][tensorbench][  info  ] H²|Ψ> version cpu3 d | total time 1.3257 s
     *
     *
     * [2020-10-05 15:23:06][tensorbench][  info  ] H²|Ψ> version cpu3 m | psi dimensions {64, 256, 256} | iter 1/3 |  time   9.7150 s | GOp/s 515.3747
     * [2020-10-05 15:23:16][tensorbench][  info  ] H²|Ψ> version cpu3 m | psi dimensions {64, 256, 256} | iter 2/3 |  time   9.7767 s | GOp/s 512.1229
     * [2020-10-05 15:23:25][tensorbench][  info  ] H²|Ψ> version cpu3 m | psi dimensions {64, 256, 256} | iter 3/3 |  time   9.5951 s | GOp/s 521.8139
     * [2020-10-05 15:23:25][tensorbench][  info  ] H²|Ψ> version cpu3 m | total time 29.0868 s
     * [2020-10-05 15:23:35][tensorbench][  info  ] H²|Ψ> version cpu3 d | psi dimensions {64, 256, 256} | iter 1/3 |  time   9.7333 s | GOp/s 106.3428
     * [2020-10-05 15:23:45][tensorbench][  info  ] H²|Ψ> version cpu3 d | psi dimensions {64, 256, 256} | iter 2/3 |  time   9.7678 s | GOp/s 105.9679
     * [2020-10-05 15:23:55][tensorbench][  info  ] H²|Ψ> version cpu3 d | psi dimensions {64, 256, 256} | iter 3/3 |  time   9.7425 s | GOp/s 106.2432
     * [2020-10-05 15:23:55][tensorbench][  info  ] H²|Ψ> version cpu3 d | total time 29.2436 s
     *
     * Notice how in the case d^6=64 and chi = 256 the supposed expected improvement should be ~400% which agrees with point (*) above
     *
     */
}

template<typename Scalar, LagrangeNorm lagrangeNorm>
void bfgs_variance_functor<Scalar, lagrangeNorm>::get_Hn(const VectorType &v) const {
    t_Hn->tic();
    auto v_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(v.derived().data(), dims);
    tools::common::contraction::matrix_vector_product(Hn_tensor, v_tensor, mpo, envL, envR);
    t_Hn->toc();
}

template<typename Scalar, LagrangeNorm lagrangeNorm>
void bfgs_variance_functor<Scalar, lagrangeNorm>::get_H2r(const VectorType &r) const {
    t_H2r->tic();
    auto r_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(r.derived().data(), dims);
    tools::common::contraction::matrix_vector_product(H2r_tensor, r_tensor, mpo2, env2L, env2R);
    t_H2r->toc();
}

template<typename Scalar, LagrangeNorm lagrangeNorm>
void bfgs_variance_functor<Scalar, lagrangeNorm>::compress() {
    if(readyCompress) return;

    svd::settings svd_settings;
    svd_settings.svd_lib        = SVDLib::lapacke;
    svd_settings.use_bdc        = false;
    svd_settings.threshold      = 1e-12;
    svd_settings.switchsize_bdc = 4096;
    svd::solver svd(svd_settings);
    auto        old_dimensions = mpo2.dimensions();
    //    Eigen::Tensor<Scalar, 4> mpo2_l2r;
    {
        // Compress left to right
        auto [U, S, V]                     = svd.split_mpo_l2r(mpo2);
        Eigen::Tensor<Scalar, 3> env2R_tmp = env2R.contract(V, tenx::idx({2}, {1}));
        env2R                              = env2R_tmp;
        mpo2                               = U.contract(tenx::asDiagonal(S), tenx::idx({1}, {0})).shuffle(std::array<long, 4>{0, 3, 1, 2});
    }
    {
        // Compress right to left
        auto [U, S, V]                     = svd.split_mpo_r2l(mpo2);
        Eigen::Tensor<Scalar, 3> env2L_tmp = env2L.contract(U, tenx::idx({2}, {0}));
        env2L                              = env2L_tmp;
        mpo2                               = tenx::asDiagonal(S).contract(V, tenx::idx({1}, {0}));
    }
    readyCompress = true;
    tools::log->debug("Compressed MPO dimensions {} -> {}", old_dimensions, mpo2.dimensions());
}

template<typename Scalar, LagrangeNorm lagrangeNorm>
void bfgs_variance_functor<Scalar, lagrangeNorm>::set_shift(double shift_) {
    if(readyShift) return; // This only happens once!!
    if(readyCompress) throw std::runtime_error("Cannot shift the mpo: it is already compressed!");
    shift = shift_;
    tools::log->debug("Setting MPO shift <(H²-σ)> - <(H-σ)>² with σ = {:.16f}", shift_);

    // The MPO is a rank4 tensor ijkl where the first 2 ij indices draw a simple
    // rank2 matrix, where each element is also a matrix with the size
    // determined by the last 2 indices kl.
    // When we shift an MPO, all we do is subtract a diagonal matrix from
    // the botton left corner of the ij-matrix.

    // Start with mpo (hamiltonian)
    {
        // Setup extents and handy objects
        auto                                       shape = mpo.dimensions();
        std::array<long, 4>                        offset4{shape[0] - 1, 0, 0, 0};
        std::array<long, 4>                        extent4{1, 1, shape[2], shape[3]};
        std::array<long, 2>                        extent2{shape[2], shape[3]};
        MatrixType                                 sigma_Id = shift * MatrixType::Identity(extent2[0], extent2[1]);
        Eigen::TensorMap<Eigen::Tensor<Scalar, 2>> sigma_Id_map(sigma_Id.data(), sigma_Id.rows(), sigma_Id.cols());
        mpo.slice(offset4, extent4).reshape(extent2) -= sigma_Id_map;
    }
    // Next is with MPO² (hamiltonian squared). Here we must shift²
    {
        // Setup extents and handy objects
        auto                                       shape = mpo2.dimensions();
        std::array<long, 4>                        offset4{shape[0] - 1, 0, 0, 0};
        std::array<long, 4>                        extent4{1, 1, shape[2], shape[3]};
        std::array<long, 2>                        extent2{shape[2], shape[3]};
        MatrixType                                 sigma_Id = shift * shift * MatrixType::Identity(extent2[0], extent2[1]);
        Eigen::TensorMap<Eigen::Tensor<Scalar, 2>> sigma_Id_map(sigma_Id.data(), sigma_Id.rows(), sigma_Id.cols());
        mpo2.slice(offset4, extent4).reshape(extent2) -= sigma_Id_map;
    }

    readyShift = true;
}

template class tools::finite::opt::internal::bfgs_variance_functor<real, LagrangeNorm::ON>;
template class tools::finite::opt::internal::bfgs_variance_functor<real, LagrangeNorm::OFF>;
template class tools::finite::opt::internal::bfgs_variance_functor<cplx, LagrangeNorm::ON>;
template class tools::finite::opt::internal::bfgs_variance_functor<cplx, LagrangeNorm::OFF>;
