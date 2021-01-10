//
// Created by david on 2019-07-15.
//

#include "ceres_direct_functor.h"
#include <tensors/class_tensors_finite.h>
#include <tools/common/contraction.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt-internal/opt-internal.h>

using namespace tools::finite::opt::internal;

template<typename Scalar>
ceres_direct_functor<Scalar>::ceres_direct_functor(const class_tensors_finite &tensors, const class_algorithm_status &status)
    : ceres_base_functor(tensors, status) {
    tools::log->trace("Constructing direct functor");

    if constexpr(std::is_same<Scalar, double>::value) {
        tools::log->trace("- Generating real-valued multisite components");
        mpo                 = tensors.get_multisite_mpo().real();
        mpo2                = tensors.get_multisite_mpo_squared().real();
        const auto &env_ene = tensors.get_multisite_ene_blk();
        const auto &env_var = tensors.get_multisite_var_blk();
        envL                = env_ene.L.real();
        envR                = env_ene.R.real();
        env2L               = env_var.L.real();
        env2R               = env_var.R.real();
    }

    if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
        tools::log->trace("- Generating complex-valued multisite components");
        mpo                 = tensors.get_multisite_mpo();
        mpo2                = tensors.get_multisite_mpo_squared();
        const auto &env_ene = tensors.get_multisite_ene_blk();
        const auto &env_var = tensors.get_multisite_var_blk();
        envL                = env_ene.L;
        envR                = env_ene.R;
        env2L               = env_var.L;
        env2R               = env_var.R;
    }
    tools::log->trace("- Allocating memory for matrix-vector products");
    dsizes = state.active_dimensions();
    Hv_tensor.resize(dsizes);
    H2v_tensor.resize(dsizes);
    num_parameters = static_cast<int>(dsizes[0] * dsizes[1] * dsizes[2]);
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) { num_parameters *= 2; }
}

template<typename Scalar>
bool ceres_direct_functor<Scalar>::Evaluate(const double *v_double_double, double *fx, double *grad_double_double) const {
    t_step->tic();
    Scalar ene, ene2, var;
    Scalar vHv, vH2v;
    double vv, log10var;
    int    vecSize = NumParameters();
    if constexpr(std::is_same<Scalar, std::complex<double>>::value) { vecSize = NumParameters() / 2; }
    Eigen::Map<const VectorType> v(reinterpret_cast<const Scalar *>(v_double_double), vecSize);
    vv   = v.squaredNorm();
    norm = std::sqrt(vv);

    get_H2v(v);
    get_Hv(v);

    auto Hv  = Eigen::Map<VectorType>(Hv_tensor.data(), Hv_tensor.size());
    auto H2v = Eigen::Map<VectorType>(H2v_tensor.data(), H2v_tensor.size());

    print_path = false;

    t_vHv->tic();
    vHv = v.dot(Hv);
    t_vHv->toc();
    t_vH2v->tic();
    vH2v = v.dot(H2v);
    t_vH2v->toc();

    // Do this next bit carefully to avoid negative variance when numbers are very small
    ene  = vHv / vv;
    ene2 = vH2v / vv;
    if(std::real(ene2) < 0.0) tools::log->debug("Counter = {}. ene2 is negative:  {:.16f} + i {:.16f}", counter, std::real(ene2), std::imag(ene2));
    ene2 = std::real(ene2) < 0.0 ? std::abs(ene2) : std::real(ene2);
    ene2 = std::real(ene2) == 0.0 ? std::numeric_limits<double>::epsilon() : std::real(ene2);

    var = ene2 - ene * ene;
    if(std::real(var) < 0.0) tools::log->debug("Counter = {}. var  is negative:  {:.16f} + i {:.16f}", counter, std::real(var), std::imag(var));
    var = std::real(var) < 0.0 ? std::abs(var) : std::real(var);
    var = std::real(var) == 0.0 ? std::numeric_limits<double>::epsilon() : std::real(var);

    energy                         = std::real(ene + energy_reduced);
    energy_per_site                = energy / static_cast<double>(length);
    variance                       = std::abs(var);
    variance_per_site              = variance / static_cast<double>(length);
    norm_offset                    = std::abs(vv) - 1.0;
    double epsilon                 = 1e-15;
    log10var                       = std::log10(epsilon + variance);

    if(fx != nullptr) { fx[0] = log10var; }

    Eigen::Map<VectorType> grad(reinterpret_cast<Scalar *>(grad_double_double), vecSize);
    if(grad_double_double != nullptr) {
        auto vv_1  = std::pow(vv, -1);
        auto var_1 = (1.0 / (epsilon + var) / std::log(10));
        grad       = var_1 * vv_1 * (H2v - 2.0 * ene * Hv - (ene2 - 2.0 * ene * ene) * v);
        if constexpr(std::is_same<Scalar, double>::value) { grad *= 2.0; }
    }

    if(std::isnan(log10var) or std::isinf(log10var)) {
        tools::log->warn("log₁₀ variance is invalid");
        tools::log->warn("vv              = {:.16f} + i{:.16f}", std::real(vv), std::imag(vv));
        tools::log->warn("vH2v            = {:.16f} + i{:.16f}", std::real(vH2v), std::imag(vH2v));
        tools::log->warn("vHv             = {:.16f} + i{:.16f}", std::real(vHv), std::imag(vHv));
        tools::log->warn("var             = {:.16f} + i{:.16f}", std::real(var), std::imag(var));
        tools::log->warn("ene             = {:.16f} + i{:.16f}", std::real(ene), std::imag(ene));
        tools::log->warn("log₁₀(var)      = {:.16f}", std::log10(variance));
        tools::log->warn("energy offset   = {:.16f}", energy_offset);
        tools::log->warn("norm   offset   = {:.16f}", norm_offset);
        throw std::runtime_error("Direct functor failed at counter = " + std::to_string(counter));
    }

    counter++;
    t_step->toc();
    return true;
}

template<typename Scalar>
void ceres_direct_functor<Scalar>::get_H2v(const VectorType &v) const {
    t_vH2->tic();
    auto v_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(v.derived().data(), dims);
    tools::common::contraction::matrix_vector_product(H2v_tensor, v_tensor, mpo2, env2L, env2R);
    t_vH2->toc();

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
     *      .contract(mpo, Textra::idx({0, 3}, {2, 0}))
     * is equivalent to
     *      .contract(mpo, Textra::idx({3, 0}, {0, 2}))
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

template<typename Scalar>
void ceres_direct_functor<Scalar>::get_Hv(const VectorType &v) const {
    t_vH->tic();
    auto v_tensor = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 3>>(v.derived().data(), dims);
    tools::common::contraction::matrix_vector_product(Hv_tensor, v_tensor, mpo, envL, envR);
    t_vH->toc();
}



template class tools::finite::opt::internal::ceres_direct_functor<double>;
template class tools::finite::opt::internal::ceres_direct_functor<std::complex<double>>;
