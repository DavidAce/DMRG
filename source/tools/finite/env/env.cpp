#include "../env.h"
#include "config/debug.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "EnvExpansionResult.h"
#include "math/linalg/matrix.h"
#include "math/num.h"
#include "math/svd.h"
#include "math/tenx.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/env/EnvEne.h"
#include "tensors/site/env/EnvVar.h"
#include "tensors/site/mpo/MpoSite.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/mps.h"
#include <Eigen/Eigenvalues>

namespace settings {
    static constexpr bool debug_edges     = false;
    static constexpr bool debug_expansion = false;
}

std::array<double, 2> tools::finite::env::get_optimal_mixing_factor_ene(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                        const EdgesFinite &edges, EnvExpandMode envExpandMode) {
    if(not has_flag(envExpandMode, EnvExpandMode::ENE)) return {0.0, 0.0};
    assert_edges_ene(state, model, edges);

    auto                                t    = state.get_multisite_mps(sites);
    auto                                env1 = edges.get_multisite_env_ene_blk(sites);
    auto                                mpos = model.get_mpo(sites);
    std::vector<Eigen::Tensor<cplx, 4>> H1;
    for(const auto &mpo : mpos) H1.emplace_back(mpo.get().MPO());

    auto H1t = tools::common::contraction::matrix_vector_product(t, H1, env1.L, env1.R);

    auto v0     = tenx::VectorMap(t);
    auto H1v0   = tenx::VectorMap(H1t);
    auto v0H1v0 = v0.dot(H1v0);
    // Gram schmidt
    Eigen::VectorXcd       v1     = (H1v0 - v0 * v0.dot(H1v0)).normalized(); // Normalized residual vector
    Eigen::Tensor<cplx, 3> t1     = Eigen::TensorMap<Eigen::Tensor<cplx, 3>>(v1.data(), t.dimensions());
    auto                   H1t1   = tools::common::contraction::matrix_vector_product(t1, H1, env1.L, env1.R);
    Eigen::VectorXcd       H1v1   = tenx::VectorMap(H1t1);
    auto                   v1H1v0 = v1.dot(H1v0);
    auto                   v1H1v1 = v1.dot(H1v1);
    auto                   v0H1v1 = std::conj(v1H1v0);

    auto K      = Eigen::Matrix2cd{{v0H1v0, v0H1v1}, {v1H1v0, v1H1v1}};
    auto solver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix2cd>(K, Eigen::ComputeEigenvectors);
    auto evals  = solver.eigenvalues();
    auto evecs  = solver.eigenvectors();
    tools::log->info("K(H1): \n{}\n", linalg::matrix::to_string(K.real(), 8));
    tools::log->info("evecs: \n{}\n", linalg::matrix::to_string(evecs, 8));
    tools::log->info("evals: \n{}\n", linalg::matrix::to_string(evals, 8));
    long                  minIdx  = 0;
    [[maybe_unused]] auto minEval = evals.minCoeff(&minIdx);
    Eigen::VectorXd       col     = evecs.col(minIdx).real();

    return {col.coeff(0), col.coeff(1)};

    // Calculate the residual_norm r = |Hv - Ev|
    // auto                                t    = state.get_multisite_mps(sites);
    // auto                                env1  = edges.get_multisite_env_ene_blk(sites);
    // auto                                mpos = model.get_mpo(sites);
    // std::vector<Eigen::Tensor<cplx, 4>> H1;
    // for(const auto &mpo : mpos) H1.emplace_back(mpo.get().MPO());
    //
    // auto Hv    = tools::common::contraction::matrix_vector_product(v, H, env1.L, env1.R);
    // auto HHv   = tools::common::contraction::matrix_vector_product(Hv, H, env.L, env1.R);
    // auto vHv   = tools::common::contraction::contract_mps_overlap(v, Hv);
    // auto vHHHv = tools::common::contraction::contract_mps_overlap(Hv, HHv);
    // auto E     = vHv;                                                              // Energy
    // auto r     = cplx((tenx::VectorMap(Hv) - tenx::VectorMap(v) * E).norm(), 0.0); // Residual norm
    // // auto V  =
    // cplx Er = vHHHv / (r * r); // "Energy" of the residual vector
    //
    // auto K   = Eigen::Matrix2cd{{E, r}, {r, Er}};
    // auto evs = Eigen::SelfAdjointEigenSolver<Eigen::Matrix2cd>(K, Eigen::ComputeEigenvectors).eigenvectors();
    // return evs.col(0).cwiseAbs().minCoeff();
}

std::array<double, 2> tools::finite::env::get_optimal_mixing_factor_var(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                        const EdgesFinite &edges, EnvExpandMode envExpandMode) {
    if(not has_flag(envExpandMode, EnvExpandMode::VAR)) return {0.0, 0.0};
    // Calculate the residual_norm r = |Hv - Ev|
    assert_edges_var(state, model, edges);
    auto                                t    = state.get_multisite_mps(sites);
    auto                                env2 = edges.get_multisite_env_var_blk(sites);
    auto                                mpos = model.get_mpo(sites);
    std::vector<Eigen::Tensor<cplx, 4>> H2;
    for(const auto &mpo : mpos) H2.emplace_back(mpo.get().MPO2());

    auto H2t = tools::common::contraction::matrix_vector_product(t, H2, env2.L, env2.R);
    // Define the krylov subspace
    auto v0   = tenx::VectorMap(t);
    auto H2v0 = tenx::VectorMap(H2t);

    // Gram schmidt
    Eigen::VectorXcd v1 = (H2v0 - v0 * v0.dot(H2v0)).normalized(); // Orthogonal vector to v0

    Eigen::Tensor<cplx, 3> t1     = Eigen::TensorMap<Eigen::Tensor<cplx, 3>>(v1.data(), t.dimensions());
    auto                   H2t1   = tools::common::contraction::matrix_vector_product(t1, H2, env2.L, env2.R);
    auto                   H2v1   = tenx::VectorMap(H2t1);
    auto                   v0H2v0 = v0.dot(H2v0);
    auto                   v1H2v0 = v1.dot(H2v0);
    auto                   v1H2v1 = v1.dot(H2v1);
    auto                   v0H2v1 = std::conj(v1H2v0);

    auto K = Eigen::Matrix2cd{
        {v0H2v0, v0H2v1},
        {v1H2v0, v1H2v1},
    };

    auto solver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix2cd>(K, Eigen::ComputeEigenvectors);
    auto evals  = solver.eigenvalues();
    auto evecs  = solver.eigenvectors();
    tools::log->info("K(H2): \n{}\n", linalg::matrix::to_string(K.real(), 8));
    tools::log->info("evecs: \n{}\n", linalg::matrix::to_string(evecs, 8));
    tools::log->info("evals: \n{}\n", linalg::matrix::to_string(evals, 8));
    long                  minIdx  = 0;
    [[maybe_unused]] auto minEval = evals.minCoeff(&minIdx);
    Eigen::VectorXd       col     = evecs.col(minIdx).real();
    // if(col.coeff(0) < 0.0) col.array() *= -1.0;

    return {col.coeff(0), col.coeff(1)};

    // Calculate the residual_norm r = |Hv - Ev|

    // auto env = edges.get_multisite_env_var_blk(sites);

    // auto H2H2t    = tools::common::contraction::matrix_vector_product(H2t, H2, env.L, env.R);
    // auto vH2v     = tools::common::contraction::contract_mps_overlap(t, H2t);
    // auto vH2H2H2v = tools::common::contraction::contract_mps_overlap(H2t, H2H2t);

    // auto v   = tenx::VectorMap(t);
    // auto H2v = tenx::VectorMap(H2t);
    // auto H2H2v = tenx::VectorMap(H2H2t);

    // auto vH2v     = v.dot(H2v);
    // auto vH2H2H2v = H2v.dot(H2H2v);

    // auto V  = vH2v;                            // Variance
    // auto rn = cplx((H2v - v * V).norm(), 0.0); // Residual norm
    // cplx Vr = vH2H2H2v / (rn * rn);            // "Variance" of the residual vector
    //
    // auto Kold = Eigen::Matrix2cd{{V, rn}, {rn, Vr}};
    // auto evs  = Eigen::SelfAdjointEigenSolver<Eigen::Matrix2cd>(K, Eigen::ComputeEigenvectors).eigenvectors();
    // tools::log->info("K  : \n{}\n", linalg::matrix::to_string(K.real(), 8));
    // tools::log->info("evs: \n{}\n", linalg::matrix::to_string(evs, 8));

    // auto r  = Eigen::VectorXcd((H2v - v * V).normalized()); // Residual vector
    // // auto rt = tenx::TensorMap(r, t.dimensions()); // Residual tensor
    // auto rt = Eigen::TensorMap<Eigen::Tensor<cplx,3>>(r.data(), t.dimensions()); // Residual tensor
    // auto H2rt = Eigen::Tensor<cplx,3>(t.dimensions());
    // tools::common::contraction::matrix_vector_product(H2rt, rt, H2, env.L, env.R);
    // auto H2r = tenx::VectorMap(H2rt);
    //
    // auto rH2v = r.dot(H2v);
    // auto rH2r = r.dot(H2r);
    // auto K2   = Eigen::Matrix2cd{{vH2v, rH2v}, {std::conj(rH2v), rH2r}};
    // auto evs2 = Eigen::SelfAdjointEigenSolver<Eigen::Matrix2cd>(K2, Eigen::ComputeEigenvectors).eigenvectors();
    //
    // tools::log->info("Kold  : \n{}\n", linalg::matrix::to_string(Kold.real(), 8));
    // tools::log->info("evs2: \n{}\n", linalg::matrix::to_string(evs, 8));
    // return evs.col(0).cwiseAbs().minCoeff();
}

std::tuple<double, double, double> tools::finite::env::get_optimal_mixing_factors(const std::vector<size_t> &sites, const StateFinite &state,
                                                                                  const ModelFinite &model, const EdgesFinite &edges,
                                                                                  EnvExpandMode envExpandMode) {
    assert_edges_ene(state, model, edges);
    assert_edges_var(state, model, edges);
    auto                                t0   = state.get_multisite_mps(sites);
    auto                                env1 = edges.get_multisite_env_ene_blk(sites);
    auto                                env2 = edges.get_multisite_env_var_blk(sites);
    auto                                mpos = model.get_mpo(sites);
    std::vector<Eigen::Tensor<cplx, 4>> H1;
    std::vector<Eigen::Tensor<cplx, 4>> H2;
    for(const auto &mpo : mpos) H1.emplace_back(mpo.get().MPO());
    for(const auto &mpo : mpos) H2.emplace_back(mpo.get().MPO2());

    auto H1t = tools::common::contraction::matrix_vector_product(t0, H1, env1.L, env1.R);
    auto H2t = tools::common::contraction::matrix_vector_product(t0, H2, env2.L, env2.R);
    // Define the krylov subspace
    auto v0   = tenx::VectorMap(t0);
    auto H1v0 = tenx::VectorMap(H1t);
    auto H2v0 = tenx::VectorMap(H2t);

    // Gram schmidt
    Eigen::VectorXcd v1 = (H1v0 - v0 * v0.dot(H1v0)).normalized();                     // Normalized residual vector
    Eigen::VectorXcd v2 = (H2v0 - v0 * v0.dot(H2v0) - v1 * v1.dot(H2v0)).normalized(); // Normalized residual vector

    Eigen::Tensor<cplx, 3> t1     = Eigen::TensorMap<Eigen::Tensor<cplx, 3>>(v1.data(), t0.dimensions());
    Eigen::Tensor<cplx, 3> t2     = Eigen::TensorMap<Eigen::Tensor<cplx, 3>>(v2.data(), t0.dimensions());
    auto                   H1t1   = tools::common::contraction::matrix_vector_product(t1, H1, env1.L, env1.R);
    auto                   H1t2   = tools::common::contraction::matrix_vector_product(t2, H1, env1.L, env1.R);
    auto                   H2t1   = tools::common::contraction::matrix_vector_product(t1, H2, env2.L, env2.R);
    auto                   H2t2   = tools::common::contraction::matrix_vector_product(t2, H2, env2.L, env2.R);
    auto                   H1v1   = tenx::VectorMap(H1t1);
    auto                   H1v2   = tenx::VectorMap(H1t2);
    auto                   H2v1   = tenx::VectorMap(H2t1);
    auto                   H2v2   = tenx::VectorMap(H2t2);
    auto                   H1_on  = static_cast<double>(has_flag(envExpandMode, EnvExpandMode::ENE));
    auto                   H2_on  = static_cast<double>(has_flag(envExpandMode, EnvExpandMode::VAR));
    auto                   v0H1v0 = v0.dot(H1v0);
    auto                   v1H1v0 = v1.dot(H1v0) * H1_on;
    auto                   v1H1v1 = v1.dot(H1v1) * H1_on;
    auto                   v1H1v2 = v1.dot(H1v2) * H1_on * H2_on;
    auto                   v2H1v0 = v2.dot(H1v0) * H2_on;
    auto                   v2H1v2 = v2.dot(H1v2) * H2_on * H2_on;

    auto v0H1v1 = std::conj(v1H1v0);
    auto v0H1v2 = std::conj(v2H1v0);
    auto v2H1v1 = std::conj(v1H1v2);

    auto K1 = Eigen::Matrix3cd{
        {v0H1v0, v0H1v1, v0H1v2},
        {v1H1v0, v1H1v1, v1H1v2},
        {v2H1v0, v2H1v1, v2H1v2},
    };

    auto solver1 = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd>(K1, Eigen::ComputeEigenvectors);
    auto evals1  = solver1.eigenvalues();
    auto evecs1  = solver1.eigenvectors();
    tools::log->info("K1   : \n{}\n", linalg::matrix::to_string(K1.real(), 8));
    tools::log->info("evecs: \n{}\n", linalg::matrix::to_string(evecs1, 8));
    tools::log->info("evals: \n{}\n", linalg::matrix::to_string(evals1, 8));

    auto v0H2v0 = v0.dot(H2v0);
    auto v1H2v0 = v1.dot(H2v0);
    auto v2H2v0 = v2.dot(H2v0);
    auto v1H2v1 = v1.dot(H2v1);
    auto v1H2v2 = v1.dot(H2v2);
    auto v2H2v2 = v2.dot(H2v2);

    auto v0H2v1 = std::conj(v1H2v0);
    auto v0H2v2 = std::conj(v2H2v0);
    auto v2H2v1 = std::conj(v1H2v2);

    auto K2 = Eigen::Matrix3cd{
        {v0H2v0, v0H2v1, v0H2v2},
        {v1H2v0, v1H2v1, v1H2v2},
        {v2H2v0, v2H2v1, v2H2v2},
    };

    auto solver2 = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd>(K2, Eigen::ComputeEigenvectors);
    auto evals2  = solver2.eigenvalues();
    auto evecs2  = solver2.eigenvectors();
    tools::log->info("K2   : \n{}\n", linalg::matrix::to_string(K2.real(), 8));
    tools::log->info("evecs: \n{}\n", linalg::matrix::to_string(evecs2, 8));
    tools::log->info("evals: \n{}\n", linalg::matrix::to_string(evals2, 8));

    auto K3      = Eigen::Matrix3cd(K2 - K1 * K1);
    auto solver3 = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd>(K3, Eigen::ComputeEigenvectors);
    auto evals3  = solver3.eigenvalues();
    auto evecs3  = solver3.eigenvectors();
    tools::log->info("K3   : \n{}\n", linalg::matrix::to_string(K3.real(), 8));
    tools::log->info("evecs: \n{}\n", linalg::matrix::to_string(evecs3, 8));
    tools::log->info("evals: \n{}\n", linalg::matrix::to_string(evals3, 8));

    // Decide to minimize <H²> or <H²> - E²

    long            minIdx  = 0;
    auto            minEval = evals3.minCoeff(&minIdx);
    Eigen::VectorXd col     = evecs3.col(minIdx).real();
    // if(col.coeff(0) < 0.0) col.array() *= -1.0;
    Eigen::VectorXcd       vf     = (col.coeff(0) * v0 + col.coeff(1) * v1 + col.coeff(2) * v2).normalized();
    Eigen::Tensor<cplx, 3> tf     = Eigen::TensorMap<Eigen::Tensor<cplx, 3>>(vf.data(), t0.dimensions());
    auto                   vfH1vf = tools::common::contraction::expectation_value(tf, H1, env1.L, env1.R);
    auto                   vfH2vf = tools::common::contraction::expectation_value(tf, H2, env2.L, env2.R);
    tools::log->info("final state <H> = {:.16f}, <H²> = {:.16f}, Var(H) = {:.3e}", vfH1vf, vfH2vf, vfH2vf - vfH1vf * vfH1vf);

    return {col.coeff(0), col.coeff(1), col.coeff(2)};
}

EnvExpansionResult tools::finite::env::internal::get_optimally_mixed_block_H1(const std::vector<size_t> &sites, const StateFinite &state,
                                                                              const ModelFinite &model, const EdgesFinite &edges, OptRitz ritz) {
    assert_edges_ene(state, model, edges);
    assert_edges_var(state, model, edges);
    auto res         = EnvExpansionResult();
    res.sites        = sites;
    res.dims_old     = state.get_mps_dims(sites);
    res.bond_old     = state.get_bond_dims(sites);
    res.posL         = sites.front();
    res.posR         = sites.back();
    const auto &mpsL = state.get_mps_site(res.posL);
    const auto &mpsR = state.get_mps_site(res.posR);
    res.dimL_old     = mpsL.dimensions();
    res.dimR_old     = mpsR.dimensions();
    auto t0          = state.get_multisite_mps(sites);
    auto t1          = Eigen::Tensor<cplx, 3>(t0.dimensions());
    auto t2          = Eigen::Tensor<cplx, 3>(t0.dimensions());
    auto v0          = tenx::VectorMap(t0);
    auto v1          = tenx::VectorMap(t1);
    auto v2          = tenx::VectorMap(t2);

    auto                                env1 = edges.get_multisite_env_ene_blk(sites);
    auto                                mpos = model.get_mpo(sites);
    std::vector<Eigen::Tensor<cplx, 4>> H1;
    std::vector<Eigen::Tensor<cplx, 4>> H2;
    for(const auto &mpo : mpos) H1.emplace_back(mpo.get().MPO());

    auto H1t0 = tools::common::contraction::matrix_vector_product(t0, H1, env1.L, env1.R);
    // Define the krylov subspace
    auto H1v0 = tenx::VectorMap(H1t0);

    // Gram schmidt
    v1 = (H1v0 - v0 * v0.dot(H1v0)).normalized();

    auto H1t1   = tools::common::contraction::matrix_vector_product(t1, H1, env1.L, env1.R);
    auto H1v1   = tenx::VectorMap(H1t1);
    auto v0H1v0 = v0.dot(H1v0);
    auto v1H1v0 = v1.dot(H1v0);
    auto v0H1v1 = std::conj(v1H1v0);
    auto v1H1v1 = v1.dot(H1v1);

    auto K1 = Eigen::Matrix2cd{
        {v0H1v0, v0H1v1},
        {v1H1v0, v1H1v1},
    };

    auto   solver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix2cd>(K1, Eigen::ComputeEigenvectors);
    auto   evals  = solver.eigenvalues();
    auto   evecs  = solver.eigenvectors();
    long   optIdx = 0;
    double optVal = 0;
    switch(ritz) {
        case OptRitz::SR: {
            optVal = evals.minCoeff(&optIdx);
            break;
        }
        case OptRitz::SM: {
            optVal = evals.cwiseAbs().minCoeff(&optIdx);
            break;
        }
        case OptRitz::LR: {
            optVal = evals.maxCoeff(&optIdx);
            break;
        }
        case OptRitz::LM: {
            optVal = evals.cwiseAbs().maxCoeff(&optIdx);
            break;
        }
        default: {
            // Take the closest to the current energy
            optVal = (evals.array() - std::real(v0H1v0)).cwiseAbs().minCoeff(&optIdx);
            optVal = evals.coeff(optIdx);
            break;
        }
    }

    Eigen::VectorXd col = evecs.col(optIdx).normalized().real();
    if(col.hasNaN()) {
        res.mixed_blk = t0;
        res.alpha_mps = 1.0;
        res.alpha_ene = 0.0;
        res.alpha_var = 0.0;
    } else {
        res.mixed_blk = Eigen::Tensor<cplx, 3>(t0.dimensions());
        res.alpha_mps = col.coeff(0);
        res.alpha_ene = col.coeff(1);
        res.alpha_var = 0.0;
        auto vf       = tenx::VectorMap(res.mixed_blk);
        vf            = (res.alpha_mps * v0 + res.alpha_ene * v1).normalized();
    }

    tools::log->debug("mixed state result: <H> = {:.16f}", optVal);

    return res;
}
EnvExpansionResult tools::finite::env::internal::get_optimally_mixed_block_H2(const std::vector<size_t> &sites, const StateFinite &state,
                                                                              const ModelFinite &model, const EdgesFinite &edges, OptRitz ritz) {
    assert_edges_ene(state, model, edges);
    assert_edges_var(state, model, edges);
    auto res         = EnvExpansionResult();
    res.sites        = sites;
    res.dims_old     = state.get_mps_dims(sites);
    res.bond_old     = state.get_bond_dims(sites);
    res.posL         = sites.front();
    res.posR         = sites.back();
    const auto &mpsL = state.get_mps_site(res.posL);
    const auto &mpsR = state.get_mps_site(res.posR);
    res.dimL_old     = mpsL.dimensions();
    res.dimR_old     = mpsR.dimensions();
    auto t0          = state.get_multisite_mps(sites);
    auto t2          = Eigen::Tensor<cplx, 3>(t0.dimensions());
    auto v0          = tenx::VectorMap(t0);
    auto v2          = tenx::VectorMap(t2);

    auto                                env1 = edges.get_multisite_env_ene_blk(sites);
    auto                                env2 = edges.get_multisite_env_var_blk(sites);
    auto                                mpos = model.get_mpo(sites);
    std::vector<Eigen::Tensor<cplx, 4>> H2;
    for(const auto &mpo : mpos) H2.emplace_back(mpo.get().MPO2());

    auto H2t0 = tools::common::contraction::matrix_vector_product(t0, H2, env2.L, env2.R);
    // Define the krylov subspace
    auto H2v0 = tenx::VectorMap(H2t0);

    // Gram schmidt
    v2 = (H2v0 - v0 * v0.dot(H2v0)).normalized();

    auto H2t2 = tools::common::contraction::matrix_vector_product(t2, H2, env2.L, env2.R);
    auto H2v2 = tenx::VectorMap(H2t2);

    auto v0H2v0 = v0.dot(H2v0);
    auto v2H2v0 = v2.dot(H2v0);
    auto v2H2v2 = v2.dot(H2v2);
    auto v0H2v2 = std::conj(v2H2v0);

    auto K2 = Eigen::Matrix2cd{
        {v0H2v0, v0H2v2},
        {v2H2v0, v2H2v2},
    };

    auto   solver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix2cd>(K2, Eigen::ComputeEigenvectors);
    auto   evals  = solver.eigenvalues();
    auto   evecs  = solver.eigenvectors();
    long   optIdx = 0;
    double optVal = 0;
    switch(ritz) {
        case OptRitz::SR: {
            optVal = evals.minCoeff(&optIdx);
            break;
        }
        case OptRitz::SM: {
            optVal = evals.cwiseAbs().minCoeff(&optIdx);
            break;
        }
        case OptRitz::LR: {
            optVal = evals.maxCoeff(&optIdx);
            break;
        }
        case OptRitz::LM: {
            optVal = evals.cwiseAbs().maxCoeff(&optIdx);
            break;
        }
        default: {
            // Take the closest to the current H2
            optVal = (evals.array() - std::real(v0H2v0)).cwiseAbs().minCoeff(&optIdx);
            optVal = evals.coeff(optIdx);
            break;
        }
    }

    Eigen::VectorXd col = evecs.col(optIdx).normalized().real();
    if(col.hasNaN()) {
        res.mixed_blk = t0;
        res.alpha_mps = 1.0;
        res.alpha_ene = 0.0;
        res.alpha_var = 0.0;
    } else {
        res.mixed_blk = Eigen::Tensor<cplx, 3>(t0.dimensions());
        res.alpha_mps = col.coeff(0);
        res.alpha_ene = 0.0;
        res.alpha_var = col.coeff(1);
        auto vf       = tenx::VectorMap(res.mixed_blk);
        vf            = (res.alpha_mps * v0 + res.alpha_var * v2).normalized();
    }
    tools::log->debug("mixed state result: <H²> = {:.16f}", optVal);
    return res;
}

EnvExpansionResult tools::finite::env::internal::get_optimally_mixed_block_VarH(const std::vector<size_t> &sites, const StateFinite &state,
                                                                                const ModelFinite &model, const EdgesFinite &edges, OptRitz ritz) {
    assert_edges_ene(state, model, edges);
    assert_edges_var(state, model, edges);
    auto res         = EnvExpansionResult();
    res.sites        = sites;
    res.dims_old     = state.get_mps_dims(sites);
    res.bond_old     = state.get_bond_dims(sites);
    res.posL         = sites.front();
    res.posR         = sites.back();
    const auto &mpsL = state.get_mps_site(res.posL);
    const auto &mpsR = state.get_mps_site(res.posR);
    res.dimL_old     = mpsL.dimensions();
    res.dimR_old     = mpsR.dimensions();
    auto t0          = state.get_multisite_mps(sites);
    auto t1          = Eigen::Tensor<cplx, 3>(t0.dimensions());
    auto t2          = Eigen::Tensor<cplx, 3>(t0.dimensions());
    auto v0          = tenx::VectorMap(t0);
    auto v1          = tenx::VectorMap(t1);
    auto v2          = tenx::VectorMap(t2);

    auto                                env1 = edges.get_multisite_env_ene_blk(sites);
    auto                                env2 = edges.get_multisite_env_var_blk(sites);
    auto                                mpos = model.get_mpo(sites);
    std::vector<Eigen::Tensor<cplx, 4>> H1;
    std::vector<Eigen::Tensor<cplx, 4>> H2;
    for(const auto &mpo : mpos) H1.emplace_back(mpo.get().MPO());
    for(const auto &mpo : mpos) H2.emplace_back(mpo.get().MPO2());

    auto H1t0 = tools::common::contraction::matrix_vector_product(t0, H1, env1.L, env1.R);
    auto H2t0 = tools::common::contraction::matrix_vector_product(t0, H2, env2.L, env2.R);
    // Define the krylov subspace
    auto H1v0 = tenx::VectorMap(H1t0);
    auto H2v0 = tenx::VectorMap(H2t0);

    // Gram schmidt
    v1 = (H1v0 - v0 * v0.dot(H1v0)).normalized();
    v2 = (H2v0 - v0 * v0.dot(H2v0) - v1 * v1.dot(H2v0)).normalized();

    auto H1t1   = tools::common::contraction::matrix_vector_product(t1, H1, env1.L, env1.R);
    auto H1t2   = tools::common::contraction::matrix_vector_product(t2, H1, env1.L, env1.R);
    auto H2t1   = tools::common::contraction::matrix_vector_product(t1, H2, env2.L, env2.R);
    auto H2t2   = tools::common::contraction::matrix_vector_product(t2, H2, env2.L, env2.R);
    auto H1v1   = tenx::VectorMap(H1t1);
    auto H1v2   = tenx::VectorMap(H1t2);
    auto H2v1   = tenx::VectorMap(H2t1);
    auto H2v2   = tenx::VectorMap(H2t2);
    auto v0H1v0 = v0.dot(H1v0);
    auto v1H1v0 = v1.dot(H1v0);
    auto v1H1v1 = v1.dot(H1v1);
    auto v1H1v2 = v1.dot(H1v2);
    auto v2H1v0 = v2.dot(H1v0);
    auto v2H1v2 = v2.dot(H1v2);

    auto v0H1v1 = std::conj(v1H1v0);
    auto v0H1v2 = std::conj(v2H1v0);
    auto v2H1v1 = std::conj(v1H1v2);

    auto K1 = Eigen::Matrix3cd{
        {v0H1v0, v0H1v1, v0H1v2},
        {v1H1v0, v1H1v1, v1H1v2},
        {v2H1v0, v2H1v1, v2H1v2},
    };

    auto v0H2v0 = v0.dot(H2v0);
    auto v1H2v0 = v1.dot(H2v0);
    auto v2H2v0 = v2.dot(H2v0);
    auto v1H2v1 = v1.dot(H2v1);
    auto v1H2v2 = v1.dot(H2v2);
    auto v2H2v2 = v2.dot(H2v2);

    auto v0H2v1 = std::conj(v1H2v0);
    auto v0H2v2 = std::conj(v2H2v0);
    auto v2H2v1 = std::conj(v1H2v2);

    auto K2 = Eigen::Matrix3cd{
        {v0H2v0, v0H2v1, v0H2v2},
        {v1H2v0, v1H2v1, v1H2v2},
        {v2H2v0, v2H2v1, v2H2v2},
    };

    auto   K3     = Eigen::Matrix3cd(K2 - K1 * K1);
    auto   solver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3cd>(K3, Eigen::ComputeEigenvectors);
    auto   evals  = solver.eigenvalues();
    auto   evecs  = solver.eigenvectors();
    long   optIdx = 0;
    double optVal = 0;
    switch(ritz) {
        case OptRitz::SR: {
            optVal = evals.minCoeff(&optIdx);
            break;
        }
        case OptRitz::SM: {
            optVal = evals.cwiseAbs().minCoeff(&optIdx);
            break;
        }
        case OptRitz::LR: {
            optVal = evals.maxCoeff(&optIdx);
            break;
        }
        case OptRitz::LM: {
            optVal = evals.cwiseAbs().maxCoeff(&optIdx);
            break;
        }
        default: {
            // Take the closest to the current VarH
            optVal = (evals.array() - std::real(v0H2v0 - v0H1v0 * v0H1v0)).cwiseAbs().minCoeff(&optIdx);
            optVal = evals.coeff(optIdx);
            break;
        }
    }

    Eigen::VectorXd col = evecs.col(optIdx).normalized().real();
    if(col.hasNaN()) {
        res.mixed_blk = t0;
        res.alpha_mps = 1.0;
        res.alpha_ene = 0.0;
        res.alpha_var = 0.0;
    } else {
        res.mixed_blk = Eigen::Tensor<cplx, 3>(t0.dimensions());
        res.alpha_mps = col.coeff(0);
        res.alpha_ene = col.coeff(1);
        res.alpha_var = col.coeff(2);
        auto vf       = tenx::VectorMap(res.mixed_blk);
        vf            = (col.coeff(0) * v0 + col.coeff(1) * v1 + col.coeff(2) * v2).normalized();
    }

    tools::log->debug("mixed state result: <H²>-<H>² = {:.16f}", optVal);
    return res;
}

EnvExpansionResult tools::finite::env::internal::get_optimally_mixed_block_GsiH(const std::vector<size_t> &sites, const StateFinite &state,
                                                                                const ModelFinite &model, const EdgesFinite &edges, OptRitz ritz) {
    auto algo = state.get_algorithm();
    assert_edges_ene(state, model, edges);
    assert_edges_var(state, model, edges);
    auto res         = EnvExpansionResult();
    res.sites        = sites;
    res.dims_old     = state.get_mps_dims(sites);
    res.bond_old     = state.get_bond_dims(sites);
    res.posL         = sites.front();
    res.posR         = sites.back();
    const auto &mpsL = state.get_mps_site(res.posL);
    const auto &mpsR = state.get_mps_site(res.posR);
    res.dimL_old     = mpsL.dimensions();
    res.dimR_old     = mpsR.dimensions();
    auto t0          = state.get_multisite_mps(sites);
    auto t1          = Eigen::Tensor<cplx, 3>(t0.dimensions());
    auto t2          = Eigen::Tensor<cplx, 3>(t0.dimensions());
    auto v0          = tenx::VectorMap(t0);
    auto v1          = tenx::VectorMap(t1);
    auto v2          = tenx::VectorMap(t2);

    auto                                env1 = edges.get_multisite_env_ene_blk(sites);
    auto                                env2 = edges.get_multisite_env_var_blk(sites);
    auto                                mpos = model.get_mpo(sites);
    std::vector<Eigen::Tensor<cplx, 4>> H1;
    std::vector<Eigen::Tensor<cplx, 4>> H2;
    for(const auto &mpo : mpos) H1.emplace_back(mpo.get().MPO());
    for(const auto &mpo : mpos) H2.emplace_back(mpo.get().MPO2());

    auto H1t0 = tools::common::contraction::matrix_vector_product(t0, H1, env1.L, env1.R);
    auto H2t0 = tools::common::contraction::matrix_vector_product(t0, H2, env2.L, env2.R);
    // Define the krylov subspace
    auto H1v0 = tenx::VectorMap(H1t0);
    auto H2v0 = tenx::VectorMap(H2t0);

    // Gram schmidt
    v1 = (H1v0 - v0 * v0.dot(H1v0)).normalized();
    v2 = (H2v0 - v0 * v0.dot(H2v0) - v1 * v1.dot(H2v0)).normalized();

    auto H1t1   = tools::common::contraction::matrix_vector_product(t1, H1, env1.L, env1.R);
    auto H1t2   = tools::common::contraction::matrix_vector_product(t2, H1, env1.L, env1.R);
    auto H2t1   = tools::common::contraction::matrix_vector_product(t1, H2, env2.L, env2.R);
    auto H2t2   = tools::common::contraction::matrix_vector_product(t2, H2, env2.L, env2.R);
    auto H1v1   = tenx::VectorMap(H1t1);
    auto H1v2   = tenx::VectorMap(H1t2);
    auto H2v1   = tenx::VectorMap(H2t1);
    auto H2v2   = tenx::VectorMap(H2t2);
    auto v0H1v0 = v0.dot(H1v0);
    auto v1H1v0 = v1.dot(H1v0);
    auto v1H1v1 = v1.dot(H1v1);
    auto v1H1v2 = v1.dot(H1v2);
    auto v2H1v0 = v2.dot(H1v0);
    auto v2H1v2 = v2.dot(H1v2);

    auto v0H1v1 = std::conj(v1H1v0);
    auto v0H1v2 = std::conj(v2H1v0);
    auto v2H1v1 = std::conj(v1H1v2);

    auto K1 = Eigen::Matrix3cd{
        {v0H1v0, v0H1v1, v0H1v2},
        {v1H1v0, v1H1v1, v1H1v2},
        {v2H1v0, v2H1v1, v2H1v2},
    };

    auto v0H2v0 = v0.dot(H2v0);
    auto v1H2v0 = v1.dot(H2v0);
    auto v2H2v0 = v2.dot(H2v0);
    auto v1H2v1 = v1.dot(H2v1);
    auto v1H2v2 = v1.dot(H2v2);
    auto v2H2v2 = v2.dot(H2v2);

    auto v0H2v1 = std::conj(v1H2v0);
    auto v0H2v2 = std::conj(v2H2v0);
    auto v2H2v1 = std::conj(v1H2v2);

    auto K2 = Eigen::Matrix3cd{
        {v0H2v0, v0H2v1, v0H2v2},
        {v1H2v0, v1H2v1, v1H2v2},
        {v2H2v0, v2H2v1, v2H2v2},
    };

    auto   solver = Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix3cd>(K1, K2, Eigen::ComputeEigenvectors | Eigen::Ax_lBx);
    auto   evals  = solver.eigenvalues();
    auto   evecs  = solver.eigenvectors();
    long   optIdx = 0;
    double optVal = 0;
    switch(ritz) {
        case OptRitz::SR: {
            optVal = evals.minCoeff(&optIdx);
            break;
        }
        case OptRitz::SM: {
            optVal = evals.cwiseAbs().minCoeff(&optIdx);
            break;
        }
        case OptRitz::LR: {
            optVal = evals.maxCoeff(&optIdx);
            break;
        }
        case OptRitz::LM: {
            optVal = evals.cwiseAbs().maxCoeff(&optIdx);
            break;
        }
        default: {
            // Take the closest to the current GsiH
            optVal = (evals.array() - std::real(v0H1v0 / v0H2v0)).cwiseAbs().minCoeff(&optIdx);
            optVal = evals.coeff(optIdx);
            break;
        }
    }

    Eigen::VectorXd col = evecs.col(optIdx).normalized().real();
    if(col.hasNaN()) {
        res.mixed_blk = t0;
        res.alpha_mps = 1.0;
        res.alpha_ene = 0.0;
        res.alpha_var = 0.0;
    } else {
        res.mixed_blk = Eigen::Tensor<cplx, 3>(t0.dimensions());
        res.alpha_mps = col.coeff(0);
        res.alpha_ene = col.coeff(1);
        res.alpha_var = col.coeff(2);
        auto vf       = tenx::VectorMap(res.mixed_blk);
        vf            = (col.coeff(0) * v0 + col.coeff(1) * v1 + col.coeff(2) * v2).normalized();
    }
    tools::log->debug("mixed state result: <H>/<H²> = {:.16f}", optVal);

    return res;
}

EnvExpansionResult tools::finite::env::get_optimally_mixed_block(const std::vector<size_t> &sites, const StateFinite &state, const ModelFinite &model,
                                                                 const EdgesFinite &edges, OptAlgo algo, OptRitz ritz) {
    switch(algo) {
        case OptAlgo::DMRG: {
            return internal::get_optimally_mixed_block_H1(sites, state, model, edges, ritz);
        }
        case OptAlgo::DMRGX: {
            return internal::get_optimally_mixed_block_VarH(sites, state, model, edges, ritz);
        }
        case OptAlgo::HYBRID_DMRGX: {
            return internal::get_optimally_mixed_block_VarH(sites, state, model, edges, ritz);
        }
        case OptAlgo::XDMRG: {
            return internal::get_optimally_mixed_block_H2(sites, state, model, edges, ritz);
        }
        case OptAlgo::GDMRG: {
            return internal::get_optimally_mixed_block_GsiH(sites, state, model, edges, ritz);
        }
        default: return {};
    }
}

void merge_expansion_term_PR(const StateFinite &state, MpsSite &mpsL, const Eigen::Tensor<cplx, 3> &ML_PL, MpsSite &mpsR, const Eigen::Tensor<cplx, 3> &MR_PR,
                             const svd::config &svd_cfg) {
    // The expanded bond sits between mpsL and mpsR
    // We zero-pad mpsL and enrich mpsR.
    // During forward expansion -->
    //      * mpsL is AC(i), or B(i) during multisite dmrg
    //      * mpsR is B(i+1)
    // During backward expansion <--
    //      * mpsL is A(i-1) always
    //      * mpsR is AC(i) always
    // Thus, the possible situations are [AC,B] or [B,B] or [A,AC]
    // After USV = SVD(MR_PR):
    //      If [AC,B]:
    //           * ML_PL is [AC(i), PL]            <--- Note that we use full AC(i)! Not bare A(i)
    //           * MR_PR is [B(i+1), PR]^T
    //           * mpsL:  AC(i)   = ML_PL * U     <---- lost left normalization, but it is not needed during optimization
    //           * mpsL:  C(i)   = S
    //           * mpsR:  B(i+1) = V
    //      If [B,B]:
    //           * ML_PL is [B(i), PL]
    //           * MR_PR is [B(i+1), PR]^T
    //           * mpsL:  B(i)   = ML_PL * U * S  <-- lost right normalization, but it is not needed during optimization
    //           * mpsL:  Λ(i)   = S
    //           * mpsR:  B(i+1) = V
    //      If [A,AC]:
    //           * ML_PL is [A(i-1), PL]
    //           * MR_PR is [A(i), PR]^T             <--- Note that we use bare A(i)! Not AC(i)
    //           * mpsL:  A(i-1) = ML_PL * U         <--- lost left normalization, but it will not be needed during the optimization later
    //           * mpsR:  Λ(i) = S, (will become a center later)
    //           * mpsR:  A(i) = S * V (note that we replace bare A(i), not AC(i)!)
    //           * The old C(i) living in mpsR is outdated, but it incorporated into AC(i) when building MR_PR.
    //           * We can simply replace C(i) with an identity matrix. It will not be needed after we transform AC(i) into a B(i) later.
    //

    svd::solver svd;
    auto        posL = mpsL.get_position();
    auto        labL = mpsL.get_label();
    auto        labR = mpsR.get_label();
    auto [U, S, V]   = svd.schmidt_into_right_normalized(MR_PR, mpsR.spin_dim(), svd_cfg);

    if(labL == "AC" and labR == "B") {
        mpsR.set_M(V);
        mpsR.stash_U(U, posL);
        mpsR.stash_C(S, -1.0, posL); // Set a negative truncation error to ignore it.
        mpsL.set_M(ML_PL);
        mpsL.take_stash(mpsR); // normalization of mpsL is lost here.
    } else if(labL == "B" and labR == "B") {
        auto US = Eigen::Tensor<cplx, 3>(U.contract(tenx::asDiagonal(S), tenx::idx({2}, {0})));
        mpsR.set_M(V);
        mpsR.stash_U(US, posL);
        mpsR.stash_S(S, -1.0, posL); // Set a negative truncation error to ignore it.
        mpsL.set_M(ML_PL);
        mpsL.take_stash(mpsR); // normalization of mpsL is lost here.
    } else if(labL == "A" and labR == "AC") {
        auto SV = Eigen::Tensor<cplx, 3>(tenx::asDiagonal(S).contract(V, tenx::idx({1}, {1})).shuffle(tenx::array3{1, 0, 2}));
        auto C1 = Eigen::Tensor<cplx, 1>(V.dimension(2));
        C1.setConstant(cplx(1.0, 0.0)); // Replaces C with an identity
        mpsR.set_M(SV);
        mpsR.set_L(S);
        mpsR.set_LC(C1);
        mpsR.stash_U(U, posL);
        mpsL.set_M(ML_PL);
        mpsL.take_stash(mpsR);
    } else {
        throw except::runtime_error("merge_expansion_term_PR: could not match case: [{},{}]", labL, labR);
    }

    {
        // Make mpsL normalized so that later checks can succeed
        auto                   multisite_mpsL = state.get_multisite_mps({posL});
        cplx                   norm_old       = tools::common::contraction::contract_mps_norm(multisite_mpsL);
        Eigen::Tensor<cplx, 3> M_tmp          = mpsL.get_M_bare() * mpsL.get_M_bare().constant(std::pow(norm_old, -0.5)); // Rescale
        mpsL.set_M(M_tmp);
        if constexpr(settings::debug_expansion) {
            auto mpsL_final = state.get_multisite_mps({mpsL.get_position()});
            cplx norm_new   = tools::common::contraction::contract_mps_norm(mpsL_final);
            tools::log->debug("Normalized expanded mps {}({}): {:.16f} -> {:.16f}", mpsL.get_label(), mpsL.get_position(), std::abs(norm_old),
                              std::abs(norm_new));
        }
    }
}

void merge_expansion_term_PL(const StateFinite &state, MpsSite &mpsL, const Eigen::Tensor<cplx, 3> &ML_PL, MpsSite &mpsR, const Eigen::Tensor<cplx, 3> &MR_PR,
                             const svd::config &svd_cfg) {
    // The expanded bond sits between mpsL and mpsR.
    // During forward expansion <--
    //      * mpsL is A(i-1) always
    //      * mpsR is AC(i), or A(i) during multisite dmrg
    // During backward expansion -->
    //      * mpsL is AC(i) always
    //      * mpsR is B(i+1) always
    // Thus, the possible situations are  [A, AC], [A,A] or [AC,B]
    // After USV = SVD(ML_PL):
    //      If [A, AC]:
    //           * ML_PL is [A(i-1), PL]
    //           * MR_PR is [A(i), PR]^T             <--- Note that we use bare A(i)! Not AC(i)
    //           * mpsL:  A(i-1) = U
    //           * mpsR:  Λ(i)   = S                <--- takes stash S
    //           * mpsR:  A(i)   = S * V * MR_PR    <--- takes stash S,V and loses left-right normalization
    //           * mpsR:  C(i)                      <--- does not change
    //      If [A,A]:
    //           * ML_PL is [A(i-1), PL]^T
    //           * MR_PR is [A(i), 0]^T             <--- Note that we use bare A(i)! Not AC(i)
    //           * mpsL:  A(i-1) = U
    //           * mpsR:  Λ(i)   = S (takes stash S)
    //           * mpsR:  A(i) = S * V * MR_PR (takes stash S,SV and loses left normalization)
    //      If [AC,B]:
    //           * ML_PL is [AC(i), PL]^T
    //           * MR_PR is [B(i+1), 0]^T
    //           * mpsL:  A(i) = U
    //           * mpsL:  C(i) = S
    //           * mpsR:  B(i+1) = V * MR_PR (loses right normalization, but that is not be needed during the next optimization)
    //

    svd::solver svd;
    auto        posR = mpsR.get_position();
    auto        labL = mpsL.get_label();
    auto        labR = mpsR.get_label();
    auto [U, S, V]   = svd.schmidt_into_left_normalized(ML_PL, mpsL.spin_dim(), svd_cfg);

    if(labL == "A" and labR == "AC") {
        auto SV = Eigen::Tensor<cplx, 3>(tenx::asDiagonal(S).contract(V, tenx::idx({1}, {1})).shuffle(tenx::array3{1, 0, 2}));
        mpsL.set_M(U);
        mpsL.stash_S(S, -1.0, posR); // Set a negative truncation error to ignore it.
        mpsL.stash_V(SV, posR);
        mpsR.set_M(MR_PR);
        mpsR.take_stash(mpsL); // normalization of mpsR is lost here
    } else if(labL == "A" and labR == "A") {
        auto SV = Eigen::Tensor<cplx, 3>(tenx::asDiagonal(S).contract(V, tenx::idx({1}, {1})).shuffle(tenx::array3{1, 0, 2}));
        mpsL.set_M(U);
        mpsL.stash_S(S, -1.0, posR); // Set a negative truncation error to ignore it.
        mpsL.stash_V(SV, posR);
        mpsR.set_M(MR_PR);
        mpsR.take_stash(mpsL); // normalization of mpsR is lost here
    } else if(labL == "AC" and labR == "B") {
        mpsL.set_M(U);
        mpsL.set_LC(S);
        mpsL.stash_V(V, posR);
        mpsR.set_M(MR_PR);
        mpsR.take_stash(mpsL);
    } else {
        throw except::runtime_error("merge_expansion_term_PL: could not match case: [{},{}]", labL, labR);
    }

    {
        // Make mpsR normalized so that later checks can succeed
        auto                   multisite_mpsR = state.get_multisite_mps({posR});
        cplx                   norm_old       = tools::common::contraction::contract_mps_norm(multisite_mpsR);
        Eigen::Tensor<cplx, 3> M_tmp          = mpsR.get_M_bare() * mpsR.get_M_bare().constant(std::pow(norm_old, -0.5)); // Rescale by the norm
        mpsR.set_M(M_tmp);
        if constexpr(settings::debug_expansion) {
            auto mpsR_final = state.get_multisite_mps({mpsR.get_position()});
            cplx norm_new   = tools::common::contraction::contract_mps_norm(mpsR_final);
            tools::log->debug("Normalized expanded mps {}({}): {:.16f} -> {:.16f}", mpsR.get_label(), mpsR.get_position(), std::abs(norm_old),
                              std::abs(norm_new));
        }
    }
}

void merge_expanded_terms(const StateFinite &state, MpsSite &mpsL, MpsSite &mpsR, const Eigen::Tensor<cplx, 4> &MP, const svd::config &svd_cfg) {
    // The expanded bond sits between mpsL and mpsR.
    // The possible situations are  [A, AC], [A,A], [AC,B] or [B,B]
    // After USV = SVD(MP), where MP = ML_PL * MR_PR
    //      If [A, AC]:
    //           * ML_PL is [A(i-1), PL]
    //           * MR_PR is [AC(i), PR]^T           <--- Note that we use AC(i)! Not bare A(i)
    //           * mpsL:  A(i-1) = U
    //           * mpsR:  Λ(i)   = S                <--- takes stash S
    //           * mpsR:  AC(i)  = S * V            <--- takes stash S,V and loses left normalization
    //           * mpsR:  C(i)   = I                <--- does not change, however, it is already baked into AC, so we can set this to identity
    //      If [A,A]:
    //           * ML_PL is [A(i-1), PL]
    //           * MR_PR is [A(i), PR]^T             <--- Note that we use bare A(i)! Not AC(i)
    //           * mpsL:  A(i-1) = U
    //           * mpsR:  Λ(i)   = S                 <--- takes stash S
    //           * mpsR:  A(i) = S * V               <--- takes stash S,V and loses left normalization
    //      If [AC,B]:
    //           * ML_PL is [AC(i), PL]
    //           * MR_PR is [B(i+1), PR]^T
    //           * mpsL:  A(i) = U
    //           * mpsL:  C(i) = S
    //           * mpsR:  B(i+1) = V
    //      If [B,B]:
    //           * ML_PL is [B(i), PL]
    //           * MR_PR is [B(i+1), PR]^T
    //           * mpsL:  B(i)   = U * S  <-- loses right normalization, but it is not needed during optimization
    //           * mpsL:  Λ(i)   = S
    //           * mpsR:  B(i+1) = V
    //

    svd::solver svd;
    auto        posL = mpsL.get_position();
    auto        posR = mpsR.get_position();
    auto        labL = mpsL.get_label();
    auto        labR = mpsR.get_label();
    auto [U, S, V]   = svd.schmidt(MP, svd_cfg);

    if(labL == "A" and labR == "AC") {
        auto SV = Eigen::Tensor<cplx, 3>(tenx::asDiagonal(S).contract(V, tenx::idx({1}, {1})).shuffle(tenx::array3{1, 0, 2}));
        auto C1 = Eigen::Tensor<cplx, 1>(V.dimension(2));
        C1.setConstant(cplx(1.0, 0.0)); // Replaces C with an identity
        mpsL.set_M(U);
        mpsL.stash_S(S, -1.0, posR); // Set a negative truncation error to ignore it.
        mpsR.set_M(SV);              // normalization of mpsR is lost here
        mpsR.take_stash(mpsL);
        mpsR.set_LC(C1);
    } else if(labL == "A" and labR == "A") {
        auto SV = Eigen::Tensor<cplx, 3>(tenx::asDiagonal(S).contract(V, tenx::idx({1}, {1})).shuffle(tenx::array3{1, 0, 2}));
        mpsL.set_M(U);
        mpsR.set_M(SV); // normalization of mpsR is lost here
        mpsR.set_L(S, -1.0);
    } else if(labL == "AC" and labR == "B") {
        mpsL.set_M(U);
        mpsL.set_LC(S);
        mpsR.set_M(V);
    } else if(labL == "B" and labR == "B") {
        auto US = Eigen::Tensor<cplx, 3>(U.contract(tenx::asDiagonal(S), tenx::idx({2}, {0})));
        mpsL.set_M(US);
        mpsL.set_L(S, -1.0);
        mpsR.set_M(V);

    } else {
        throw except::runtime_error("merge_expansion_term_PL: could not match case: [{},{}]", labL, labR);
    }

    {
        // Make mpsR normalized so that later checks can succeed
        auto                   multisite_mpsR = state.get_multisite_mps({posR});
        cplx                   norm_old       = tools::common::contraction::contract_mps_norm(multisite_mpsR);
        Eigen::Tensor<cplx, 3> M_tmp          = mpsR.get_M_bare() * mpsR.get_M_bare().constant(std::pow(norm_old, -0.5)); // Rescale by the norm
        mpsR.set_M(M_tmp);
        if constexpr(settings::debug or settings::debug_expansion) {
            auto mpsR_final = state.get_multisite_mps({mpsR.get_position()});
            cplx norm_new   = tools::common::contraction::contract_mps_norm(mpsR_final);
            tools::log->debug("Normalized expanded mps {}({}): {:.16f} -> {:.16f}", mpsR.get_label(), mpsR.get_position(), std::abs(norm_old),
                              std::abs(norm_new));
        }
    }
}

/*!
    Follows the environment (aka subspace) expansion technique as explained in https://link.aps.org/doi/10.1103/PhysRevB.91.155115
    Note that in this convention we expand/enrich the trailing bond during a DMRG sweep.
    Going left-to-right with active site == i, we update:

        - AC(i)   --> SVD( [AC(i)  , PL(i)] ) = U*S*V, to update AC(i) = U, LC = S and V is contracted onto B(i+1) below
        - B(i+1)  --> V * [B(i+1) , 0   ]^T  : This site loses right-normalization when we contract V from the SVD above.

    where PL(i) = alpha * ENVL(i) * AC(i) * MPO(i) (dimensions d, AC.dimension(2), AC.dimension(0)*MPO(i).dimension(1))
    Then, immediately after calling this function, we should move the current center site i --> i+1

    Similarly, going right-to-left with active site == i, we update:

        - AC(i)   --> SVD( [AC(i) , PR(i) ]^T ) = U*S*V, to update AC(i)=S*V, and U is contracted onto A(i-1) below
        - A(i-1)  -->      [A(i-1), 0    ] * U  : This one loses left-normalization when we contract U from the SVD above.

    where PR(i) = alpha * ENVR(i) * AC(i) * MPO(i) (dimensions d, AC.dimension(1), AC.dimension(2) * MPO(i).dimension(0))
    Then, immediately after calling this function, we should move the current center site i --> i-1
*/
EnvExpansionResult tools::finite::env::expand_environment_backward(StateFinite &state, const ModelFinite &model, EdgesFinite &edges,
                                                                   EnvExpandMode envExpandMode, svd::config svd_cfg) {
    if(state.get_position<long>() < 0) { return {}; }

    auto                pos = state.get_position<size_t>();
    std::vector<size_t> pos_expanded;
    if(state.get_direction() > 0 and pos == std::clamp<size_t>(pos, 0, state.get_length<size_t>() - 2)) pos_expanded = {pos, pos + 1};
    if(state.get_direction() < 0 and pos == std::clamp<size_t>(pos, 1, state.get_length<size_t>() - 1)) pos_expanded = {pos - 1, pos};
    if(pos_expanded.empty()) {
        return {}; // No expansion
    }

    // Define the left and right mps that will get modified
    state.clear_cache();
    state.clear_measurements();

    auto  posL     = pos_expanded.front();
    auto  posR     = pos_expanded.back();
    auto &mpsL     = state.get_mps_site(posL);
    auto &mpsR     = state.get_mps_site(posR);
    auto  dimL_old = mpsL.dimensions();
    auto  dimR_old = mpsR.dimensions();

    // Save expansion metadata
    auto res     = EnvExpansionResult();
    res.posL     = posL;
    res.posR     = posR;
    res.dimL_old = dimL_old;
    res.dimR_old = dimR_old;

    // Stop here if it's just a query
    if(envExpandMode == EnvExpandMode::NONE) return res; // Position query
    res.ene_old                                           = tools::finite::measure::energy(state, model, edges);
    res.var_old                                           = tools::finite::measure::energy_variance(state, model, edges);
    std::tie(res.alpha_mps, res.alpha_ene, res.alpha_var) = get_optimal_mixing_factors(pos_expanded, state, model, edges, envExpandMode);
    if(res.alpha_ene == 0 and res.alpha_var == 0) {
        tools::log->debug("Expansion canceled: {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(),
                          mpsR.get_position(), res.alpha_ene, res.alpha_var);
        return res;
    }
    tools::log->debug("Expanding {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(), mpsR.get_position(),
                      res.alpha_ene, res.alpha_var);

    // Set up the SVD
    svd_cfg.truncation_limit = svd_cfg.truncation_limit.value_or(settings::precision::svd_truncation_min);
    svd_cfg.rank_max         = svd_cfg.rank_max.value_or(settings::get_bond_max(state.get_algorithm()));
    long        bond_lim     = svd_cfg.rank_max.value();
    svd::solver svd(svd_cfg);

    if(state.get_direction() > 0 and posL > 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is left-to-right:
        //      * mpsL is "AC(i)" and will get moved into envL later
        //      * mpsR is "B(i+1)" and will become the active site after moving center position
        auto &mpoL = model.get_mpo(posL);
        // auto                  &mpoR   = model.get_mpo(posR);
        Eigen::Tensor<cplx, 3> ML_PL = mpsL.get_M(); // mpsL is going to be optimized, enriched with PL
        // Eigen::Tensor<cplx, 3> MR_PR  = mpsR.get_M(); // mpsL is going to be optimized, enriched with PL
        long PLdim2 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene > 0) {
            auto &envL = edges.get_env_eneL(posL);
            // auto                  &envR      = edges.get_env_eneR(posR);
            long                   chi_max = 2 * bond_lim / (1 + mpoL.MPO().dimension(1));
            Eigen::Tensor<cplx, 3> PL      = envL.get_expansion_term(mpsL, mpoL, res.alpha_ene, chi_max);
            // Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, 0.0, -1);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            // Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL = std::move(ML_PL_tmp);
            // MR_PR                            = std::move(MR_PR_tmp);
            PLdim2 += PL.dimension(2);
            res.env_ene = envL; // Saves a reference to the expanded environment
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var > 0) {
            auto &envL = edges.get_env_varL(posL);
            // auto                  &envR      = edges.get_env_varR(posR);
            long                   chi_max = 2 * bond_lim / (1 + mpoL.MPO2().dimension(1));
            Eigen::Tensor<cplx, 3> PL      = envL.get_expansion_term(mpsL, mpoL, res.alpha_var, chi_max);
            // Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, 0.0, -1);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            // Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL = std::move(ML_PL_tmp);
            // MR_PR                            = std::move(MR_PR_tmp);
            PLdim2 += PL.dimension(2);
            res.env_var = envL; // Saves a reference to the expanded environment
        }
        if(PLdim2 == 0) { tools::log->error("expand_environment_backward: PLdim2 == 0: this is likely an error"); }
        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsR.spin_dim(), PLdim2, mpsR.get_chiR());
        Eigen::Tensor<cplx, 3> MR_P0 = mpsR.get_M().concatenate(P0, 1); // mpsR is going into the environment, padded with zeros
        merge_expansion_term_PL(state, mpsL, ML_PL, mpsR, MR_P0, svd_cfg);
        tools::log->debug("Environment expansion backward pos {} | {} | αₑ:{:.2e} αᵥ:{:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), dimL_old[2],
                          ML_PL.dimension(2), mpsL.get_chiR());
        if(dimL_old[1] != mpsL.get_chiL()) throw except::runtime_error("mpsL changed chiL during environment expansion: {} -> {}", dimL_old, mpsL.dimensions());
        if(dimR_old[2] != mpsR.get_chiR()) throw except::runtime_error("mpsR changed chiR during environment expansion: {} -> {}", dimR_old, mpsR.dimensions());
    }
    if(state.get_direction() < 0) {
        // The expanded bond sits between mpsL and mpsR. In 1-site DMRG, when the direction is left-to-right:
        //      * mpsL is "A(i-1)" and will become the active site after moving center position
        //      * mpsR is "AC(i)" and is the active site that will get moved into envR later
        // auto                  &mpoL   = model.get_mpo(posL);
        auto &mpoR = model.get_mpo(posR);
        // Eigen::Tensor<cplx, 3> ML_PL  = mpsL.get_M(); // [ AC  P ]^T  including LC
        Eigen::Tensor<cplx, 3> MR_PR  = mpsR.get_M(); // [ AC  P ]^T  including LC
        long                   PRdim1 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene > 0) {
            // auto                  &envL       = edges.get_env_eneL(posL);
            auto &envR    = edges.get_env_eneR(posR);
            long  chi_max = 2 * bond_lim / (1 + mpoR.MPO().dimension(0));
            // Eigen::Tensor<cplx, 3> PL         = envL.get_expansion_term(mpsL, mpoL, 0.0, chi_max);
            Eigen::Tensor<cplx, 3> PR = envR.get_expansion_term(mpsR, mpoR, res.alpha_ene, chi_max);
            // Eigen::Tensor<cplx, 3> ML_PL_temp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_temp = MR_PR.concatenate(PR, 1);
            // ML_PL                             = std::move(ML_PL_temp);
            MR_PR = std::move(MR_PR_temp);
            PRdim1 += PR.dimension(1);
            res.env_ene = envR;
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var > 0) {
            // auto                  &envL       = edges.get_env_varL(posL);
            auto &envR    = edges.get_env_varR(posR);
            long  chi_max = 2 * bond_lim / (1 + mpoR.MPO2().dimension(0));
            // Eigen::Tensor<cplx, 3> PL         = envL.get_expansion_term(mpsL, mpoL, 0.0, chi_max);
            Eigen::Tensor<cplx, 3> PR = envR.get_expansion_term(mpsR, mpoR, res.alpha_var, chi_max);
            // Eigen::Tensor<cplx, 3> ML_PL_temp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_temp = MR_PR.concatenate(PR, 1);
            // ML_PL                             = std::move(ML_PL_temp);
            MR_PR = std::move(MR_PR_temp);
            PRdim1 += PR.dimension(1);
            res.env_var = envR;
        }
        if(PRdim1 == 0) { tools::log->error("expand_environment_backward: PRdim1 == 0: this is likely an error"); }
        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsL.spin_dim(), mpsL.get_chiL(), PRdim1);
        Eigen::Tensor<cplx, 3> ML_P0 = mpsL.get_M().concatenate(P0, 2); // [ A   0 ]
        merge_expansion_term_PR(state, mpsL, ML_P0, mpsR, MR_PR, svd_cfg);

        tools::log->debug("Environment expansion backward pos {} | {} | αₑ:{:.2e} αᵥ:{:.3e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), dimL_old[2],
                          MR_PR.dimension(1), mpsL.get_chiR());
    }
    if constexpr(settings::debug_expansion) mpsL.assert_normalized();
    if constexpr(settings::debug_expansion) mpsR.assert_normalized();
    env::rebuild_edges(state, model, edges);
    res.dimL_new = mpsL.dimensions();
    res.dimR_new = mpsR.dimensions();
    res.ene_new  = tools::finite::measure::energy(state, model, edges);
    res.var_new  = tools::finite::measure::energy_variance(state, model, edges);
    res.ok       = true;
    return res;
}

EnvExpansionResult tools::finite::env::expand_environment_forward_nsite(StateFinite &state, const ModelFinite &model, EdgesFinite &edges, OptAlgo algo,
                                                                        OptRitz ritz, svd::config svd_cfg) {
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("expand_environment_forward_nsite: All lengths not equal: state {} | model {} | edges {}", state.get_length(),
                                    model.get_length(), edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("expand_environment_forward_nsite: All active sites are not equal: state {} | model {} | edges {}", state.active_sites,
                                    model.active_sites, edges.active_sites);
    if(state.active_sites.empty()) throw except::logic_error("No active sites for environment expansion");

    // Determine which bond to expand
    std::vector<size_t> pos_active_and_expanded;
    if(state.get_direction() > 0) {
        auto pos = state.active_sites.back();
        // if(pos == std::clamp<size_t>(pos, 0, state.get_length<size_t>() - 2)) pos_expanded = {pos, pos + 1};
        if(pos != std::clamp<size_t>(pos, 0, state.get_length<size_t>() - 2)) {
            tools::log->debug("Expansion canceled: No sites to expand");
            return {};
        }
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), state.active_sites.begin(), state.active_sites.end());
        pos_active_and_expanded.emplace_back(pos + 1);
    }
    if(state.get_direction() < 0) {
        auto pos = state.active_sites.front();
        // if(pos == std::clamp<size_t>(pos, 1, state.get_length<size_t>() - 1)) pos_expanded = {pos - 1, pos};
        if(pos != std::clamp<size_t>(pos, 1, state.get_length<size_t>() - 1)) {
            tools::log->debug("Expansion canceled: No sites to expand");
            return {};
        }
        pos_active_and_expanded.emplace_back(pos - 1);
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), state.active_sites.begin(), state.active_sites.end());
    }

    // Define the left and right mps that will get modified
    auto res = get_optimally_mixed_block(pos_active_and_expanded, state, model, edges, algo, ritz);

    const auto &mpsL = state.get_mps_site(res.posL);
    const auto &mpsR = state.get_mps_site(res.posR);
    res.ene_old      = tools::finite::measure::energy(state, model, edges);
    res.var_old      = tools::finite::measure::energy_variance(state, model, edges);

    if(res.alpha_ene == 0 and res.alpha_var == 0) {
        tools::log->debug("Expansion canceled: {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(),
                          mpsR.get_position(), res.alpha_ene, res.alpha_var);
        return res;
    }
    tools::log->debug("Expanding {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(), mpsR.get_position(),
                      res.alpha_ene, res.alpha_var);

    mps::merge_multisite_mps(state, res.mixed_blk, pos_active_and_expanded, state.get_position<long>(), MergeEvent::EXP, svd_cfg);

    res.dims_new = state.get_mps_dims(pos_active_and_expanded);
    res.bond_new = state.get_bond_dims(pos_active_and_expanded);

    tools::log->debug("Environment expansion forward pos {} | {} {} | αₑ:{:.2e} αᵥ:{:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {}", pos_active_and_expanded,
                      enum2sv(algo), enum2sv(ritz), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), res.bond_old,
                      res.bond_new);
    state.clear_cache();
    state.clear_measurements();
    for(const auto &mps : state.get_mps(pos_active_and_expanded)) mps.get().assert_normalized();
    env::rebuild_edges(state, model, edges);

    res.dimL_new = mpsL.dimensions();
    res.dimR_new = mpsR.dimensions();
    res.ene_new  = tools::finite::measure::energy(state, model, edges);
    res.var_new  = tools::finite::measure::energy_variance(state, model, edges);
    res.ok       = true;
    return res;
}

/*!
    This is a modified (aka subspace) expansion technique, compared to the one explained in https://link.aps.org/doi/10.1103/PhysRevB.91.155115
    In this convention we expand/enrich the forward bond during a DMRG sweep, i.e. the bond one step ahead in the sweep direction.

    In 1-site DMRG going left-to-right with active site == i, we update:

        - AC(i), LC(i) -->      [AC(i)     ,  0 ]  * U, S(i)      : AC(i) is bare, without LC(i) (which is used below). Loses left-normalization due to U.
        - LC*B(i+1)    --> SVD( [LC*B(i+1) , PR(i+1) ]^T ) = U*S*V: update B(i+1) = V, LC(i) = S and U is contracted onto AC(i) above

    where PR(i+1) = alpha * ENVR(i+1) * B(i+1) * MPO(i+1) (dimensions d(i-1), B(i+1).dimension(1), B(i+1).dimension(2) * MPO(i+1).dimension(0))
    Then, immediately after calling this function, we should optimize the current site AC(i) with a DMRG step.

    Similarly, in 1-site DMRG going right-to-left with active site == i, we update:

         - A(i-1)       --> SVD( [A(i-1), PL(i-1)] ) = U*S*V : update A(i-1) = U, L(i) = S and S*V is contracted onto A(i) below.
         - L(i), A(i)   --> S, S*V*[ A(i)    ,  0      ]^T     : A(i) (which is AC(i) bare), loses left-normalization due to SV.

    where PL(i-1) = alpha * ENVL(i-1) * A(i-1) * MPO(i-1) (dimensions d(i-1), A(i-1).dimension(2), A(i-1).dimension(0)*MPO(i-1).dimension(1))
    Then, immediately after calling this function, we should optimize the current site AC(i) with a DMRG step.


    Thus, one step of the DMRG algorithm proceeds as

    1. expand between sites i, i+1
    2. optimize $\Psi^{i} = A_C^{i} \Lambda_C^{i}$.
    3. split $\Psi^{i} \stackrel{\text{SVD}}{\rightarrow}A_C^{i}  \Lambda_C^{i} V^\dagger B^{i+1}$ normally and update the MPS on sites $i,i+1$.
    4. move: $A_C^i \Lambda_C^i \rightarrow A^i$ and $\Lambda_C^i B^{i+1} \rightarrow \Lambda^i A_C^{i+1} \Lambda_C^{i+1}$, and repeat from step 1.
    5. update $\alpha$: decrease by $0.01/L$ if **if the variance decreased** in the last step, else increase by $10/L$, bounded by $\alpha \in [\alpha_{\min},
   \alpha_{\max}]$, where:

      - $\alpha_{\min} = \text{Var}(H) \cdot 10^{-3}$,
      - $\alpha_{\max}= \min{(10^{-4}, \text{Var}(H))}$.

    The main differences compared with the original:
    - The SVD in the expansion can be done with high precision and bond dimension: the real truncation happens after the optimization in step 3.
    - Therefore, the optimization step sees a much richer environment, which speeds up convergence.
    - Backwards expansion pushes a non-optimal/degraded site-tensor into the trailing environment. Forward expansion optimizes it first.

    Note that this expansion works similarly for multisite dmrg, with the expanded site chosen accordingly.


*/
EnvExpansionResult tools::finite::env::expand_environment_forward_1site(StateFinite &state, const ModelFinite &model, EdgesFinite &edges,
                                                                        EnvExpandMode envExpandMode, svd::config svd_cfg) {
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("expand_environment_forward_1site: All lengths not equal: state {} | model {} | edges {}", state.get_length(),
                                    model.get_length(), edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("expand_environment_forward_1site: All active sites are not equal: state {} | model {} | edges {}", state.active_sites,
                                    model.active_sites, edges.active_sites);
    if(state.active_sites.empty()) throw except::logic_error("No active sites for environment expansion");

    // Determine which bond to expand
    std::vector<size_t> pos_expanded;
    std::vector<size_t> pos_active_and_expanded;
    if(state.get_direction() > 0) {
        auto pos = state.active_sites.back();
        if(pos == std::clamp<size_t>(pos, 0, state.get_length<size_t>() - 2)) pos_expanded = {pos, pos + 1};
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), state.active_sites.begin(), state.active_sites.end());
        pos_active_and_expanded.emplace_back(pos + 1);
    }
    if(state.get_direction() < 0) {
        auto pos = state.active_sites.front();
        if(pos == std::clamp<size_t>(pos, 1, state.get_length<size_t>() - 1)) pos_expanded = {pos - 1, pos};
        pos_active_and_expanded.emplace_back(pos - 1);
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), state.active_sites.begin(), state.active_sites.end());
    }
    if(pos_expanded.empty()) {
        return {}; // No update
    }

    // Define the left and right mps that will get modified
    state.clear_cache();
    state.clear_measurements();
    auto  posL     = pos_expanded.front();
    auto  posR     = pos_expanded.back();
    auto &mpsL     = state.get_mps_site(posL);
    auto &mpsR     = state.get_mps_site(posR);
    auto  dimL_old = mpsL.dimensions();
    auto  dimR_old = mpsR.dimensions();

    // Save expansion metadata
    auto res     = EnvExpansionResult();
    res.posL     = posL;
    res.posR     = posR;
    res.dimL_old = dimL_old;
    res.dimR_old = dimR_old;

    // Stop here if it's just a query
    if(envExpandMode == EnvExpandMode::NONE) return res; // Position query
    res.ene_old = tools::finite::measure::energy(state, model, edges);
    res.var_old = tools::finite::measure::energy_variance(state, model, edges);

    // res.alpha_ene = get_optimal_mixing_factor_ene_new(pos_active_and_expanded, state, model, edges, envExpandMode);
    // res.alpha_var = get_optimal_mixing_factor_var_new(pos_active_and_expanded, state, model, edges, envExpandMode);
    std::tie(res.alpha_mps, res.alpha_ene, res.alpha_var) = get_optimal_mixing_factors(pos_expanded, state, model, edges, envExpandMode);
    if(res.alpha_ene == 0 and res.alpha_var == 0) {
        tools::log->debug("Expansion canceled: {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(),
                          mpsR.get_position(), res.alpha_ene, res.alpha_var);
        return res;
    }
    tools::log->debug("Expanding {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(), mpsR.get_position(),
                      res.alpha_ene, res.alpha_var);

    // Set up the SVD
    // Bond dimension can't grow faster than x spin_dim, but we can generate a highly enriched environment here for optimization,
    // and let the proper truncation happen after optimization instead.
    auto bond_lim            = std::min(mpsL.spin_dim() * mpsL.get_chiL(), mpsR.spin_dim() * mpsR.get_chiR());
    svd_cfg.truncation_limit = svd_cfg.truncation_limit.value_or(settings::precision::svd_truncation_min);
    svd_cfg.rank_max         = std::min(bond_lim, svd_cfg.rank_max.value_or(bond_lim));
    // svd_cfg.svd_rtn          = svd::rtn::gersvd;
    // svd_cfg.loglevel         = 0;
    svd::solver svd;

    if(state.get_direction() > 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is left-to-right:
        // We have
        //      * mpsL is AC(i), or B(i) during multisite dmrg
        //      * mpsR is B(i+1) and belongs in envR later in the optimization step
        auto                  &mpoL   = model.get_mpo(posL);
        auto                  &mpoR   = model.get_mpo(posR);
        Eigen::Tensor<cplx, 3> ML_PL  = mpsL.get_M();                                                   // mpsR is going into the environment, enriched with PR.
        Eigen::Tensor<cplx, 3> MR_PR  = mpsR.get_M() * mpsR.get_M().constant(cplx(res.alpha_mps, 0.0)); // mpsR is going into the environment, enriched with PR.
        long                   PRdim1 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene != 0) {
            auto                  &envL      = edges.get_env_eneL(posL);
            auto                  &envR      = edges.get_env_eneR(posR);
            long                   chi_max   = 2 * bond_lim / (1 + mpoR.MPO().dimension(0));
            Eigen::Tensor<cplx, 3> PL        = envL.get_expansion_term(mpsL, mpoL, std::sqrt(std::abs(res.alpha_ene)), -1);
            Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, std::sqrt(std::abs(res.alpha_ene)), -1);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL                            = std::move(ML_PL_tmp);
            MR_PR                            = std::move(MR_PR_tmp);
            PRdim1 += PR.dimension(1);
            res.env_ene = envR;
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var != 0) {
            auto                  &envL      = edges.get_env_varL(posL);
            auto                  &envR      = edges.get_env_varR(posR);
            long                   chi_max   = 2 * bond_lim / (1 + mpoR.MPO2().dimension(0));
            Eigen::Tensor<cplx, 3> PL        = envL.get_expansion_term(mpsL, mpoL, std::sqrt(std::abs(res.alpha_var)), -1);
            Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, std::sqrt(std::abs(res.alpha_var)), -1);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL                            = std::move(ML_PL_tmp);
            MR_PR                            = std::move(MR_PR_tmp);
            PRdim1 += PR.dimension(1);
            res.env_var = envR;
        }
        // Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsL.spin_dim(), mpsL.get_chiL(), PRdim1);
        // Eigen::Tensor<cplx, 3> ML_P0 = mpsL.get_M().concatenate(P0, 2); // mpsL is going to be optimized, zero-padded with P0
        // merge_expansion_term_PR(state, mpsL, ML_PL, mpsR, MR_PR, svd_cfg);
        Eigen::Tensor<cplx, 4> MP = ML_PL.contract(MR_PR, tenx::idx({2}, {1}));
        merge_expanded_terms(state, mpsL, mpsR, MP, svd_cfg);

        state.clear_measurements();
        state.clear_cache();

        tools::log->debug("Environment expansion forward pos {} | {} | αₑ:{:.2e} αᵥ:{:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), dimL_old[2],
                          MR_PR.dimension(1), mpsL.get_chiR());
    }
    if(state.get_direction() < 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is right-to-left:
        //      * mpsL is A(i-1) and belongs in envL later in the optimization step
        //      * mpsR is A(i) (or AC(i)) and is the active site
        auto &mpoL = model.get_mpo(posL);
        // auto                  &mpoR   = model.get_mpo(posR);
        Eigen::Tensor<cplx, 3> ML_PL  = mpsL.get_M() * mpsL.get_M().constant(cplx(res.alpha_mps, 0.0));
        Eigen::Tensor<cplx, 3> MR_PR  = mpsR.get_M();
        long                   PLdim2 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene != 0) {
            auto &envL = edges.get_env_eneL(posL);
            // auto &envR                       = edges.get_env_eneR(posR);
            long chi_max              = 2 * bond_lim / (1 + mpoL.MPO().dimension(1));
            chi_max                   = -1;
            Eigen::Tensor<cplx, 3> PL = envL.get_expansion_term(mpsL, mpoL, res.alpha_ene, chi_max);
            // Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, 1.0, chi_max);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            // Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL = std::move(ML_PL_tmp);
            // MR_PR                            = std::move(MR_PR_tmp);
            PLdim2 += PL.dimension(2);
            res.env_ene = envL;
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var != 0) {
            auto &envL = edges.get_env_varL(posL);
            // auto &envR                       = edges.get_env_varR(posR);
            long chi_max              = 2 * bond_lim / (1 + mpoL.MPO2().dimension(1));
            chi_max                   = -1;
            Eigen::Tensor<cplx, 3> PL = envL.get_expansion_term(mpsL, mpoL, res.alpha_var, chi_max);
            // Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, 1.0, chi_max);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            // Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL = std::move(ML_PL_tmp);
            // MR_PR                            = std::move(MR_PR_tmp);
            PLdim2 += PL.dimension(2);
            res.env_var = envL;
        }
        if(PLdim2 == 0) { tools::log->error("expand_environment_forward_1site: PRdim1 == 0: this is likely an error: mpsL was needlessly modified! "); }

        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsR.spin_dim(), PLdim2, mpsR.get_chiR());
        Eigen::Tensor<cplx, 3> MR_P0 = mpsR.get_M_bare().concatenate(P0, 1);
        merge_expansion_term_PL(state, mpsL, ML_PL, mpsR, MR_P0, svd_cfg);
        // Eigen::Tensor<cplx, 4> MP = ML_PL.contract(MR_PR, tenx::idx({2}, {1}));
        // merge_expanded_terms(state, mpsL, mpsR, MP, svd_cfg);

        state.clear_cache();
        state.clear_measurements();

        // if constexpr(settings::debug) {
        // state.clear_cache();
        // state.clear_measurements();
        // auto newS = tools::finite::measure::entanglement_entropy(mpsR.get_L());
        // if(std::abs(newS - oldS) > std::max(1e-4, alpha)) {
        //     tools::log->warn("Entanglement entropy changed by too much: {:.16f} -> {:.16f}, diff = {:.3e}", oldS, newS, newS - oldS);
        // }
        // }

        tools::log->debug("Environment expansion forward pos {} | {} | αₑ:{:.2e} αᵥ:{:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), dimR_old[1],
                          ML_PL.dimension(2), mpsR.get_chiL());
    }

    if(dimL_old[1] != mpsL.get_chiL()) throw except::runtime_error("mpsL changed chiL during environment expansion: {} -> {}", dimL_old, mpsL.dimensions());
    if(dimR_old[2] != mpsR.get_chiR()) throw except::runtime_error("mpsR changed chiR during environment expansion: {} -> {}", dimR_old, mpsR.dimensions());
    if constexpr(settings::debug_expansion) mpsL.assert_normalized();
    if constexpr(settings::debug_expansion) mpsR.assert_normalized();
    state.clear_cache();
    state.clear_measurements();
    env::rebuild_edges(state, model, edges);
    res.dimL_new = mpsL.dimensions();
    res.dimR_new = mpsR.dimensions();
    res.ene_new  = tools::finite::measure::energy(state, model, edges);
    res.var_new  = tools::finite::measure::energy_variance(state, model, edges);
    res.ok       = true;
    return res;
}

EnvExpansionResult tools::finite::env::expand_environment_reversed(StateFinite &state, const ModelFinite &model, EdgesFinite &edges,
                                                                   EnvExpandMode envExpandMode, svd::config svd_cfg) {
    if(state.get_position<long>() < 0) { return {}; }

    auto                pos = state.get_position<size_t>();
    std::vector<size_t> pos_expanded;
    if(state.get_direction() > 0 and pos == std::clamp<size_t>(pos, 1, state.get_length<size_t>() - 1)) pos_expanded = {pos - 1, pos};
    if(state.get_direction() < 0 and pos == std::clamp<size_t>(pos, 0, state.get_length<size_t>() - 2)) pos_expanded = {pos, pos + 1};
    if(pos_expanded.empty()) {
        return {}; // No expansion
    }

    // Define the left and right mps that will get modified
    state.clear_cache();
    state.clear_measurements();

    auto  posL     = pos_expanded.front();
    auto  posR     = pos_expanded.back();
    auto &mpsL     = state.get_mps_site(posL);
    auto &mpsR     = state.get_mps_site(posR);
    auto  dimL_old = mpsL.dimensions();
    auto  dimR_old = mpsR.dimensions();

    // Save expansion metadata
    auto res     = EnvExpansionResult();
    res.posL     = posL;
    res.posR     = posR;
    res.dimL_old = dimL_old;
    res.dimR_old = dimR_old;

    // Stop here if it's just a query
    if(envExpandMode == EnvExpandMode::NONE) return res; // Position query
    res.ene_old = tools::finite::measure::energy(state, model, edges);
    res.var_old = tools::finite::measure::energy_variance(state, model, edges);
    // res.alpha_ene = get_optimal_mixing_factor_ene(pos_expanded, state, model, edges, envExpandMode);
    // res.alpha_var = get_optimal_mixing_factor_var(pos_expanded, state, model, edges, envExpandMode);
    std::tie(res.alpha_mps, res.alpha_ene, res.alpha_var) = get_optimal_mixing_factors(pos_expanded, state, model, edges, envExpandMode);

    if(!(res.alpha_ene > 0) and !(res.alpha_var > 0)) {
        tools::log->debug("Expansion canceled: {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(),
                          mpsR.get_position(), res.alpha_ene, res.alpha_var);
        return res;
    }
    tools::log->debug("Expanding {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(), mpsR.get_position(),
                      res.alpha_ene, res.alpha_var);

    // Set up the SVD
    auto bond_lim = std::min(mpsL.spin_dim() * mpsL.get_chiL(), mpsR.spin_dim() * mpsR.get_chiR()); // Bond dimension can't grow faster than x spin_dim
    svd_cfg.truncation_limit = svd_cfg.truncation_limit.value_or(settings::precision::svd_truncation_min);
    svd_cfg.rank_max         = std::min(bond_lim, svd_cfg.rank_max.value_or(bond_lim));
    svd::solver svd(svd_cfg);

    if(state.get_direction() > 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is left-to-right:
        //      * mpsL is "A(i-1)" Gets expanded
        //      * mpsR is "AC(i)" Gets optimized
        auto &mpoL = model.get_mpo(posL);
        // auto                  &mpoR   = model.get_mpo(posR);
        Eigen::Tensor<cplx, 3> ML_PL = mpsL.get_M() * mpsL.get_M().constant(cplx(res.alpha_mps, 0.0)); // mpsL is going to be optimized, enriched with PL
        // Eigen::Tensor<cplx, 3> MR_PR  = mpsR.get_M(); // mpsL is going to be optimized, enriched with PL
        long PLdim2 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene > 0) {
            auto                  &envL    = edges.get_env_eneL(posL);
            auto                  &envR    = edges.get_env_eneR(posR);
            long                   chi_max = 2 * bond_lim / (1 + mpoL.MPO().dimension(1));
            Eigen::Tensor<cplx, 3> PL      = envL.get_expansion_term(mpsL, mpoL, res.alpha_ene, chi_max);
            // Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR,  res.alpha_ene, chi_max);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            // Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL = std::move(ML_PL_tmp);
            // MR_PR                            = std::move(MR_PR_tmp);
            PLdim2 += PL.dimension(2);
            res.env_ene = envL; // Saves a reference to the expanded environment
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var > 0) {
            auto                  &envL    = edges.get_env_varL(posL);
            auto                  &envR    = edges.get_env_varR(posR);
            long                   chi_max = 2 * bond_lim / (1 + mpoL.MPO2().dimension(1));
            Eigen::Tensor<cplx, 3> PL      = envL.get_expansion_term(mpsL, mpoL, res.alpha_var, chi_max);
            // Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, res.alpha_var, chi_max);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            // Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL = std::move(ML_PL_tmp);
            // MR_PR                            = std::move(MR_PR_tmp);
            PLdim2 += PL.dimension(2);
            res.env_var = envL; // Saves a reference to the expanded environment
        }
        if(PLdim2 == 0) { tools::log->error("expand_environment_reversed: PLdim2 == 0: this is likely an error"); }
        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsR.spin_dim(), PLdim2, mpsR.get_chiR());
        Eigen::Tensor<cplx, 3> MR_PR = mpsR.get_M_bare().concatenate(P0, 1);

        merge_expansion_term_PL(state, mpsL, ML_PL, mpsR, MR_PR, svd_cfg);
        tools::log->debug("Environment expansion reversed pos {} | {} | αₑ:{:.2e} αᵥ:{:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), dimL_old[2],
                          ML_PL.dimension(2), mpsL.get_chiR());
        if(dimL_old[1] != mpsL.get_chiL()) throw except::runtime_error("mpsL changed chiL during environment expansion: {} -> {}", dimL_old, mpsL.dimensions());
        if(dimR_old[2] != mpsR.get_chiR()) throw except::runtime_error("mpsR changed chiR during environment expansion: {} -> {}", dimR_old, mpsR.dimensions());
    }
    if(state.get_direction() < 0) {
        // The expanded bond sits between mpsL and mpsR. In 1-site DMRG, when the direction is right-to-left:
        //      * mpsL is "AC(i)" Gets optimized
        //      * mpsR is "B(i+1)" Gets expanded
        // auto                  &mpoL   = model.get_mpo(posL);
        auto &mpoR = model.get_mpo(posR);
        // Eigen::Tensor<cplx, 3> ML_PL  = mpsL.get_M(); // [ AC  P ]^T  including LC
        Eigen::Tensor<cplx, 3> MR_PR  = mpsR.get_M(); // [ AC  P ]^T  including LC
        long                   PRdim1 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene > 0) {
            // auto                  &envL       = edges.get_env_eneL(posL);
            auto &envR    = edges.get_env_eneR(posR);
            long  chi_max = 2 * bond_lim / (1 + mpoR.MPO().dimension(0));
            // Eigen::Tensor<cplx, 3> PL         = envL.get_expansion_term(mpsL, mpoL, res.alpha_ene, chi_max);
            Eigen::Tensor<cplx, 3> PR = envR.get_expansion_term(mpsR, mpoR, res.alpha_ene, chi_max);
            // Eigen::Tensor<cplx, 3> ML_PL_temp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_temp = MR_PR.concatenate(PR, 1);
            // ML_PL                             = std::move(ML_PL_temp);
            MR_PR = std::move(MR_PR_temp);
            PRdim1 += PR.dimension(1);
            res.env_ene = envR;
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var > 0) {
            // auto                  &envL       = edges.get_env_varL(posL);
            auto &envR    = edges.get_env_varR(posR);
            long  chi_max = 2 * bond_lim / (1 + mpoR.MPO2().dimension(0));
            // Eigen::Tensor<cplx, 3> PL         = envL.get_expansion_term(mpsL, mpoL, res.alpha_var, chi_max);
            Eigen::Tensor<cplx, 3> PR = envR.get_expansion_term(mpsR, mpoR, res.alpha_var, chi_max);
            // Eigen::Tensor<cplx, 3> ML_PL_temp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_temp = MR_PR.concatenate(PR, 1);
            // ML_PL                             = std::move(ML_PL_temp);
            MR_PR = std::move(MR_PR_temp);
            PRdim1 += PR.dimension(1);
            res.env_var = envR;
        }
        if(PRdim1 == 0) { tools::log->error("expand_environment_reversed: PRdim1 == 0: this is likely an error"); }
        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsL.spin_dim(), mpsL.get_chiL(), PRdim1);
        Eigen::Tensor<cplx, 3> ML_PL = mpsL.get_M().concatenate(P0, 2); // mpsL is going to be optimized, zero-padded with P0

        merge_expansion_term_PR(state, mpsL, ML_PL, mpsR, MR_PR, svd_cfg);

        tools::log->debug("Environment expansion reversed pos {} | {} | αₑ:{:.2e} αᵥ:{:.3e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), dimL_old[2],
                          MR_PR.dimension(1), mpsL.get_chiR());
    }
    if constexpr(settings::debug_expansion) mpsL.assert_normalized();
    if constexpr(settings::debug_expansion) mpsR.assert_normalized();
    env::rebuild_edges(state, model, edges);
    res.dimL_new = mpsL.dimensions();
    res.dimR_new = mpsR.dimensions();
    res.ene_new  = tools::finite::measure::energy(state, model, edges);
    res.var_new  = tools::finite::measure::energy_variance(state, model, edges);
    res.ok       = true;
    return res;
}

EnvExpansionResult tools::finite::env::expand_environment_forward_active(StateFinite &state, const ModelFinite &model, EdgesFinite &edges,
                                                                         EnvExpandMode envExpandMode, svd::config svd_cfg) {
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("expand_environment_forward_1site: All lengths not equal: state {} | model {} | edges {}", state.get_length(),
                                    model.get_length(), edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("expand_environment_forward_1site: All active sites are not equal: state {} | model {} | edges {}", state.active_sites,
                                    model.active_sites, edges.active_sites);
    if(state.active_sites.empty()) throw except::logic_error("No active sites for environment expansion");

    // Determine which bond to expand
    std::vector<size_t> pos_expanded;
    std::vector<size_t> pos_active_and_expanded;
    if(state.get_direction() > 0 and state.active_sites.back() + 1 < state.get_length<size_t>()) {
        auto pos     = state.active_sites.back();
        pos_expanded = {pos + 1};
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), state.active_sites.begin(), state.active_sites.end());
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), pos_expanded.begin(), pos_expanded.end());
    }
    if(state.get_direction() < 0 and state.active_sites.front() >= 1) {
        auto pos     = state.active_sites.front();
        pos_expanded = {pos - 1};
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), pos_expanded.begin(), pos_expanded.end());
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), state.active_sites.begin(), state.active_sites.end());
    }
    if(pos_expanded.empty()) {
        return {}; // No update
    }

    // Define the left and right mps that will get modified
    state.clear_cache();
    state.clear_measurements();
    // auto posL = pos_active_and_expanded.front();
    // auto posR = pos_active_and_expanded.back();
    // auto &mpsL     = state.get_mps_site(posL);
    // auto &mpsR     = state.get_mps_site(posR);
    // auto  dimL_old = mpsL.dimensions();
    // auto  dimR_old = mpsR.dimensions();

    // Save expansion metadata
    auto res = EnvExpansionResult();
    // res.posL     = posL;
    // res.posR     = posR;
    // res.dimL_old = dimL_old;
    // res.dimR_old = dimR_old;

    // Stop here if it's just a query
    if(envExpandMode == EnvExpandMode::NONE) return res; // Position query
    res.ene_old = tools::finite::measure::energy(state, model, edges);
    res.var_old = tools::finite::measure::energy_variance(state, model, edges);

    std::tie(res.alpha_mps, res.alpha_ene, res.alpha_var) = get_optimal_mixing_factors(pos_active_and_expanded, state, model, edges, envExpandMode);

    if(!(res.alpha_ene > 0) and !(res.alpha_var > 0)) {
        tools::log->debug("Expansion canceled: {} | αₑ:{:.2e} αᵥ:{:.2e}", pos_active_and_expanded, res.alpha_ene, res.alpha_var);
        return res;
    }
    tools::log->debug("Expanding {} | αₑ:{:.2e} αᵥ:{:.2e}", pos_active_and_expanded, res.alpha_ene, res.alpha_var);

    // Set up the SVD
    // Bond dimension can't grow faster than x spin_dim, but we can generate a highly enriched environment here for optimization,
    // and let the proper truncation happen after optimization instead.
    // auto bond_lim            = std::min(mpsL.spin_dim() * mpsL.get_chiL(), mpsR.spin_dim() * mpsR.get_chiR());
    svd_cfg.truncation_limit = svd_cfg.truncation_limit.value_or(settings::precision::svd_truncation_min);
    svd_cfg.rank_max         = svd_cfg.rank_max.value_or(settings::get_bond_max(state.get_algorithm()));
    svd::solver svd;

    if(state.get_direction() > 0) {
        // We expand the bond between active_sites.back() and the site to its right.
        auto  mpsL     = state.get_multisite_mps();
        auto &mpsR     = state.get_mps_site(pos_active_and_expanded.back());
        auto &mpoR     = model.get_mpo(pos_active_and_expanded.back());
        auto  dimL_old = mpsL.dimensions();
        auto  dimR_old = mpsR.dimensions();

        Eigen::Tensor<cplx, 3> ML_PL  = mpsL;         // mpsR is going into the environment, enriched with PR.
        Eigen::Tensor<cplx, 3> MR_PR  = mpsR.get_M(); // mpsR is going into the environment, enriched with PR.
        long                   PRdim1 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene > 0) {
            const auto            &mpoL      = model.get_multisite_mpo();
            const auto            &envL      = edges.get_ene_active().L;
            const auto            &envR      = edges.get_env_eneR(pos_active_and_expanded.back());
            Eigen::Tensor<cplx, 3> PL        = envL.get_expansion_term(mpsL, mpoL, 1.0);
            Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, res.alpha_ene, -1);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL                            = std::move(ML_PL_tmp);
            MR_PR                            = std::move(MR_PR_tmp);
            PRdim1 += PR.dimension(1);
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var > 0) {
            const auto            &mpoL      = model.get_multisite_mpo();
            const auto            &envL      = edges.get_var_active().L;
            const auto            &envR      = edges.get_env_varR(pos_active_and_expanded.back());
            Eigen::Tensor<cplx, 3> PL        = envL.get_expansion_term(mpsL, mpoL, 1.0);
            Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, res.alpha_var, -1);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL                            = std::move(ML_PL_tmp);
            MR_PR                            = std::move(MR_PR_tmp);
            PRdim1 += PR.dimension(1);
        }
        // merge_expansion_term_PR(state, mpsL, ML_PL, mpsR, MR_PR, svd_cfg);
        auto [U, S, V] = svd.schmidt_into_right_normalized(MR_PR, mpsR.spin_dim(), svd_cfg);

        if(mpsR.get_label() == "B") {
            Eigen::Tensor<cplx, 2> US =
                U.contract(tenx::asDiagonal(S), tenx::idx({2}, {0})).reshape(tenx::array2{U.dimension(0) * U.dimension(1), S.dimension(0)});
            Eigen::Tensor<cplx, 3> mpsL_US = mpsL.contract(US, tenx::idx({2}, {0}));
            mpsR.set_M(V);
            mps::merge_multisite_mps(state, mpsL_US, state.active_sites, state.get_position<long>(), MergeEvent::NORM, svd_cfg);
        } else {
            throw except::runtime_error("Expansion not implemented for label mpsR={}", mpsR.get_label());
        }

        // Eigen::Tensor<cplx, 4> MP = ML_PL.contract(MR_PR, tenx::idx({2}, {1}));
        // merge_expanded_terms(state, mpsL, mpsR, MP, svd_cfg);

        state.clear_measurements();
        state.clear_cache();

        tools::log->debug("Environment expansion forward pos {} | {} | αₑ:{:.2e} αᵥ:{:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}",
                          pos_active_and_expanded, flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(),
                          svd_cfg.rank_max.value(), dimL_old[2], PRdim1, mpsR.get_chiL());
    }
    if(state.get_direction() < 0) {
        auto  posL     = pos_active_and_expanded.front();
        auto  posR     = pos_active_and_expanded.back();
        auto &mpsL     = state.get_mps_site(posL);
        auto &mpsR     = state.get_mps_site(posR);
        auto  dimL_old = mpsL.dimensions();
        auto  dimR_old = mpsR.dimensions();
        // The expanded bond sits between mpsL and mpsR. When direction is right-to-left:
        //      * mpsL is A(i-1) and belongs in envL later in the optimization step
        //      * mpsR is A(i) (or AC(i)) and is the active site
        auto                  &mpoL   = model.get_mpo(posL);
        auto                  &mpoR   = model.get_mpo(posR);
        Eigen::Tensor<cplx, 3> ML_PL  = mpsL.get_M();
        Eigen::Tensor<cplx, 3> MR_PR  = mpsR.get_M();
        long                   PLdim2 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene > 0) {
            auto                  &envL      = edges.get_env_eneL(posL);
            auto                  &envR      = edges.get_env_eneR(posR);
            long                   chi_max   = 2 * svd_cfg.rank_max.value() / (1 + mpoL.MPO().dimension(1));
            Eigen::Tensor<cplx, 3> PL        = envL.get_expansion_term(mpsL, mpoL, res.alpha_ene, chi_max);
            Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, 1.0, chi_max);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL                            = std::move(ML_PL_tmp);
            MR_PR                            = std::move(MR_PR_tmp);
            PLdim2 += PL.dimension(2);
            res.env_ene = envL;
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var > 0) {
            auto                  &envL      = edges.get_env_varL(posL);
            auto                  &envR      = edges.get_env_varR(posR);
            long                   chi_max   = 2 * svd_cfg.rank_max.value() / (1 + mpoL.MPO2().dimension(1));
            Eigen::Tensor<cplx, 3> PL        = envL.get_expansion_term(mpsL, mpoL, res.alpha_var, chi_max);
            Eigen::Tensor<cplx, 3> PR        = envR.get_expansion_term(mpsR, mpoR, 1.0, chi_max);
            Eigen::Tensor<cplx, 3> ML_PL_tmp = ML_PL.concatenate(PL, 2);
            Eigen::Tensor<cplx, 3> MR_PR_tmp = MR_PR.concatenate(PR, 1);
            ML_PL                            = std::move(ML_PL_tmp);
            MR_PR                            = std::move(MR_PR_tmp);
            PLdim2 += PL.dimension(2);
            res.env_var = envL;
        }
        if(PLdim2 == 0) { tools::log->error("expand_environment_forward_1site: PRdim1 == 0: this is likely an error: mpsL was needlessly modified! "); }

        // Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsR.spin_dim(), PLdim2, mpsR.get_chiR());
        // Eigen::Tensor<cplx, 3> MR_P0 = mpsR.get_M_bare().concatenate(P0, 1);
        merge_expansion_term_PL(state, mpsL, ML_PL, mpsR, MR_PR, svd_cfg);

        state.clear_cache();
        state.clear_measurements();

        // if constexpr(settings::debug) {
        // state.clear_cache();
        // state.clear_measurements();
        // auto newS = tools::finite::measure::entanglement_entropy(mpsR.get_L());
        // if(std::abs(newS - oldS) > std::max(1e-4, alpha)) {
        //     tools::log->warn("Entanglement entropy changed by too much: {:.16f} -> {:.16f}, diff = {:.3e}", oldS, newS, newS - oldS);
        // }
        // }

        tools::log->debug("Environment expansion forward pos {} | {} | αₑ:{:.2e} αᵥ:{:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), dimR_old[1],
                          PLdim2, mpsR.get_chiL());
    }

    // if(dimL_old[1] != mpsL.get_chiL()) throw except::runtime_error("mpsL changed chiL during environment expansion: {} -> {}", dimL_old, mpsL.dimensions());
    // if(dimR_old[2] != mpsR.get_chiR()) throw except::runtime_error("mpsR changed chiR during environment expansion: {} -> {}", dimR_old, mpsR.dimensions());
    // if constexpr(settings::debug_expansion) mpsL.assert_normalized();
    // if constexpr(settings::debug_expansion) mpsR.assert_normalized();
    state.clear_cache();
    state.clear_measurements();
    env::rebuild_edges(state, model, edges);
    // res.dimL_new = mpsL.dimensions();
    // res.dimR_new = mpsR.dimensions();
    res.ene_new = tools::finite::measure::energy(state, model, edges);
    res.var_new = tools::finite::measure::energy_variance(state, model, edges);
    res.ok      = true;
    return res;
}

EnvExpansionResult tools::finite::env::expand_environment_forward_old(StateFinite &state, const ModelFinite &model, EdgesFinite &edges,
                                                                      EnvExpandMode envExpandMode, svd::config svd_cfg) {
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("expand_environment_forward_old: All lengths not equal: state {} | model {} | edges {}", state.get_length(),
                                    model.get_length(), edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("expand_environment_forward_old: All active sites are not equal: state {} | model {} | edges {}", state.active_sites,
                                    model.active_sites, edges.active_sites);
    if(state.active_sites.empty()) throw except::logic_error("No active sites for environment expansion");

    // Determine which bond to expand
    std::vector<size_t> pos_expanded;
    std::vector<size_t> pos_active_and_expanded;
    if(state.get_direction() > 0) {
        auto pos = state.active_sites.back();
        if(pos == std::clamp<size_t>(pos, 0, state.get_length<size_t>() - 2)) pos_expanded = {pos, pos + 1};
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), state.active_sites.begin(), state.active_sites.end());
        pos_active_and_expanded.emplace_back(pos + 1);
    }
    if(state.get_direction() < 0) {
        auto pos = state.active_sites.front();
        if(pos == std::clamp<size_t>(pos, 1, state.get_length<size_t>() - 1)) pos_expanded = {pos - 1, pos};
        pos_active_and_expanded.emplace_back(pos - 1);
        pos_active_and_expanded.insert(pos_active_and_expanded.end(), state.active_sites.begin(), state.active_sites.end());
    }
    if(pos_expanded.empty()) {
        return {}; // No update
    }

    // Define the left and right mps that will get modified
    state.clear_cache();
    state.clear_measurements();
    auto  posL     = pos_expanded.front();
    auto  posR     = pos_expanded.back();
    auto &mpsL     = state.get_mps_site(posL);
    auto &mpsR     = state.get_mps_site(posR);
    auto  dimL_old = mpsL.dimensions();
    auto  dimR_old = mpsR.dimensions();

    // Save expansion metadata
    auto res     = EnvExpansionResult();
    res.posL     = posL;
    res.posR     = posR;
    res.dimL_old = dimL_old;
    res.dimR_old = dimR_old;

    // Stop here if it's just a query
    if(envExpandMode == EnvExpandMode::NONE) return res; // Position query
    res.ene_old = tools::finite::measure::energy(state, model, edges);
    res.var_old = tools::finite::measure::energy_variance(state, model, edges);

    std::tie(res.alpha_mps, res.alpha_ene, res.alpha_mps) = get_optimal_mixing_factors(pos_active_and_expanded, state, model, edges, envExpandMode);
    // res.alpha_ene = get_optimal_mixing_factor_ene(pos_active_and_expanded, state, model, edges, envExpandMode);
    // res.alpha_var = get_optimal_mixing_factor_var(pos_active_and_expanded, state, model, edges, envExpandMode);
    if(!(res.alpha_ene > 0) and !(res.alpha_var > 0)) {
        tools::log->debug("Expansion canceled: {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(),
                          mpsR.get_position(), res.alpha_ene, res.alpha_var);
        return res;
    }
    tools::log->debug("Expanding {}({}) - {}({}) | αₑ:{:.2e} αᵥ:{:.2e}", mpsL.get_label(), mpsL.get_position(), mpsR.get_label(), mpsR.get_position(),
                      res.alpha_ene, res.alpha_var);

    // Set up the SVD
    // Bond dimension can't grow faster than x spin_dim, but we can generate a highly enriched environment here for optimization,
    // and let the proper truncation happen after optimization instead.
    auto bond_lim            = std::min(mpsL.spin_dim() * mpsL.get_chiL(), mpsR.spin_dim() * mpsR.get_chiR());
    svd_cfg.truncation_limit = svd_cfg.truncation_limit.value_or(settings::precision::svd_truncation_min);
    svd_cfg.rank_max         = std::min(bond_lim, svd_cfg.rank_max.value_or(bond_lim));
    svd::solver svd;

    if(state.get_direction() > 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is left-to-right:
        // On single-site DMRG
        //      * mpsL is AC(i), or B(i) during multisite dmrg
        //      * mpsR is B(i+1) and belongs in envR later in the optimization step
        Eigen::Tensor<cplx, 3> MR_PR  = mpsR.get_M(); // mpsR is going into the environment, enriched with PR.
        auto                  &mpoR   = model.get_mpo(posR);
        long                   PRdim1 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene > 0) {
            auto                  &env         = edges.get_env_eneR(posR);
            long                   PR_chiL_max = 2 * bond_lim / (1 + mpoR.MPO().dimension(0));
            Eigen::Tensor<cplx, 3> PR          = env.get_expansion_term(mpsR, mpoR, res.alpha_ene, PR_chiL_max);
            Eigen::Tensor<cplx, 3> MR_PR_tmp   = MR_PR.concatenate(PR, 1);
            MR_PR                              = std::move(MR_PR_tmp);
            PRdim1 += PR.dimension(1);
            res.env_ene = env;
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var > 0) {
            auto                  &env         = edges.get_env_varR(posR);
            long                   PR_chiL_max = 2 * bond_lim / (1 + mpoR.MPO2().dimension(0));
            Eigen::Tensor<cplx, 3> PR          = env.get_expansion_term(mpsR, mpoR, res.alpha_var, PR_chiL_max);
            Eigen::Tensor<cplx, 3> MR_PR_tmp   = MR_PR.concatenate(PR, 1);
            MR_PR                              = std::move(MR_PR_tmp);
            PRdim1 += PR.dimension(1);
            res.env_var = env;
        }
        if(PRdim1 == 0) { tools::log->error("expand_environment_forward_old: PRdim1 == 0: this is likely an error: mpsR was needlessly modified! "); }
        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsL.spin_dim(), mpsL.get_chiL(), PRdim1);
        Eigen::Tensor<cplx, 3> ML_P0 = mpsL.get_M().concatenate(P0, 2); // mpsL is going to be optimized, zero-padded with P0

        merge_expansion_term_PR(state, mpsL, ML_P0, mpsR, MR_PR, svd_cfg);

        state.clear_measurements();
        state.clear_cache();

        tools::log->debug("Environment expansion forward pos {} | {} | αₑ:{:.2e} αᵥ:{:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), dimL_old[2],
                          ML_P0.dimension(2), mpsL.get_chiR());
    }
    if(state.get_direction() < 0) {
        // The expanded bond sits between mpsL and mpsR. When direction is right-to-left:
        //      * mpsL is A(i-1) and belongs in envL later in the optimization step
        //      * mpsR is A(i) (or AC(i)) and is the active site
        auto                  &mpoL   = model.get_mpo(posL);
        Eigen::Tensor<cplx, 3> ML_PL  = mpsL.get_M();
        long                   PLdim2 = 0;
        if(has_flag(envExpandMode, EnvExpandMode::ENE) and res.alpha_ene > 0) {
            auto                  &env         = edges.get_env_eneL(posL);
            long                   PL_chiR_max = 2 * bond_lim / (1 + mpoL.MPO().dimension(1));
            Eigen::Tensor<cplx, 3> PL          = env.get_expansion_term(mpsL, mpoL, res.alpha_ene, PL_chiR_max);
            Eigen::Tensor<cplx, 3> ML_PL_tmp   = ML_PL.concatenate(PL, 2);
            ML_PL                              = std::move(ML_PL_tmp);
            PLdim2 += PL.dimension(2);
            res.env_ene = env;
        }
        if(has_flag(envExpandMode, EnvExpandMode::VAR) and res.alpha_var > 0) {
            auto                  &env         = edges.get_env_varL(posL);
            long                   PL_chiR_max = 2 * bond_lim / (1 + mpoL.MPO2().dimension(1));
            Eigen::Tensor<cplx, 3> PL          = env.get_expansion_term(mpsL, mpoL, res.alpha_var, PL_chiR_max);
            Eigen::Tensor<cplx, 3> ML_PL_tmp   = ML_PL.concatenate(PL, 2);
            ML_PL                              = std::move(ML_PL_tmp);
            PLdim2 += PL.dimension(2);
            res.env_var = env;
        }
        if(PLdim2 == 0) { tools::log->error("expand_environment_forward_old: PRdim1 == 0: this is likely an error: mpsL was needlessly modified! "); }

        Eigen::Tensor<cplx, 3> P0    = tenx::TensorConstant<cplx>(cplx(0.0, 0.0), mpsR.spin_dim(), PLdim2, mpsR.get_chiR());
        Eigen::Tensor<cplx, 3> MR_P0 = mpsR.get_M_bare().concatenate(P0, 1);
        merge_expansion_term_PL(state, mpsL, ML_PL, mpsR, MR_P0, svd_cfg);

        state.clear_cache();
        state.clear_measurements();

        // if constexpr(settings::debug) {
        // state.clear_cache();
        // state.clear_measurements();
        // auto newS = tools::finite::measure::entanglement_entropy(mpsR.get_L());
        // if(std::abs(newS - oldS) > std::max(1e-4, alpha)) {
        //     tools::log->warn("Entanglement entropy changed by too much: {:.16f} -> {:.16f}, diff = {:.3e}", oldS, newS, newS - oldS);
        // }
        // }

        tools::log->debug("Environment expansion forward pos {} | {} | αₑ:{:.2e} αᵥ:{:.2e} | svd_ε {:.2e} | χlim {} | χ {} -> {} -> {}", pos_expanded,
                          flag2str(envExpandMode), res.alpha_ene, res.alpha_var, svd_cfg.truncation_limit.value(), svd_cfg.rank_max.value(), dimR_old[1],
                          MR_P0.dimension(1), mpsR.get_chiL());
    }

    if(dimL_old[1] != mpsL.get_chiL()) throw except::runtime_error("mpsL changed chiL during environment expansion: {} -> {}", dimL_old, mpsL.dimensions());
    if(dimR_old[2] != mpsR.get_chiR()) throw except::runtime_error("mpsR changed chiR during environment expansion: {} -> {}", dimR_old, mpsR.dimensions());
    if constexpr(settings::debug_expansion) mpsL.assert_normalized();
    if constexpr(settings::debug_expansion) mpsR.assert_normalized();
    state.clear_cache();
    state.clear_measurements();
    env::rebuild_edges(state, model, edges);
    res.dimL_new = mpsL.dimensions();
    res.dimR_new = mpsR.dimensions();
    res.ene_new  = tools::finite::measure::energy(state, model, edges);
    res.var_new  = tools::finite::measure::energy_variance(state, model, edges);
    res.ok       = true;
    return res;
}

void tools::finite::env::assert_edges_ene(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) throw except::logic_error("assert_edges_var: fLBIT algorithm should never assert energy edges!");
    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites we shouldn't be asserting edges.
    // For instance, the active sites are cleared after a move of center site.
    // We could always keep all edges refreshed but that would be wasteful, since the next iteration
    // may activate other sites and not end up needing those edges.
    // Instead, we force the hand of the algorithm, to only allow edge assertions with active sites defined.
    // Ideally, then, this should be done directly after activating new sites in a new iteration.
    if(edges.active_sites.empty())
        throw except::runtime_error("assert_edges_ene: no active sites.\n"
                                    "Hint:\n"
                                    " One could in principle keep edges refreshed always, but\n"
                                    " that would imply rebuilding many edges that end up not\n"
                                    " being used. Make sure to only run this assertion after\n"
                                    " activating sites.");

    long current_position = state.get_position<long>();

    // size_t posL_active      = edges.active_sites.front();
    // size_t posR_active      = edges.active_sites.back();

    // These back and front positions will seem reversed: we need extra edges for optimal subspace expansion: see the Log from 2024-07-23
    size_t posL_active = edges.active_sites.back();
    size_t posR_active = edges.active_sites.front();
    if constexpr(settings::debug_edges)
        tools::log->trace("assert_edges_ene: pos {} | dir {} | "
                          "asserting edges eneL from [{} to {}]",
                          current_position, state.get_direction(), min_pos, posL_active);

    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &ene = edges.get_env_eneL(pos);
        if(pos == 0 and not ene.has_block()) throw except::runtime_error("ene L at pos {} does not have a block", pos);
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &ene_next = edges.get_env_eneL(pos + 1);
        ene_next.assert_unique_id(ene, mps, mpo);
    }
    if constexpr(settings::debug_edges)
        tools::log->trace("assert_edges_ene: pos {} | dir {} | "
                          "asserting edges eneR from [{} to {}]",
                          current_position, state.get_direction(), posR_active, max_pos);

    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); --pos) {
        auto &ene = edges.get_env_eneR(pos);
        if(pos == state.get_length() - 1 and not ene.has_block()) throw except::runtime_error("ene R at pos {} does not have a block", pos);
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &ene_prev = edges.get_env_eneR(pos - 1);
        ene_prev.assert_unique_id(ene, mps, mpo);
    }
}

void tools::finite::env::assert_edges_var(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) throw except::logic_error("assert_edges_var: fLBIT algorithm should never assert variance edges!");
    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites we shouldn't be asserting edges.
    // For instance, the active sites are cleared after a move of center site.
    // We could always keep all edges refreshed but that would be wasteful, since the next iteration
    // may activate other sites and not end up needing those edges.
    // Instead, we force the hand of the algorithm, to only allow edge assertions with active sites defined.
    // Ideally, then, this should be done directly after activating new sites in a new iteration.
    if(edges.active_sites.empty())
        throw except::runtime_error("assert_edges_var: no active sites.\n"
                                    "Hint:\n"
                                    " One could in principle keep edges refreshed always, but\n"
                                    " that would imply rebuilding many edges that end up not\n"
                                    " being used. Make sure to only run this assertion after\n"
                                    " activating sites.");

    long current_position = state.get_position<long>();
    // size_t posL_active      = edges.active_sites.front();
    // size_t posR_active      = edges.active_sites.back();

    // These back and front positions will seem reversed: we need extra edges for optimal subspace expansion: see the Log from 2024-07-23
    size_t posL_active = edges.active_sites.back();
    size_t posR_active = edges.active_sites.front();

    if constexpr(settings::debug_edges)
        tools::log->trace("assert_edges_var: pos {} | dir {} | "
                          "asserting edges varL from [{} to {}]",
                          current_position, state.get_direction(), min_pos, posL_active);
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &var = edges.get_env_varL(pos);
        if(pos == 0 and not var.has_block()) throw except::runtime_error("var L at pos {} does not have a block", pos);
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &var_next = edges.get_env_varL(pos + 1);
        var_next.assert_unique_id(var, mps, mpo);
    }
    if constexpr(settings::debug_edges)
        tools::log->trace("assert_edges_var: pos {} | dir {} | "
                          "asserting edges varR from [{} to {}]",
                          current_position, state.get_direction(), posR_active, max_pos);
    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); --pos) {
        auto &var = edges.get_env_varR(pos);
        if(pos == state.get_length() - 1 and not var.has_block()) throw except::runtime_error("var R at pos {} does not have a block", pos);
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &mps      = state.get_mps_site(pos);
        auto &mpo      = model.get_mpo(pos);
        auto &var_prev = edges.get_env_varR(pos - 1);
        var_prev.assert_unique_id(var, mps, mpo);
    }
}

void tools::finite::env::assert_edges(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) return;
    assert_edges_ene(state, model, edges);
    assert_edges_var(state, model, edges);
}

void tools::finite::env::rebuild_edges_ene(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) throw except::logic_error("rebuild_edges_ene: fLBIT algorithm should never rebuild energy edges!");
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(), edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites,
                                    edges.active_sites);
    auto   t_reb   = tid::tic_scope("rebuild_edges_ene", tid::higher);
    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites then we can build up until current position
    /*
     * LOG:
     * - 2021-10-14:
     *      Just had a terribly annoying bug:
     *      Moving the center position clears active_sites, which caused problems when turning back from the right edge.
     *          1)  active_sites [A(L-1), AC(L)] are updated, left edge exist for A(L-1), right edge exists for AC(L)
     *          2)  move dir -1, clear active sites
     *          3)  assert_edges checks up to AC(L-1), but this site has a stale right edge.
     *      Therefore, one would have to rebuild edges between steps 2) and 3) to solve this issue
     *
     *      One solution would be to always rebuild edges up to the current position from both sides, but that would be
     *      wasteful. Instead, we could just accept that some edges are stale after moving the center-point,
     *      as long as we rebuild those when sites get activated again.
     *
     * - 2024-07-23
     *      Just found a way to calculate the optimal mixing factor for subspace expansion.
     *      In forward expansion we need H_eff including one site beyond active_sites.
     *      Therefore, we need to build more environments than we have needed previously.
     *      Examples:
     *          - Forward, direction == 1, active_sites = [5,6,7,8].
     *            Then we expand bond [8,9], and so we need envL[8] and envR[9].
     *          - Forward, direction == -1, active_sites = [5,6,7,8].
     *            Then we expand bond [4,5], and so we need envL[4] and envR[5].
     *      These environments weren't built before.
     *      Therefore, we must now rebuild
     *          - envL 0 to active_sites.back()
     *          - envR L to active_sites.front()
     */

    // If there are no active sites we shouldn't be rebuilding edges.
    // For instance, the active sites are cleared after a move of center site.
    // We could always keep all edges refreshed but that would be wasteful, since the next iteration
    // may activate other sites and not end up needing those edges.
    // Instead, we force the hand of the algorithm, to only allow edge rebuilds with active sites defined.
    // Ideally, then, this should be done directly after activating new sites in a new iteration.
    if(edges.active_sites.empty())
        throw except::runtime_error("rebuild_edges_ene: no active sites.\n"
                                    "Hint:\n"
                                    " One could in principle keep edges refreshed always, but\n"
                                    " that would imply rebuilding many edges that end up not\n"
                                    " being used. Make sure to only run this rebuild after\n"
                                    " activating sites.");

    long current_position = state.get_position<long>();
    // size_t posL_active      = edges.active_sites.front();
    // size_t posR_active      = edges.active_sites.back();

    // These back and front positions will seem reversed: we need extra edges for optimal subspace expansion: see the Log from 2024-07-23
    size_t posL_active = edges.active_sites.back();
    size_t posR_active = edges.active_sites.front();
    if constexpr(settings::debug_edges)
        tools::log->trace("rebuild_edges_ene: pos {} | dir {} | "
                          "inspecting edges eneL from [{} to {}]",
                          current_position, state.get_direction(), min_pos, posL_active);
    std::vector<size_t> env_pos_log;
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &env_here = edges.get_env_eneL(pos);
        auto  id_here  = env_here.get_unique_id();
        if(pos == 0) env_here.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        if(not env_here.has_block()) throw except::runtime_error("rebuild_edges_ene: No eneL block detected at pos {}", pos);
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &env_rght = edges.get_env_eneL(pos + 1);
        auto  id_rght  = env_rght.get_unique_id();
        env_rght.refresh(env_here, state.get_mps_site(pos), model.get_mpo(pos));
        if(id_here != env_here.get_unique_id()) env_pos_log.emplace_back(env_here.get_position());
        if(id_rght != env_rght.get_unique_id()) env_pos_log.emplace_back(env_rght.get_position());
    }
    if(not env_pos_log.empty()) tools::log->trace("rebuild_edges_ene: rebuilt eneL edges: {}", env_pos_log);

    env_pos_log.clear();
    if constexpr(settings::debug_edges)
        tools::log->trace("rebuild_edges_ene: pos {} | dir {} | "
                          "inspecting edges eneR from [{} to {}]",
                          current_position, state.get_direction(), posR_active, max_pos);
    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); --pos) {
        auto &env_here = edges.get_env_eneR(pos);
        auto  id_here  = env_here.get_unique_id();
        if(pos == state.get_length() - 1) env_here.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        if(not env_here.has_block()) throw except::runtime_error("rebuild_edges_ene: No eneR block detected at pos {}", pos);
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &env_left = edges.get_env_eneR(pos - 1);
        auto  id_left  = env_left.get_unique_id();
        env_left.refresh(env_here, state.get_mps_site(pos), model.get_mpo(pos));
        if(id_here != env_here.get_unique_id()) env_pos_log.emplace_back(env_here.get_position());
        if(id_left != env_left.get_unique_id()) env_pos_log.emplace_back(env_left.get_position());
    }
    std::reverse(env_pos_log.begin(), env_pos_log.end());
    if(not env_pos_log.empty()) tools::log->trace("rebuild_edges_ene: rebuilt eneR edges: {}", env_pos_log);
    if(not edges.get_env_eneL(posL_active).has_block()) throw except::logic_error("rebuild_edges_ene: active env eneL has undefined block");
    if(not edges.get_env_eneR(posR_active).has_block()) throw except::logic_error("rebuild_edges_ene: active env eneR has undefined block");
}

void tools::finite::env::rebuild_edges_var(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) throw except::logic_error("rebuild_edges_var: fLBIT algorithm should never rebuild variance edges!");
    if(not num::all_equal(state.get_length(), model.get_length(), edges.get_length()))
        throw except::runtime_error("rebuild_edges_var: All lengths not equal: state {} | model {} | edges {}", state.get_length(), model.get_length(),
                                    edges.get_length());
    if(not num::all_equal(state.active_sites, model.active_sites, edges.active_sites))
        throw except::runtime_error("rebuild_edges_var: All active sites are not equal: state {} | model {} | edges {}", state.active_sites, model.active_sites,
                                    edges.active_sites);
    auto t_reb = tid::tic_scope("rebuild_edges_var", tid::level::higher);

    size_t min_pos = 0;
    size_t max_pos = state.get_length() - 1;

    // If there are no active sites we shouldn't be rebuilding edges.
    // For instance, the active sites are cleared after a move of center site.
    // We could always keep all edges refreshed but that would be wasteful, since the next iteration
    // may activate other sites and not end up needing those edges.
    // Instead, we force the hand of the algorithm, to only allow edge rebuilds with active sites defined.
    // Ideally, then, this should be done directly after activating new sites in a new iteration.
    if(edges.active_sites.empty())
        throw except::runtime_error("rebuild_edges_var: no active sites.\n"
                                    "Hint:\n"
                                    " One could in principle keep edges refreshed always, but\n"
                                    " that would imply rebuilding many edges that end up not\n"
                                    " being used. Make sure to only run this assertion after\n"
                                    " activating sites.");

    long current_position = state.get_position<long>();
    // size_t posL_active      = edges.active_sites.front();
    // size_t posR_active      = edges.active_sites.back();

    // These back and front positions will seem reversed: we need extra edges for optimal subspace expansion: see the Log from 2024-07-23
    size_t posL_active = edges.active_sites.back();
    size_t posR_active = edges.active_sites.front();
    if constexpr(settings::debug_edges) {
        tools::log->trace("rebuild_edges_var: pos {} | dir {} | "
                          "inspecting edges varL from [{} to {}]",
                          current_position, state.get_direction(), min_pos, posL_active);
    }

    std::vector<size_t> env_pos_log;
    for(size_t pos = min_pos; pos <= posL_active; pos++) {
        auto &env_here = edges.get_env_varL(pos);
        auto  id_here  = env_here.get_unique_id();
        if(pos == 0) env_here.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        if(not env_here.has_block()) throw except::runtime_error("rebuild_edges_var: No varL block detected at pos {}", pos);
        if(pos >= std::min(posL_active, state.get_length() - 1)) continue;
        auto &env_rght = edges.get_env_varL(pos + 1);
        auto  id_rght  = env_rght.get_unique_id();
        env_rght.refresh(env_here, state.get_mps_site(pos), model.get_mpo(pos));
        if(id_here != env_here.get_unique_id()) env_pos_log.emplace_back(env_here.get_position());
        if(id_rght != env_rght.get_unique_id()) env_pos_log.emplace_back(env_rght.get_position());
    }

    if(not env_pos_log.empty()) tools::log->trace("rebuild_edges_var: rebuilt varL edges: {}", env_pos_log);
    env_pos_log.clear();
    if constexpr(settings::debug_edges) {
        tools::log->trace("rebuild_edges_var: pos {} | dir {} | "
                          "inspecting edges varR from [{} to {}]",
                          current_position, state.get_direction(), posR_active, max_pos);
    }

    for(size_t pos = max_pos; pos >= posR_active and pos < state.get_length(); --pos) {
        auto &env_here = edges.get_env_varR(pos);
        auto  id_here  = env_here.get_unique_id();
        if(pos == state.get_length() - 1) env_here.set_edge_dims(state.get_mps_site(pos), model.get_mpo(pos));
        if(not env_here.has_block()) throw except::runtime_error("rebuild_edges_var: No varR block detected at pos {}", pos);
        if(pos <= std::max(posR_active, 0ul)) continue;
        auto &env_left = edges.get_env_varR(pos - 1);
        auto  id_left  = env_left.get_unique_id();
        env_left.refresh(env_here, state.get_mps_site(pos), model.get_mpo(pos));
        if(id_here != env_here.get_unique_id()) env_pos_log.emplace_back(env_here.get_position());
        if(id_left != env_left.get_unique_id()) env_pos_log.emplace_back(env_left.get_position());
    }
    std::reverse(env_pos_log.begin(), env_pos_log.end());
    if(not env_pos_log.empty()) tools::log->trace("rebuild_edges_var: rebuilt varR edges: {}", env_pos_log);
    if(not edges.get_env_varL(posL_active).has_block()) throw except::logic_error("rebuild_edges_var: active env varL has undefined block");
    if(not edges.get_env_varR(posR_active).has_block()) throw except::logic_error("rebuild_edges_var: active env varR has undefined block");
}

void tools::finite::env::rebuild_edges(const StateFinite &state, const ModelFinite &model, EdgesFinite &edges) {
    if(state.get_algorithm() == AlgorithmType::fLBIT) return;
    rebuild_edges_ene(state, model, edges);
    rebuild_edges_var(state, model, edges);
}
