#include "../lbit.h"
#include "../spin.h"
#include "config/debug.h"
#include "config/enums.h"
// #include "config/settings.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "io/fmt.h"
#include "io/spdlog.h"
#include "math/fit.h"
#include "math/linalg.h"
#include "math/num.h"
#include "math/rnd.h"
#include "math/stat.h"
#include "math/svd.h"
#include "math/tenx.h"
#include "tensors/site/mps/MpsSite.h"
#include "tensors/state/StateFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/mps.h"
#include "tools/finite/ops.h"
#include <algorithm>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

namespace settings {
    inline constexpr bool debug_circuit = false;
}

using cplx = qm::cplx;

qm::lbit::UnitaryGateProperties::UnitaryGateProperties(const std::vector<double> &h)
    : sites(settings::model::model_size), depth(settings::model::lbit::u_depth), fmix(settings::model::lbit::u_fmix), tstd(settings::model::lbit::u_tstd),
      cstd(settings::model::lbit::u_cstd), tgw8(settings::model::lbit::u_tgw8), cgw8(settings::model::lbit::u_cgw8), hmean(settings::model::lbit::J1_mean),
      hwdth(settings::model::lbit::J1_wdth), hdist(settings::model::lbit::distribution), hvals(h) {
    if(tgw8 == UnitaryGateWeight::EXPDECAY and hvals.empty()) throw except::logic_error("No onsite fields h given for t UnitaryGateWeight::EXPDECAY");
    if(cgw8 == UnitaryGateWeight::EXPDECAY and hvals.empty()) throw except::logic_error("No onsite fields h given for c UnitaryGateWeight::EXPDECAY");
}

std::string qm::lbit::UnitaryGateProperties::string() const {
    return fmt::format("depth {} | fmix {:.3f} | tstd {:.3f} | cstd {:.3f} | tw8 {} | cw8 {} | hvals {::.3e}", depth, fmix, tstd, cstd, enum2sv(tgw8),
                       enum2sv(cgw8), hvals);
}

void qm::lbit::UnitaryGateProperties::randomize_hvals() const { hvals = rnd::random(hdist, hmean, hwdth, sites); }

// std::vector<qm::Gate> qm::lbit::get_unitary_2gate_layer(size_t sites, double fmix) {
//     /*! Returns a set of unitary two site operators used to transform between physical and l-bit representations
//      *
//      *
//      @verbatim
//                    0      1                0
//                    |      |                |
//                  [ exp(-ifM) ]  ---> [ exp(-ifM) ]
//                    |      |               |
//                    2      3               1
//
//      @endverbatim
//
//      Where M is defined as
//         θ³n[i]n[i+1] +
//         θ²n[i](1-n[i+1]) +
//         θ¹(1-n[i])n[i+1] +
//         θ⁰(1-n[i])(1-n[i+1]) +
//         c  σ+σ- +
//         c* σ-σ+
//
//      and n[i] = 0.5(σz[i] + 1) acts on site i.
//      Alternatively, we could write this in terms of σz operators, as
//         θ³/4 (1+σz[i])(1+σz[i+1]) +
//         θ²/4 (1+σz[i])(1-σz[i+1]) +
//         θ¹/4 (1-σz[i])(1+σz[i+1]) +
//         θ⁰/4 (1-σz[i])(1-σz[i+1]) +
//         c  σ+σ- +
//         c* σ-σ+
//
//     */
//
//     tools::log->trace("Generating twosite unitaries");
//     auto                  SZ        = qm::spin::half::gen_twobody_spins(qm::spin::half::sz); // We use these as matrices
//     auto                  SP        = qm::spin::half::gen_twobody_spins(qm::spin::half::sp); // We use these as matrices
//     auto                  SM        = qm::spin::half::gen_twobody_spins(qm::spin::half::sm); // We use these as matrices
//     auto                  ID        = qm::spin::half::gen_twobody_spins(qm::spin::half::id); // We use these as matrices
//     auto                  N         = std::vector<Eigen::Matrix4cd>{0.5 * (ID[0] + SZ[0]), 0.5 * (ID[1] + SZ[1])};
//     auto                  spin_dims = std::vector<long>{2l, 2l};
//     std::vector<qm::Gate> unitaries;
//     unitaries.reserve(sites - 1);
//     for(size_t idx = 0; idx < sites - 1; idx++) {
//         double               th0 = rnd::normal(0, 1);
//         double               th1 = rnd::normal(0, 1);
//         double               th2 = rnd::normal(0, 1);
//         double               th3 = rnd::normal(0, 1);
//         std::complex<double> c(rnd::normal(0, 1), rnd::normal(0, 1));
//         auto                 indices = std::vector<size_t>{idx, idx + 1};
//         Eigen::Matrix4cd     H       = th3 * N[0] * N[1] + th2 * N[1] * (ID[0] - N[0]) + th1 * N[0] * (ID[1] - N[1]) + th0 * (ID[0] - N[0]) * (ID[1] - N[1])
//         +
//                              SP[0] * SM[1] * c + SP[1] * SM[0] * std::conj(c);
//
//         // Here we shuffle to get the correct underlying index pattern: Sites are contracted left-to right, but
//         // the kronecker product that generated two-site gates above has indexed right-to-left
//         //         0                   1      0              0      1                0
//         //         |                   |      |              |      |                |
//         //   [ exp(-ifM) ]  --->    [ exp(-ifM) ]   --->  [ exp(-ifM) ]  --->  [ exp(-ifM) ]
//         //         |                   |      |              |      |                |
//         //         1                   3      2              2      3                1
//         Eigen::Tensor<cplx, 2> M_shuffled = tenx::TensorMap(H, 2, 2, 2, 2).shuffle(tenx::array4{1, 0, 3, 2}).reshape(tenx::array2{4, 4});
//         Eigen::MatrixXcd       expifM     = (imn * fmix * tenx::MatrixMap(M_shuffled)).exp();
//         unitaries.emplace_back(tenx::TensorMap(expifM), indices, spin_dims);
//     }
//     if constexpr(settings::debug) {
//         // Sanity check
//         for(const auto &u : unitaries)
//             if(not tenx::MatrixMap(u.op).isUnitary()) throw except::logic_error("u is not unitary!");
//     }
//
//     return unitaries;
// }

std::vector<qm::Gate> qm::lbit::get_unitary_2gate_layer(const qm::lbit::UnitaryGateProperties &u) {
    /*! Returns a set of unitary two site operators used to transform between physical and l-bit representations
     *
     *
     @verbatim
                   0      1                0
                   |      |                |
                 [ exp(-ifM) ]  ---> [ exp(-ifM) ]
                   |      |               |
                   2      3               1
     @endverbatim

     Where M is defined as
        θ³n[i]n[i+1] +
        θ²n[i](1-n[i+1]) +
        θ¹(1-n[i])n[i+1] +
        θ⁰(1-n[i])(1-n[i+1]) +
        c  σ+σ- +
        c* σ-σ+

     and n[i] = 0.5(σz[i] + 1) acts on site i.
     Alternatively, we could write this in terms of σz operators, as
        θ⁰/4 (1+σz[i])(1+σz[i+1]) +
        θ¹/4 (1+σz[i])(1-σz[i+1]) +
        θ²/4 (1-σz[i])(1+σz[i+1]) +
        θ³/4 (1-σz[i])(1-σz[i+1]) +
        c  σ+σ- +
        c* σ-σ+
    */

    tools::log->trace("Generating twosite unitaries");
    if(u.tgw8 == UnitaryGateWeight::EXPDECAY or u.cgw8 == UnitaryGateWeight::EXPDECAY) {
        if(u.hvals.empty() or u.hvals.size() != u.sites) throw except::logic_error("uprop.hvals.size() {} != sites {}", u.hvals.size(), u.sites);
    }
    constexpr bool        kroneckerSwap = false;
    auto                  SZ            = qm::spin::half::gen_twobody_spins(qm::spin::half::sz, kroneckerSwap); // We use these as matrices
    auto                  SP            = qm::spin::half::gen_twobody_spins(qm::spin::half::sp, kroneckerSwap); // We use these as matrices
    auto                  SM            = qm::spin::half::gen_twobody_spins(qm::spin::half::sm, kroneckerSwap); // We use these as matrices
    auto                  ID            = qm::spin::half::gen_twobody_spins(qm::spin::half::id, kroneckerSwap); // We use these as matrices
    auto                  N             = std::vector<Eigen::Matrix4cd>{0.5 * (ID[0] + SZ[0]), 0.5 * (ID[1] + SZ[1])};
    auto                  spin_dims     = std::vector<long>{2l, 2l};
    std::vector<qm::Gate> gates;
    gates.reserve(u.sites - 1);
    for(size_t idx = 0; idx < u.sites - 1; idx++) {
        // This 2-site gate connects sites idx and idx+1
        // EXPDECAY correspond to squared exponential decay with adjacent field differences
        auto tw = u.tgw8 == UnitaryGateWeight::EXPDECAY ? std::exp(-2.0 * std::abs(u.hvals[idx] - u.hvals[idx + 1])) : 1.0;
        auto cw = u.cgw8 == UnitaryGateWeight::EXPDECAY ? std::exp(-2.0 * std::abs(u.hvals[idx] - u.hvals[idx + 1])) : 1.0;

        auto             th0 = tw * rnd::normal(0, u.tstd);
        auto             th1 = tw * rnd::normal(0, u.tstd);
        auto             th2 = tw * rnd::normal(0, u.tstd);
        auto             th3 = tw * rnd::normal(0, u.tstd);
        auto             c   = cw * std::complex<double>(rnd::normal(0, u.cstd), rnd::normal(0, u.cstd)); // (weighted) complex normal random variable
        Eigen::Matrix4cd M   = th3 * N[0] * N[1] + th2 * N[1] * (ID[0] - N[0]) + th1 * N[0] * (ID[1] - N[1]) + th0 * (ID[0] - N[0]) * (ID[1] - N[1]) +
                             SP[0] * SM[1] * c + SP[1] * SM[0] * std::conj(c);

        // Here we shuffle to get the correct underlying index pattern: Sites are contracted left-to right, but
        // the kronecker product that generated two-site gates above has indexed right-to-left
        //         0                   1      0              0      1                0
        //         |                   |      |              |      |                |
        //   [ exp(-ifM) ]  --->    [ exp(-ifM) ]   --->  [ exp(-ifM) ]  --->  [ exp(-ifM) ]
        //         |                   |      |              |      |                |
        //         1                   3      2              2      3                1
        Eigen::Tensor<cplx, 2> M_shuffled = tenx::TensorMap(M, 2, 2, 2, 2).shuffle(tenx::array4{1, 0, 3, 2}).reshape(tenx::array2{4, 4});
        Eigen::MatrixXcd       expifM     = (imn * u.fmix * tenx::MatrixMap(M_shuffled)).exp();
        gates.emplace_back(tenx::TensorMap(expifM), std::vector<size_t>{idx, idx + 1}, spin_dims);
    }
    if constexpr(settings::debug) {
        // Sanity check
        for(const auto &g : gates)
            if(not tenx::MatrixMap(g.op).isUnitary()) throw except::logic_error("u is not unitary!");
    }

    return gates;
}

/*! \brief Make MPO's out of a layer of unitary 2-site gates.
 *
 * Step 1: Reshape the first unitary gate, introducing a dummy 0 index.
 *         The numbers 01 denotes the site indices in the gate.
 * @verbatim
 *          0     1            1     2          1     3
 *          |     |            |     |          |     |
 *         [ g(01) ]   =   0--[ g(01) ]  =  0--[ g(01) ]
 *          |     |            |     |          |     |
 *          2     3            3     4          2     4
 * @endverbatim
 *
 * Step 2: Make a matrix merging indices (012,34) and split the gate using SVD into an MPO and a gate:
 * @verbatim
 *            1     3                                   2               1
 *            |     |   SVD                             |               |
 *        0--[ g(01) ]   =  US^0.5 S^0.5V^T  =  0--[ mpo(0) ]--1   0--[g(1)]
 *            |     |                                   |               |
 *            2     4                                   3               2
 * @endverbatim
 *
 * Step 3: Connect g(1) to the next gate. Connect from under if g is odd, or from above if g is even:
 * @verbatim
 *                0     1
 *                |     |
 *               [ g(12) ]             2     3               1     3
 *                |     |              |     |               |     |
 *                2     3      =   0--[ g(12) ]     =    0--[ g(12) ]
 *                1                    |     |               |     |
 *                |                    1     4               2     4
 *          0--[g(1)]
 *                |
 *                2
 * or
 *                1
 *                |
 *          0--[g(1)]
 *                |                    1     2               1     3
 *                2                    |     |               |     |
 *                0     1      =   0--[ g(12) ]     =    0--[ g(12) ]
 *                |     |              |     |               |     |
 *               [ g(12) ]             3     4               2     4
 *                |     |
 *                2     3
 * @endverbatim
 *
 * Repeat from step 2, replacing g(01) with g(12), and so on, until there are no more gates.
 *
 * Step 4: Add a dummy index with dimension 1 to the remaining right-most gate, to convert it to an mpo.
 */
std::vector<Eigen::Tensor<cplx, 4>> qm::lbit::get_unitary_mpo_layer(const std::vector<qm::Gate> &ulayer) {
    if(ulayer.empty()) return {};
    if(ulayer.size() <= 1) throw except::logic_error("Can't make MPO's with less than two gates");
    auto mpos = std::vector<Eigen::Tensor<cplx, 4>>(ulayer.size() + 1); // L-1 gates in the unitary circuit
    auto svd  = svd::solver();
    // Start Step 1 with the first gate
    const auto &ugate0 = ulayer.front();                                                  // The left-most unitary gate
    auto        dims4  = ugate0.shape<4>();                                               // Original shape of the gate
    auto        dims5  = std::array<long, 5>{1, dims4[0], dims4[1], dims4[2], dims4[3]};  // Shape including the dummy index
    auto        shuf5  = std::array<long, 5>{0, 1, 3, 2, 4};                              // Shuffling to prepare for SVD
    auto        gate5  = Eigen::Tensor<cplx, 5>(ugate0.op.reshape(dims5).shuffle(shuf5)); // The gate including dummy index
    auto        gate3  = Eigen::Tensor<cplx, 3>();                                        // Temporary used in step 2

    // There are L-1 gates on a single layer of the unitary circuit
    for(size_t gidx = 0; gidx < ulayer.size() /* == L-1 */; ++gidx) {
        // Step 2
        std::tie(mpos[gidx], gate3) = svd.split_mpo_gate(gate5);
        if(gidx + 1 >= ulayer.size()) break; // We can't add more gates. Now gate 3 has the mpo corresponding to site L
        // Step 3
        const auto &ugate = ulayer[gidx + 1];
        dims4             = ugate.shape<4>(); // Original shape of the gate
        if((gidx % 2) == 0 /* even */) {
            gate5 = gate3.contract(ugate.op.reshape(dims4), tenx::idx({1}, {2})).shuffle(std::array<long, 5>{0, 2, 1, 3, 4});
        } else {
            gate5 = gate3.contract(ugate.op.reshape(dims4), tenx::idx({2}, {0})).shuffle(std::array<long, 5>{0, 1, 3, 2, 4});
        }
    }
    if(gate3.size() != 0) {
        // Step 4:
        mpos.back() = gate3.reshape(std::array<long, 4>{gate3.dimension(0), 1, gate3.dimension(1), gate3.dimension(2)});
    }

    for(const auto &[i, mpo] : iter::enumerate(mpos)) tools::log->info("mpo {}: {}", i, mpo.dimensions());
    return mpos;
}

std::vector<Eigen::Tensor<cplx, 4>> qm::lbit::get_unitary_mpo_layer(const UnitaryGateProperties &u) {
    return get_unitary_mpo_layer(get_unitary_2gate_layer(u));
}

/*! \brief Merge two MPO layers into a single one using SVD.
 *
 * Step 1:
 * @verbatim
 *                                        2                    2
 *                                        |                    |
 *        1        4              0--[ mpo_up_L ]--1   0--[ mpo_up_R ]--1
 *        |        |                      |                    |
 *   0--[ mpo_updn_L ]--3   =             3                    3
 *        |        |                      2                    2
 *        2        5                      |                    |
 *                               0--[ mpo_dn_L ]--1   0--[ mpo_dn_R ]--1
 *                                        |                    |
 *                                        3                    3
 *
 * @endverbatim
 *
 * Step 2:
 * @verbatim
 *          1       4                       2                    2
 *          |       |       SVD             |                    |
 *     0--[ mpo_up_LR ]--3   =      0--[ mpo(i) ]--1   0--[ mpo_updn_tmp ]--1
 *          |       |                       |                    |
 *          2       5                       3                    3
 * @endverbatim
 *
 * Step 3:
 * @verbatim
 *
 *
 *                                           2
 *                                           |
 *                2                  0--[ mpo_up_R ]--1
 *                |                          |
 *      0--[ mpo_updn_tmp ]--1               3
 *                |                          2
 *                3                          |
 *                                   0--[ mpo_dn_R ]--1
 *                                           |
 *                                           3
 *
 * @endverbatim
 */
std::vector<Eigen::Tensor<cplx, 4>> qm::lbit::merge_unitary_mpo_layers(const std::vector<Eigen::Tensor<cplx, 4>> &mpos_dn,
                                                                       const std::vector<Eigen::Tensor<cplx, 4>> &mpos_up) {
    if(mpos_dn.size() != mpos_up.size()) except::logic_error("size mismatch: {} != {}", mpos_dn.size(), mpos_up.size());
    tools::log->info("Merging mpos");
    auto mpos     = std::vector<Eigen::Tensor<cplx, 4>>(mpos_dn.size());
    auto ddn      = mpos_dn[0].dimensions();
    auto dup      = mpos_up[0].dimensions();
    auto shf6     = std::array<long, 6>{0, 3, 1, 4, 5, 2};
    auto rsh4     = std::array<long, 4>{ddn[0] * dup[0], ddn[1] * dup[1], dup[2], ddn[3]};
    auto mpo_updn = Eigen::Tensor<cplx, 4>(mpos_dn[0].contract(mpos_up[0], tenx::idx({2}, {3})).shuffle(shf6).reshape(rsh4));

    auto cfg           = svd::config();
    cfg.rank_max       = 128;
    cfg.truncation_lim = 1e-12;
    cfg.svd_lib        = svd::lib::lapacke;
    cfg.switchsize_bdc = 16;
    cfg.use_bdc        = true;
    auto   svd         = svd::solver(cfg);
    auto   S           = Eigen::Tensor<cplx, 1>();
    auto   VT          = Eigen::Tensor<cplx, 2>();
    auto   Rm          = Eigen::MatrixXcd(); // The R in matrix form
    auto   Qm          = Eigen::MatrixXcd(); // The Q in matrix form
    auto   qr          = Eigen::ColPivHouseholderQR<Eigen::MatrixXcd>();
    auto   mpo_updn_R  = Eigen::Tensor<cplx, 4>();
    double t_conR      = 0;
    double t_svd       = 0;
    double t_conSVR    = 0;
    for(size_t idx = 0; idx < mpos.size() - 1; ++idx) {
        auto dud = mpo_updn.dimensions();
        {
            auto t     = tid::ur();
            ddn        = mpos_dn[idx + 1].dimensions();
            dup        = mpos_up[idx + 1].dimensions();
            rsh4       = std::array<long, 4>{ddn[0] * dup[0], ddn[1] * dup[1], dup[2], ddn[3]};
            mpo_updn_R = Eigen::Tensor<cplx, 4>(mpos_dn[idx + 1].contract(mpos_up[idx + 1], tenx::idx({2}, {3})).shuffle(shf6).reshape(rsh4));
            t_conR     = t.get_last_interval();
        }
        {
            auto t                     = tid::ur();
            std::tie(mpos[idx], S, VT) = svd.split_mpo_l2r(mpo_updn, cfg);
            t_svd                      = t.get_last_interval();
            mpo_updn                   = tenx::asDiagonal(S).contract(VT, tenx::idx({1}, {0})).contract(mpo_updn_R, tenx::idx({1}, {0}));
            t_conSVR                   = t.get_last_interval() - t_svd;
        }
        tools::log->info("split svd mpo {}: {} --> mpo {} + S {} VT {} | mpoR {} | trunc {:.4e} | time svd {:.3e} conR {:.3e} + conSVR {:.3e} = tot {:.3e}",
                         idx, dud, mpos[idx].dimensions(), S.dimensions(), VT.dimensions(), mpo_updn_R.dimensions(), svd.get_truncation_error(), t_svd, t_conR,
                         t_conSVR, t_svd + t_conR);
    }
    mpos.back() = mpo_updn;

    return mpos;
}

/*! \brief Merge multiple MPO layers into a single one using SVD.
 */
std::vector<Eigen::Tensor<cplx, 4>> qm::lbit::merge_unitary_mpo_layers(const std::vector<std::vector<Eigen::Tensor<cplx, 4>>> &mpo_layers) {
    auto mpo_merged = std::vector<Eigen::Tensor<cplx, 4>>(mpo_layers.front());
    for(size_t idx = 1; idx < mpo_layers.size(); ++idx) { mpo_merged = merge_unitary_mpo_layers(mpo_merged, mpo_layers[idx]); }
    return mpo_merged;
}

Eigen::Tensor<cplx, 2> qm::lbit::get_time_evolution_operator(cplx delta_t, const Eigen::Tensor<cplx, 2> &H) {
    // Given a matrix H, this returns exp(delta_t * H)
    // For time evolution, just make sure delta_t = -i*d,  where d is a (small) real positive number.
    return tenx::TensorCast((delta_t * tenx::MatrixMap(H)).exp());
}

std::vector<qm::Gate> qm::lbit::get_time_evolution_gates(cplx delta_t, const std::vector<qm::Gate> &hams_nsite, double id_threshold) {
    std::vector<Gate> time_evolution_gates;
    time_evolution_gates.reserve(hams_nsite.size());
    for(auto &h : hams_nsite) {
        time_evolution_gates.emplace_back(h.exp(imn * delta_t)); // exp(-i * delta_t * h)
        if(tenx::isIdentity(time_evolution_gates.back().op, id_threshold)) {
            tools::log->trace("get_time_evolution_gates: ignoring time evolution swap gate {} == I +- {:.2e}", time_evolution_gates.back().pos, id_threshold);
            time_evolution_gates.pop_back(); // Skip this gate if it is just an identity.
        }
    }
    if constexpr(settings::debug) {
        for(auto &t : time_evolution_gates)
            if(not t.isUnitary(Eigen::NumTraits<double>::dummy_precision() * static_cast<double>(t.op.dimension(0)))) {
                throw except::runtime_error("Time evolution operator at pos {} is not unitary:\n{}", t.pos, linalg::tensor::to_string(t.op));
            }
    }
    time_evolution_gates.shrink_to_fit();
    return time_evolution_gates;
}

std::vector<qm::SwapGate> qm::lbit::get_time_evolution_swap_gates(cplx delta_t, const std::vector<qm::SwapGate> &hams_nsite, double id_threshold) {
    std::vector<SwapGate> time_evolution_swap_gates;
    time_evolution_swap_gates.reserve(hams_nsite.size());
    size_t count_ignored = 0;
    for(auto &h : hams_nsite) {
        time_evolution_swap_gates.emplace_back(h.exp(imn * delta_t)); // exp(-i * delta_t * h)
        if(tenx::isIdentity(time_evolution_swap_gates.back().op, id_threshold)) {
            count_ignored++;
            tools::log->trace("get_time_evolution_swap_gates: ignoring time evolution swap gate {} == I +- {:.2e}", time_evolution_swap_gates.back().pos,
                              id_threshold);
            time_evolution_swap_gates.pop_back(); // Skip this gate if it is just an identity.
        } else {
            time_evolution_swap_gates.back().generate_swap_sequences();
        }
    }
    if constexpr(settings::debug) {
        for(auto &t : time_evolution_swap_gates)
            if(not t.isUnitary(Eigen::NumTraits<double>::dummy_precision() * static_cast<double>(t.op.dimension(0)))) {
                throw except::runtime_error("Time evolution operator at pos {} is not unitary:\n{}", t.pos, linalg::tensor::to_string(t.op));
            }
    }
    time_evolution_swap_gates.shrink_to_fit();
    tools::log->debug("get_time_evolution_swap_gates: ignored {} time evolution swap gates: I +- {:.2e}", count_ignored, id_threshold);
    return time_evolution_swap_gates;
}

std::vector<Eigen::Tensor<cplx, 2>> qm::lbit::get_time_evolution_operators_2site(size_t sites, cplx delta_t,
                                                                                 const std::vector<Eigen::Tensor<cplx, 2>> &hams_2site) {
    // In l-bit systems we are aldready in a diagonal basis, so h_{j,j+1} and h_{j+1,j+2} commute. Therefore we can immediately use the relation
    //      exp(-i*dt *[h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}]) =  exp(-i*dt [h_{j,j+1}]) * exp(-i*dt*[h_{j+1,j+2}]) * ... * exp(-i*dt*[h_{L-2, L-1}])
    // without passing through the Suzuki-Trotter decomposition.
    // Here we expect "hams_2site" to contain terms like  h_{j,j+1} + h_{j+1,j+2} + ... + h_{L-2, L-1}.

    if(hams_2site.size() != sites - 1) throw except::logic_error("Wrong number of twosite hamiltonians: {}. Expected {}", hams_2site.size(), sites - 1);

    std::vector<Eigen::Tensor<cplx, 2>> time_evolution_operators;
    time_evolution_operators.reserve(sites - 1);
    for(const auto &h : hams_2site) time_evolution_operators.emplace_back(get_time_evolution_operator(delta_t, h));
    return time_evolution_operators;
}

std::vector<Eigen::Tensor<cplx, 2>> qm::lbit::get_time_evolution_operators_3site(size_t sites, cplx delta_t,
                                                                                 const std::vector<Eigen::Tensor<cplx, 2>> &hams_3site) {
    // In l-bit systems we are aldready in a diagonal basis, so h_{i,j,k} and h_{l,m,n} commute. Therefore we can immediately use the relation
    // exp(A + B) = exp(A)exp(B)
    // without passing through the Suzuki-Trotter decomposition.

    if(hams_3site.size() != sites - 2) throw except::logic_error("Wrong number of three-site hamiltonians: {}. Expected {}", hams_3site.size(), sites - 2);

    std::vector<Eigen::Tensor<cplx, 2>> time_evolution_operators;
    time_evolution_operators.reserve(sites - 1);
    for(const auto &h : hams_3site) time_evolution_operators.emplace_back(get_time_evolution_operator(delta_t, h));
    return time_evolution_operators;
}

qm::cplx qm::lbit::get_lbit_exp_value3(const std::vector<std::vector<qm::Gate>> &unitary_layers, const Eigen::Matrix2cd &szi, size_t pos_szi,
                                       const Eigen::Matrix2cd &szj, size_t pos_szj, long sites) {
    /*! \brief Calculates the operator overlap O(i,j) = Tr(𝜌'_i σ^z_j) / Tr(𝜌'_i) = Tr(τ^z_i σ^z_j) / 2^L

        Consider first the density matrix

             ρ_i = 1/2^L  (σ^z_i + 1 ) ⊗ I_i',

        where I_i' is the identity matrix on sites j != i.
        Then ρ_i represents a state |psi><psi| where site i is at level 0 with probability 1,
        and the others sites are in a cat-state of level 0 or 1. E.g. ρ could be a state with
        magnetization 1 at site i, and 0 elsewhere, or ρ could be a state with 1 particle at
        site i, and cat state of 0 and 1 particles elsewhere.

        Applying the unitary transformation does not change the trace, since it is just a
        change of basis. We get

             ρ'_i = U† ρ_i U
                  = 1/2^L (U† σ^z_i U + 1)
                  = 1/2^L (τ^z_i + 1) ,

        where the l-bit τ^z_i = U† σ^z_i U acts non-trivially on all sites. Carrying out the
        trace gives us

            Tr(ρ'_i σ^z_j) / Tr(ρ'_i) = 1/2^L Tr(τ^z_i σ^z_j).

        where we used Tr(1) = 2^L, Tr(ρ_i) = Tr(ρ'_i) = 1 and Tr(τ^z_i) = 0.

        We expect an l-bit fully localized at site i to give O(i,i) = 1,  and a fully
        delocalized l-bit to give O(i,j) = 0.5.

        Note that when the light-cone from site i can't reach site j in the unitary circuit, then
        τ^z_i has no support where σ^z_j connects.  Therefore, we get

               Tr(τ^z_i σ^z_j) = Tr(τ^z_i ⊗ σ^z_j) = Tr(τ^z_i) Tr(σ^z_j) = 0

        Read more here:
        https://link.aps.org/doi/10.1103/PhysRevB.91.085425
        https://onlinelibrary.wiley.com/doi/10.1002/andp.201600322
    */

    // Generate gates for the operators
    if constexpr(settings::debug_circuit) tools::log->trace("Computing Tr (τ_{} σ_{}) / Tr(σ_{})", pos_szi, pos_szj, pos_szj);
    auto t_olap = tid::ur();
    t_olap.tic();

    auto result          = cplx(0, 0);
    auto szi_gate        = qm::Gate(szi, {pos_szi}, {2l}); //
    auto szj_gate        = qm::Gate(szj, {pos_szj}, {2l});
    auto intersection    = qm::get_lightcone_intersection(unitary_layers, pos_szi, pos_szj);
    auto unitary_slayers = qm::get_lightcone_gate_selection(unitary_layers, intersection, false); // Selected gates in each layer
    auto is_disconnected = std::any_of(intersection.begin(), intersection.end(), [](auto &layer) { return layer.empty(); });
    if(is_disconnected) {
        if constexpr(settings::debug_circuit) tools::log->trace("σzi:{} and σzj:{} are disconnected -> result = {:.6f}", pos_szi, pos_szj, result);
        return result;
    }
    // Setup debug printing of unitary circuits
    size_t                  uw        = 7;                                // Width of a unitary 2-site gate box
    size_t                  op        = 3;                                // Overlap of a unitary 2-site gate box
    size_t                  tw        = 6;                                // Tag width
    size_t                  mp        = static_cast<size_t>(sites) - 1ul; // Max pos
    auto                    empty_str = fmt::format("{0:^{1}}", " ", tw + mp * (uw - op) + op);
    std::deque<std::string> net(2 + unitary_slayers.size(), empty_str); // Great for debugging
    std::deque<std::string> log(2 + unitary_slayers.size());            // Great for debugging

    auto all_empty = [](const auto &slayers) {
        return std::all_of(slayers.begin(), slayers.end(), [](const auto &slayer) -> bool { return slayer.empty(); });
    };
    auto num_gates_L = [](const auto &slayers, size_t idx = 0) {
        size_t num = 0;
        for(size_t i = idx + 1; i < slayers.size(); i++) {
            if(slayers.at(i - 1).empty()) return num;
            if(slayers.at(i).empty()) return num;
            size_t posLL = slayers.at(i).front().pos.front();
            size_t posL  = slayers.at(i - 1).front().pos.front();
            if(posLL + 1 == posL)
                num++;
            else
                break;
        }
        return num;
    };
    auto num_gates_R = [](const auto &slayers, size_t idx = 0) {
        size_t num = 0;
        for(size_t i = idx + 1; i < slayers.size(); i++) {
            if(slayers.at(i - 1).empty()) return num;
            if(slayers.at(i).empty()) return num;
            size_t posRR = slayers.at(i).back().pos.back();
            size_t posR  = slayers.at(i - 1).back().pos.back();
            if(posR + 1 == posRR)
                num++;
            else
                break;
        }
        return num;
    };
    auto slayer_width = [&]() -> size_t {
        size_t width = 0;
        for(const auto &s : unitary_slayers) {
            size_t w = 0;
            for(const auto &l : s) w += l.pos.size();
            width = std::max(width, w);
        }
        return width;
    };
    auto slayer_diagl = [&]() -> size_t {
        size_t lb = -1ul; // left bottom position
        size_t rb = -1ul; // right bottom position
        size_t lt = -1ul; // left  top position
        size_t rt = -1ul; // right top position
        for(const auto &us : unitary_slayers) {
            if(not us.empty()) {
                if(lb == -1ul) lb = us.front().pos.front();
                if(rb == -1ul) rb = us.front().pos.back();
                lb = std::min(lb, us.front().pos.front());
            }
        }
        for(const auto &us : iter::reverse(unitary_slayers)) {
            if(not us.empty()) {
                if(rt == -1ul) rt = us.front().pos.back();
                if(lt == -1ul) lt = us.front().pos.front();
                rt = std::max(rt, us.back().pos.back());
            }
        }
        return std::max(rb - lb, rt - lt) + 1;
    };

    auto slayer_diagr = [&]() -> size_t {
        size_t lb = -1ul; // left bottom position
        size_t rb = -1ul; // right bottom position
        size_t lt = -1ul; // left  top position
        size_t rt = -1ul; // right top position
        for(const auto &us : unitary_slayers) {
            if(not us.empty()) {
                if(lb == -1ul) lb = us.back().pos.front();
                if(rb == -1ul) rb = us.back().pos.back();
                rb = std::max(rb, us.back().pos.back());
            }
        }
        for(const auto &us : iter::reverse(unitary_slayers)) {
            if(not us.empty()) {
                if(rt == -1ul) rt = us.back().pos.back();
                if(lt == -1ul) lt = us.back().pos.front();
                lt = std::min(lt, us.front().pos.front());
            }
        }
        return std::max(rb - lb, rt - lt) + 1;
    };

    auto width = slayer_width();
    auto diagr = slayer_diagr();
    auto diagl = slayer_diagl();

    bool go_diagonal = std::min(diagl, diagr) < width; // Avoid long thin diagonals.

    // Start contracting selected unitary gates bottom to top.
    auto g = szi_gate; // The gate that accumulates everything. Starts at σzi position
    if constexpr(settings::debug_circuit) {
        szi_gate.draw_pos(net.front(), "szi  :");
        log.front().append(fmt::format("insert σzi{}", szi_gate.pos));
    }
    if(go_diagonal) {
        // Decide to take the left or right diagonal (whatever is cheapest)
        bool go_right = diagr <= diagl;
        while(not all_empty(unitary_slayers)) {
            for(auto &&[idx_slayer, slayer] : iter::enumerate(unitary_slayers)) {
                if(slayer.empty()) continue; // Already applied the whole slayer
                auto &layer_str = net.at(idx_slayer + 1);
                auto &story_str = log.at(idx_slayer + 1);
                auto  numL      = num_gates_L(unitary_slayers, idx_slayer); // number of gates remaining to take going left
                auto  numR      = num_gates_R(unitary_slayers, idx_slayer); // number of gates remaining to take going right

                std::vector<size_t> pos_out;
                if(go_right) {
                    //                    tools::log->trace("-> insert u[{}]:{}", idx_slayer, slayer.back().pos);
                    if constexpr(settings::debug_circuit) {
                        slayer.back().draw_pos(layer_str, fmt::format("u[{:^2}]:", idx_slayer));
                        story_str.append(fmt::format("insert u{} ", slayer.back().pos));
                    }
                    g       = g.insert(slayer.back());
                    pos_out = slayer.back().pos_difference(intersection.at(idx_slayer + 1));
                    slayer.pop_back();
                } else {
                    //                    tools::log->trace("<- insert u[{}]:{}", idx_slayer, slayer.front().pos);
                    if constexpr(settings::debug_circuit) {
                        slayer.front().draw_pos(layer_str, fmt::format("u[{:^2}]:", idx_slayer));
                        story_str.append(fmt::format("insert u{} ", slayer.front().pos));
                    }
                    g       = g.insert(slayer.front());
                    pos_out = slayer.front().pos_difference(intersection.at(idx_slayer + 1));
                    slayer.pop_front();
                }

                // Collect positions that could be traced
                if(not pos_out.empty()) {                                                      // Trace positions outside the light cone intersection
                    pos_out.erase(std::unique(pos_out.begin(), pos_out.end()), pos_out.end()); // Keep unique elements
                    //                    tools::log->trace("trace[{}]:{}", idx_slayer, pos_out);
                    g    = g.trace_pos(pos_out);
                    g.op = g.op / g.op.constant(std::pow(2, pos_out.size())); // Normalize by dividing the trace of each 2x2 identity.
                    if constexpr(settings::debug_circuit) story_str.append(fmt::format("trace{} ", pos_out));
                }

                if((numR == 0 and go_right) or (numL == 0 and not go_right)) {
                    if constexpr(settings::debug_circuit) story_str.append(fmt::format("break {} ", g.pos));
                    break; // Start from the bottom of the circuit
                }
            }
        }
    } else {
        for(auto &&[idx_slayer, slayer] : iter::enumerate(unitary_slayers)) {
            auto &layer_str = net.at(idx_slayer + 1);
            auto &story_str = log.at(idx_slayer + 1);
            for(auto &sgate : slayer) {
                if(g.has_pos(sgate.pos)) {
                    g = g.insert(sgate);
                    //                    tools::log->trace("insert u[{}]:{}", idx_slayer, sgate.pos);
                    if constexpr(settings::debug_circuit) {
                        sgate.draw_pos(layer_str, fmt::format("u[{:^2}]:", idx_slayer));
                        story_str.append(fmt::format("insert u{} ", sgate.pos));
                    }
                    auto pos_out = sgate.pos_difference(intersection.at(idx_slayer + 1));
                    if(not pos_out.empty()) {
                        // Trace positions outside the light cone intersection
                        g    = g.trace_pos(pos_out);
                        g.op = g.op / g.op.constant(std::pow(2, pos_out.size())); // Normalize by dividing the trace of each 2x2 identity.
                        if constexpr(settings::debug_circuit) story_str.append(fmt::format("trace{} ", pos_out));
                    }
                }
            }
            if constexpr(settings::debug_circuit) story_str.append(fmt::format("now{} ", g.pos));
        }
    }

    if(g.has_pos(szj_gate.pos)) {
        // Connect σ^z_j at the top
        szj_gate.draw_pos(net.back(), "szj  :");
        log.back().append(fmt::format("insert σzj{} ", szj_gate.pos));
        g = szj_gate.connect_above(g);
    }
    // In the last step we trace everything down to a cplx
    result = g.trace();
    result /= std::pow(2, g.pos.size()); // Normalize by dividing the trace of each 2x2 identity.
    t_olap.toc();
    if constexpr(settings::debug_circuit) {
        log.back().append(fmt::format("result = {:.6f} | time {:.3e}", result, t_olap.get_time()));
        for(const auto &[idx, layer] : iter::enumerate_reverse(net)) tools::log->debug("{} | {}", layer, log[idx]);
    }
    return result;
}

qm::cplx qm::lbit::get_lbit_exp_value4(const std::vector<Eigen::Tensor<cplx, 4>> &mpo_layer, const Eigen::Matrix2cd &szi, size_t pos_szi,
                                       const Eigen::Matrix2cd &szj, size_t pos_szj) {
    /*! \brief Calculates the operator overlap O(i,j) = Tr(𝜌'_i σ^z_j) / Tr(𝜌'_i) = Tr(τ^z_i σ^z_j) / 2^L

        Consider first the density matrix

             ρ_i = 1/2^L  (σ^z_i + 1 ) ⊗ I_i',

        where I_i' is the identity matrix on sites j != i.
        Then ρ_i represents a state |psi><psi| where site i is at level 0 with probability 1,
        and the others sites are in a cat-state of level 0 or 1. E.g. ρ could be a state with
        magnetization 1 at site i, and 0 elsewhere, or ρ could be a state with 1 particle at
        site i, and cat state of 0 and 1 particles elsewhere.

        Applying the unitary transformation does not change the trace, since it is just a
        change of basis. We get

             ρ'_i = U† ρ_i U
                  = 1/2^L (U† σ^z_i U + 1)
                  = 1/2^L (τ^z_i + 1) ,

        where the l-bit τ^z_i = U† σ^z_i U acts non-trivially on all sites. Carrying out the
        trace gives us

            Tr(ρ'_i σ^z_j) / Tr(ρ'_i) = 1/2^L Tr(τ^z_i σ^z_j).

        where we used Tr(1) = 2^L, Tr(ρ_i) = Tr(ρ'_i) = 1 and Tr(τ^z_i) = 0.

        We expect an l-bit fully localized at site i to give O(i,i) = 1,  and a fully
        delocalized l-bit to give O(i,j) = 0.5.

        Note that when the light-cone from site i can't reach site j in the unitary circuit, then
        τ^z_i has no support where σ^z_j connects.  Therefore, we get

               Tr(τ^z_i σ^z_j) = Tr(τ^z_i ⊗ σ^z_j) = Tr(τ^z_i) Tr(σ^z_j) = 0

        Read more here:
        https://link.aps.org/doi/10.1103/PhysRevB.91.085425
        https://onlinelibrary.wiley.com/doi/10.1002/andp.201600322


        In this implementation we have contracted the unitary circuit into a series of MPOs.

        Step 1: Contract the top operator (I or sigma_j) on top of the upper mpo

        @verbatim
                          0
                          |
                        [opj]
                          |                       3
                          1                       |
                          2          =    0--[ mpo_up_opj ]--1
                          |                       |
                   0--[ mpo_up ]--1               2
                          |
                          3

        @endverbatim

        Step 2: Similarly, contract the middle operator (I or sigma_i) on top of the bottom mpo^dagger
        @verbatim
                           0                                      0
                           |                                      |
                         [opi]                                  [opi]
                           |                                      |                            3                                   2
                           1                   shuffle            1                            |              shuffle              |
                  [        2         ]^dagger    =                3             =     0--[ mpo_dn_opi  ]--1     =         0--[ mpo_dn_opi  ]--1
                  [        |         ]                            |                            |                                   |
                  [ 0--[ mpo_dn ]--1 ]                     0--[ mpo_dn*]--1                    2                                   3
                  [        |         ]                            |
                  [        3         ]                            2

        We can save two shuffles by contracting on the opposite spin index immediately instead
                           0
                           |
                         [opi]
                           |                                      2
                           1 ------------                         |
                           2            |                0--[ mpo_dn_opi  ]--1
                           |            |       =                 |
                    0--[ mpo_dn* ]--1   |                         3
                           |            |
                           3 <-----------

        Then contract this object upside down using leg 3 when conecting to mpo_up_opj.
        @endverbatim

     Step 3: Contract like


        |



    */

    // Generate gates for the operators
    if constexpr(settings::debug_circuit) tools::log->trace("Computing Tr (τ_{} σ_{}) / Tr(σ_{})", pos_szi, pos_szj, pos_szj);
    auto t_olap = tid::ur();
    t_olap.tic();

    auto id         = tenx::TensorCast(qm::spin::half::id);
    auto si         = tenx::TensorCast(szi);
    auto sj         = tenx::TensorCast(szj);
    auto result     = Eigen::Tensor<cplx, 2>(1, 1);
    auto temp       = Eigen::Tensor<cplx, 2>();
    auto mpo_dn_opi = Eigen::Tensor<cplx, 4>();
    auto mpo_up_opj = Eigen::Tensor<cplx, 4>();
    result.setConstant(1.0);
    for(size_t idx = 0; idx < mpo_layer.size(); ++idx) {
        auto opi   = idx == pos_szi ? si : id;
        auto opj   = idx == pos_szj ? sj : id;
        mpo_up_opj = Eigen::Tensor<cplx, 4>(mpo_layer[idx].contract(opj, tenx::idx({2}, {1}))); //.shuffle(std::array<long, 4>{0, 1, 3, 2}));
        mpo_dn_opi = Eigen::Tensor<cplx, 4>(mpo_layer[idx].conjugate().contract(opi, tenx::idx({3}, {1})));

        // Start contracting with the result
        temp   = result.contract(mpo_dn_opi, tenx::idx({0}, {0})).contract(mpo_up_opj, tenx::idx({0, 3}, {0, 2})).trace(std::array<long, 2>{1, 3});
        result = temp / temp.constant(2.0); // Divide by two for each trace
    }
    //    tools::log->info("result {:+.16f}{:+.16f}i | time tot {:.3e}", std::real(result.coeff(0)), std::imag(result.coeff(0)), t_olap.get_last_interval());
    return result.coeff(0);
}

qm::cplx qm::lbit::get_lbit_correlator(StateFinite &state1, StateFinite &state2, const Eigen::Matrix2cd &szi, size_t pos_szi, const Eigen::Matrix2cd &szj,
                                       size_t pos_szj, long len) {
    state1.get_mps_site(pos_szi).apply_mpo(tenx::TensorCast(szi));
    state2.get_mps_site(pos_szj).apply_mpo(tenx::TensorCast(szj));
    auto overlap = tools::finite::ops::overlap(state1, state2);
    tools::log->info("overlap²: {:.3e}", std::pow(overlap, 2));
    return std::pow(overlap, 2);
}

Eigen::Tensor<qm::cplx, 2> qm::lbit::get_lbit_support(const std::vector<std::vector<qm::Gate>> &unitary_layers, size_t sites) {
    /*! \brief Calculates the operator overlap O(i,j) = Tr(𝜌_i σ^z_j) / Tr(𝜌_j)
                                                      = Tr(τ^z_i σ^z_j) / 2^L
                                                      = Tr(U†σ^z_iU σ^z_j) / 2^L
        Where
            𝜌_i   = 1/2^L (1 + τ^z_i)   a density matrix describing a state localized around site i
            τ^z_i = U† σ^z_i  U,        is an l-bit operator at site i
            U                           is a unitary transformation in the form of a finite-depth circuit.
        The second equality comes from the fact that Tr(τ^z_i) = 0 (l-bit operators are traceless) and Tr(1) = 2^L.
        An l-bit fully localized at site i gives O(i,i) = 1,  and when fully delocalized one gets O(i,j) = 0.5

        Read more here:
        https://link.aps.org/doi/10.1103/PhysRevB.91.085425
        https://onlinelibrary.wiley.com/doi/10.1002/andp.201600322
    */
    auto ssites       = static_cast<long>(sites);
    auto lbit_overlap = Eigen::Tensor<qm::cplx, 2>(ssites, ssites);
    auto state        = StateFinite(AlgorithmType::fLBIT, sites, 0, 2);
    auto pauli        = qm::spin::half::sz;
    tools::finite::mps::init::set_product_state_aligned(state, StateInitType::REAL, "+x");
    for(const auto &layer : unitary_layers) { tools::finite::mps::apply_gates(state, layer, true, GateMove::AUTO); }

#pragma omp parallel for collapse(2) schedule(guided, 4)
    for(long j = 0; j < ssites; j++) {
        for(long i = 0; i < ssites; i++) {
            auto state_j = StateFinite(state);
            auto state_i = StateFinite(state);
            lbit_overlap(i, j) =
                //                qm::lbit::get_lbit_exp_value3(unitary_layers, qm::spin::half::sz, static_cast<size_t>(i), qm::spin::half::sz,
                //                static_cast<size_t>(j), ssites);
                qm::lbit::get_lbit_correlator(state_i, state_j, pauli, static_cast<size_t>(i), pauli, static_cast<size_t>(j), ssites);
        }
    }
    // We require that lbit_overlap(i,j) has rows that sum up to 1
    auto sums_rowwise = tenx::MatrixMap(lbit_overlap).rowwise().sum();
    if(not sums_rowwise.cwiseAbs().isOnes(1e-4)) {
        tools::log->error("lbit overlap rows do not sum to one. Perhaps normalization is wrong.\n"
                          "lbit_overlap: \n{}\nsums\n{}\n",
                          linalg::tensor::to_string(lbit_overlap.real(), 16), linalg::matrix::to_string(sums_rowwise, 6));
        //        throw except::logic_error("lbit overlap rows do not sum to one. Perhaps normalization is wrong");
    }
    return lbit_overlap;
}

Eigen::Tensor<qm::cplx, 2> qm::lbit::get_lbit_support(const std::vector<Eigen::Tensor<cplx, 4>> &mpo_layer) {
    auto                       ssites = static_cast<long>(mpo_layer.size());
    Eigen::Tensor<qm::cplx, 2> lbit_overlap(ssites, ssites);

    /*! \brief Calculates the operator overlap O(i,j) = Tr(𝜌_i σ^z_j) / Tr(𝜌_j)
                                                      = Tr(τ^z_i σ^z_j) / 2^L
                                                      = Tr(U†σ^z_iU σ^z_j) / 2^L
        Where
            𝜌_i   = 1/2^L (1 + τ^z_i)   a density matrix describing a state localized around site i
            τ^z_i = U† σ^z_i  U,        is an l-bit operator at site i
            U                           is a unitary transformation in the form of a finite-depth circuit.
        The second equality comes from the fact that Tr(τ^z_i) = 0 (l-bit operators are traceless) and Tr(1) = 2^L.
        An l-bit fully localized at site i gives O(i,i) = 1,  and when fully delocalized one gets O(i,j) = 0.5

        Read more here:
        https://link.aps.org/doi/10.1103/PhysRevB.91.085425
        https://onlinelibrary.wiley.com/doi/10.1002/andp.201600322
    */
#pragma omp parallel for collapse(2) schedule(guided, 4)
    for(long j = 0; j < ssites; j++) {
        for(long i = 0; i < ssites; i++) {
            lbit_overlap(i, j) =
                qm::lbit::get_lbit_exp_value4(mpo_layer, qm::spin::half::sz, static_cast<size_t>(i), qm::spin::half::sz, static_cast<size_t>(j));
        }
    }
    // We require that lbit_overlap(i,j) has rows that sum up to 1
    auto sums_rowwise = tenx::MatrixMap(lbit_overlap).rowwise().sum();
    if(not sums_rowwise.cwiseAbs().isOnes(1e-4)) {
        tools::log->error("lbit overlap rows do not sum to one. Perhaps normalization is wrong\n"
                          "lbit_overlap: \n{}\nsums\n{}\n",
                          linalg::tensor::to_string(lbit_overlap, 16), linalg::matrix::to_string(sums_rowwise, 6));
        //        throw except::logic_error("lbit overlap rows do not sum to one. Perhaps normalization is wrong");
    }
    return lbit_overlap;
}

std::tuple<double, double, std::vector<double>, size_t> qm::lbit::get_characteristic_length_scale(const Eigen::Tensor<double, 2> &lbit_permuted_disorder_avg) {
    // Average along each column to get an estimate of the lbit
    Eigen::Tensor<double, 1> lbit_permuted_disorder_avg_site_avg = lbit_permuted_disorder_avg.real().mean(std::array<long, 1>{0}); // Average over each l-bit
    Eigen::Tensor<double, 1> lbit_permuted_disorder_avg_site_avg_abs_log =
        lbit_permuted_disorder_avg_site_avg.abs().log(); // Abs Log for fitting                                 // Make all values positive, and take log
    //    tools::log->info("lbit_permuted_disorder_avg_site_avg    : \n{}\n", linalg::tensor::to_string(lbit_permuted_disorder_avg_site_avg, 18, 20));
    //    tools::log->info("lbit_permuted_disorder_avg_site_avg_log: \n{}\n", linalg::tensor::to_string(lbit_permuted_disorder_avg_site_avg_abs_log, 18,
    //    20));

    // Data becomes noisy if the exponential has decayed, so find a cutoff to get the slope using only the first part of the curve
    auto y      = std::vector<double>(lbit_permuted_disorder_avg_site_avg.data(),
                                 lbit_permuted_disorder_avg_site_avg.data() + lbit_permuted_disorder_avg_site_avg.size());
    auto v      = stat::find_last_valid_point(y);
    auto c      = std::count_if(y.begin(), y.begin() + static_cast<long>(v), [](auto &val) { return std::abs(val) > std::numeric_limits<double>::epsilon(); });
    auto x_full = num::range<double>(0, c);
    auto x_tail = num::range<double>(c / 2, c);
    auto y_full = tenx::span(lbit_permuted_disorder_avg_site_avg.data(), c);
    auto y_tail = tenx::span(y_full.data() + c / 2, y_full.data() + c);
    auto ylog_full                = tenx::span(lbit_permuted_disorder_avg_site_avg_abs_log.data(), c);
    auto ylog_tail                = tenx::span(ylog_full.data() + c / 2, ylog_full.data() + c);
    auto [slope_full, res_full]   = stat::slope(x_full, ylog_full);
    auto [slope_tail, res_tail]   = stat::slope(x_tail, ylog_tail);
    double cls_full               = 1.0 / std::abs(slope_full);
    double cls_tail               = 1.0 / std::abs(slope_tail);
    auto   fit_log_stretched_full = fit::log_stretched(x_full, ylog_full, {0.9, 1.0, 1.0});
    auto   fit_log_stretched_tail = fit::log_stretched(x_tail, ylog_tail, {0.9, 1.0, 1.0});
    auto   fit_log_full           = fit::log(x_full, ylog_full, {1.0, 1.0});
    auto   fit_log_tail           = fit::log(x_tail, ylog_tail, {1.0, 1.0});
    auto   fit_exp_stretched_full = fit::exp_stretched(x_full, y_full, {0.9, 1.0, 1.0});
    auto   fit_exp_stretched_tail = fit::exp_stretched(x_tail, y_tail, {0.9, 1.0, 1.0});
    auto   fit_exp_full           = fit::exp(x_full, y_full, {1.0, 1.0});
    auto   fit_exp_tail           = fit::exp(x_tail, y_tail, {1.0, 1.0});

    tools::log->info("Computed fit full slope    cls    {:>8.6f} | sse {:>8.6f} | using {} points: {::8.6f}", cls_full, res_full, ylog_full.size(), ylog_full);
    tools::log->info("Computed fit tail slope    cls    {:>8.6f} | sse {:>8.6f} | using {} points: {::8.6f}", cls_tail, res_tail, ylog_tail.size(), ylog_tail);
    //    tools::log->info("Computed fit full log_stretched coeffs {::>8.6f} | status {}", fit_log_stretched_full.coeffs, fit_log_stretched_full.status);
    //    tools::log->info("Computed fit tail log_stretched coeffs {::>8.6f} | status {}", fit_log_stretched_tail.coeffs, fit_log_stretched_tail.status);
    //    tools::log->info("Computed fit full log           coeffs {::>8.6f} | status {}", fit_log_full.coeffs, fit_log_full.status);
    //    tools::log->info("Computed fit tail log           coeffs {::>8.6f} | status {}", fit_log_tail.coeffs, fit_log_tail.status);
    //    tools::log->info("Computed fit full exp_stretched coeffs {::>8.6f} | status {}", fit_exp_stretched_full.coeffs, fit_exp_stretched_full.status);
    //    tools::log->info("Computed fit tail exp_stretched coeffs {::>8.6f} | status {}", fit_exp_stretched_tail.coeffs, fit_exp_stretched_tail.status);
    tools::log->info("Computed fit full exp           coeffs {::>8.6f} | status {}", fit_exp_full.coeffs, fit_exp_full.status);
    tools::log->info("Computed fit tail exp           coeffs {::>8.6f} | status {}", fit_exp_tail.coeffs, fit_exp_tail.status);
    return {cls_full, res_full, y, c};
}

std::pair<Eigen::Tensor<double, 2>, Eigen::Tensor<double, 2>> qm::lbit::get_lbit_support_stats(const std::vector<Eigen::Tensor<cplx, 2>> &lbit_overlap_vec) {
    Eigen::Tensor<double, 2> avg, err;
    if(not lbit_overlap_vec.empty()) {
        long                rows = lbit_overlap_vec.front().dimension(0); // dim 0
        long                cols = lbit_overlap_vec.front().dimension(1); // dim 1
        size_t              reps = lbit_overlap_vec.size();               // dim 2
        std::vector<double> slice(reps);
        avg.resize(rows, cols);
        err.resize(rows, cols);
        for(long c = 0; c < cols; c++) {
            for(long r = 0; r < rows; r++) {
                for(const auto &[i, elem] : iter::enumerate(lbit_overlap_vec)) {
                    if(std::abs(std::imag(elem(r, c))) > 1e-12)
                        tools::log->warn("elem[{}]({},{}) has imaginary component : |Im({})| > 1e-12", i, r, c, elem(r, c));
                    slice[i] = elem(r, c).real();
                }
                avg(r, c) = stat::mean(slice);
                err(r, c) = stat::sterr(slice);
            }
        }
    }
    return {avg, err};
}

template<typename Scalar>
auto get_permuted(const Eigen::Tensor<Scalar, 2> &in) {
    // First, subtract the center position of each lbit, so we get L lbits centered around zero.
    // In practice, we make a cyclic permutation of the rows of lbit_support
    // In addition, we mirror the lbit along its vertical, so that we can average its left and right half together
    long rows = in.dimension(0);
    long cols = in.dimension(1);
    auto out  = Eigen::Tensor<Scalar, 2>(cols, rows);
    out.setZero();
    for(long j = 0; j < cols; j++) {
        for(long i = 0; i < rows; i++) {
            // If i < j, this corresponds to the left side of an lbit.
            // Consider i == j to be "origo" for an lbit, so mod(i+j,cols) becomes the index starting from that origo.
            // In addition, we fold the left side of an lbit back onto the right side.
            // To get the correct average, add just half of the value;
            long j_twin   = j;
            long distance = std::abs(i - j);
            long j_perm   = distance;
            if(j >= i) {
                j_twin = i - distance;
                if(j_twin < 0) j_twin = j;
            } else {
                j_twin = i + distance;
                if(j_twin >= cols) j_twin = j;
            }
            out(i, j_perm) = 0.5 * (in(i, j) + in(i, j_twin));
        }
    }
    return out;
}

std::pair<Eigen::Tensor<double, 2>, Eigen::Tensor<double, 2>> qm::lbit::get_lbit_permute_stats(const std::vector<Eigen::Tensor<cplx, 2>> &lbit_support_vec) {
    auto lbit_permute_vec = std::vector<Eigen::Tensor<cplx, 2>>();
    for(const auto &lbit_support : lbit_support_vec) lbit_permute_vec.emplace_back(get_permuted(lbit_support));
    return get_lbit_support_stats(lbit_permute_vec);
}

std::vector<Eigen::Tensor<cplx, 2>> qm::lbit::get_lbit_supports(const UnitaryGateProperties &uprop, size_t reps, bool randomize_fields) {
    auto lbit_supp_vec = std::vector<Eigen::Tensor<cplx, 2>>(reps);
    for(auto &lbit_supp : lbit_supp_vec) {
        if(randomize_fields) { uprop.randomize_hvals(); }
        auto ulayers = std::vector<std::vector<qm::Gate>>();
        for(size_t idx = 0; idx < uprop.depth; ++idx) { ulayers.emplace_back(qm::lbit::get_unitary_2gate_layer(uprop)); }
        lbit_supp = qm::lbit::get_lbit_support(ulayers, uprop.sites);
    }
    return lbit_supp_vec;
}

/* clang-format off */
qm::lbit::lbitSupportAnalysis qm::lbit::get_lbit_support_analysis(const UnitaryGateProperties      u_defaults,
                                                                  size_t reps,
                                                                  bool                             randomize_fields,
                                                                  std::vector<size_t            >  u_depths,
                                                                  std::vector<double            >  u_fmixs,
                                                                  std::vector<double            >  u_tstds,
                                                                  std::vector<double            >  u_cstds,
                                                                  std::vector<UnitaryGateWeight >  u_tgw8s,
                                                                  std::vector<UnitaryGateWeight >  u_cgw8s) {
    auto t_lbit_analysis = tid::tic_scope("lbit_analysis");
    if(u_depths.empty()) u_depths = {u_defaults.depth};
    if(u_fmixs.empty())  u_fmixs  = {u_defaults.fmix};
    if(u_tstds.empty())  u_tstds  = {u_defaults.tstd};
    if(u_cstds.empty())  u_cstds  = {u_defaults.cstd};
    if(u_tgw8s.empty())  u_tgw8s  = {u_defaults.tgw8};
    if(u_cgw8s.empty())  u_cgw8s  = {u_defaults.cgw8};

    lbitSupportAnalysis lbitSA(u_depths.size(), u_fmixs.size(), u_tstds.size(), u_cstds.size(),u_tgw8s.size(),u_cgw8s.size() , reps, u_defaults.sites);
    std::array<long, 7> offset7{}, extent7{};
    std::array<long, 9> offset9{}, extent9{};
    auto i_width = static_cast<long>(u_defaults.sites);
    extent9 = {1, 1, 1, 1, 1, 1, 1 , i_width, i_width};

    for (const auto & [i_depth, u_depth] : iter::enumerate<long>(u_depths))
    for (const auto & [i_fmix, u_fmix]   : iter::enumerate<long>(u_fmixs))
    for (const auto & [i_cstd, u_cstd]   : iter::enumerate<long>(u_cstds))
    for (const auto & [i_tstd, u_tstd]   : iter::enumerate<long>(u_tstds))
    for (const auto & [i_tgw8, u_tgw8]   : iter::enumerate<long>(u_tgw8s))
    for (const auto & [i_cgw8, u_cgw8]   : iter::enumerate<long>(u_cgw8s))
    {
        auto uprop = u_defaults;
        uprop.depth = u_depth;
        uprop.fmix = u_fmix;
        uprop.tstd = u_tstd;
        uprop.cstd = u_cstd;
        uprop.tgw8 = u_tgw8;
        uprop.cgw8 = u_cgw8;
        tools::log->info("Computing lbit supports | rand h {} | {}", randomize_fields, uprop.string());
        auto lbit_support_vec = get_lbit_supports(uprop, reps, randomize_fields);
        for (const auto & [i_reps, lbit_supp] : iter::enumerate<long>(lbit_support_vec)){
            offset9 = {i_depth, i_fmix, i_tstd, i_cstd, i_tgw8, i_cgw8, i_reps , 0, 0};
            lbitSA.support.slice(offset9, extent9) = lbit_supp.real().reshape(extent9);
            lbitSA.permute.slice(offset9, extent9) = get_permuted(lbit_supp).real().reshape(extent9);
        }
//        auto [lbit_support_avg, lbit_support_err] = qm::lbit::get_lbit_support_stats(lbit_support_vec);
        auto [lbit_permute_avg, lbit_permute_err] = qm::lbit::get_lbit_permute_stats(lbit_support_vec);
        auto [cls, sse, y, c] = qm::lbit::get_characteristic_length_scale(lbit_permute_avg);
        tools::log->info("Computed lbit {} | threads {} | time {:8.3f} s | cls {:>8.6f} | sse {:>8.6f} | decay {:2} sites: {::8.2e}",
                         uprop.string(), omp_get_max_threads(), t_lbit_analysis->restart_lap(),cls, sse, c, y);

        lbitSA.cls_avg(i_depth, i_fmix, i_tstd, i_cstd, i_tgw8, i_cgw8) = cls;
        lbitSA.sse_avg(i_depth, i_fmix, i_tstd, i_cstd, i_tgw8, i_cgw8) = sse;
        offset7                              = {i_depth, i_fmix, i_tstd, i_cstd,i_tgw8, i_cgw8, 0};
        extent7                              = {1, 1, 1, 1, 1, 1, static_cast<long>(y.size())};
        lbitSA.decay_avg.slice(offset7, extent7) = Eigen::TensorMap<Eigen::Tensor<double, 7>>(y.data(), extent7);
    }
    /* clang-format on */
    return lbitSA;
}
