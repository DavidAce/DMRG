#include "mpo.h"
#include "config/debug.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/float.h"
#include "math/svd.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/mpo/MpoSite.h"
#include "tid/tid.h"

std::pair<Eigen::Tensor<cplx, 4>, Eigen::Tensor<cplx, 4>> tools::finite::mpo::swap_mpo(const Eigen::Tensor<cplx, 4> &mpoL, const Eigen::Tensor<cplx, 4> &mpoR) {
    /* The swap operation takes two neighboring sites
     *
     *            (2)dL                           (2)dR
     *              |                               |
     *   (0)mL---[ mpsL ]---mC(1)       mC(0)---[ mpsR ]---(1)mR
     *              |                               |
     *           (3)dL                            (3)dR
     *
     *
     * and joins them into
     *
     *          (1)dL (4)dR
     *             |    |
     *   (0)mL ---[tensor]--- (3)mR
     *             |    |
     *          (2)dL (5)dR
     *
     *
     * Then we swap the physical indices
     *
     *
     *         (4)dR     (1)dL
     *            |        |
     *             \      /
     *               \  /
     *                /\      <---- Swap
     *              /    \
     *            /        \
     *            |        |
     *   (0)mL---[  tensor  ]---(3)mR
     *            |        |
     *             \      /
     *               \  /
     *                /\      <---- Swap
     *              /    \
     *            /        \
     *            |        |
     *         (5)dR     (2)dL
     *
     *
     * This is accomplished by the shuffle:
     *       shuffle(0, 4, 5, 3, 1, 2);
     *
     *
     * We can now split the tensor back into two MPO's, by first reshaping the tensor into a matrix with merged indices (012) x (345)
     *
     *         (4)dR     (1)dL                             (1)dR                                      (1)dL
     *            |        |                                 |                                          |
     *   (0)mL---[  tensor  ]---(3)mR     --->    (0)mL---[ mpoL ]---mC(3)  (0)---[S]---(1)  mC(0)---[ mpoR ]---(1)mR
     *            |        |                                 |                                          |
     *         (5)dR     (2)dL                            (2)dR                                       (2)dL
     *
     *
     * The square root of S can then be multiplied into both left and right MPO's, on the mC index.
     * The left mpo can be shuffled back to standard form with
     *   mpoL: shuffle(0,3,1,2)
     *
     */

    Eigen::Tensor<cplx, 6> swapped_mpo = mpoL.contract(mpoR, tenx::idx({1}, {0})).shuffle(tenx::array6{0, 4, 5, 3, 1, 2}); // swap
    auto                   svd_cfg     = svd::config();
    svd_cfg.truncation_limit           = 1e-16;
    svd_cfg.svd_lib                    = svd::lib::lapacke;
    svd_cfg.svd_rtn                    = svd::rtn::gejsv;

    auto svd = svd::solver(svd_cfg);
    return svd.split_mpo_pair(swapped_mpo, svd_cfg);
}

void tools::finite::mpo::swap_sites(ModelFinite &model, size_t posL, size_t posR, std::vector<size_t> &sites) {
    auto t_swap = tid::tic_scope("swap");
    if(posR != posL + 1) throw except::logic_error("Expected posR == posL+1. Got posL {} and posR {}", posL, posR);
    if(posR != std::clamp(posR, 0ul, model.get_length() - 1ul)) throw except::logic_error("Expected posR in [0,{}]. Got {}", model.get_length() - 1, posR);
    if(posL != std::clamp(posL, 0ul, model.get_length() - 1ul)) throw except::logic_error("Expected posL in [0,{}]. Got {}", model.get_length() - 1, posL);

    auto &mpo_posL = model.get_mpo(posL);
    auto &mpo_posR = model.get_mpo(posR);

    auto [mpoL, mpoR]   = swap_mpo(mpo_posL.MPO(), mpo_posR.MPO());
    auto [mpoL2, mpoR2] = swap_mpo(mpo_posL.MPO2(), mpo_posR.MPO2());
    tools::log->debug("mpo({}) : {} -> {}", mpo_posL.get_position(), mpo_posL.MPO().dimensions(), mpoL.dimensions());
    tools::log->debug("mpo({}) : {} -> {}", mpo_posR.get_position(), mpo_posR.MPO().dimensions(), mpoR.dimensions());
    tools::log->debug("mpo²({}): {} -> {}", mpo_posL.get_position(), mpo_posL.MPO2().dimensions(), mpoL2.dimensions());
    tools::log->debug("mpo²({}): {} -> {}", mpo_posR.get_position(), mpo_posR.MPO2().dimensions(), mpoR2.dimensions());
    mpo_posL.set_mpo(mpoL);
    mpo_posR.set_mpo(mpoR);
    mpo_posL.set_mpo_squared(mpoL2);
    mpo_posR.set_mpo_squared(mpoR2);

    std::swap(sites[posL], sites[posR]);
}

std::vector<Eigen::Tensor<cplx, 4>> tools::finite::mpo::get_mpos_with_edges(const std::vector<Eigen::Tensor<cplx, 4>> &mpos,
                                                                            const Eigen::Tensor<cplx, 1> &Ledge, const Eigen::Tensor<cplx, 1> &Redge) {
    auto mpos_with_edge = mpos;
    auto &threads = tenx::threads::get();

    /* We can prepend edgeL to the first mpo to reduce the size of subsequent operations.
     * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
     *                        2               2
     *                        |               |
     *    0---[L]---(1)(0)---[M]---1 =  0---[LM]---1
     *                        |               |
     *                        3               3
     */
    const auto &mpoL_src = mpos.front();
    auto       &mpoL_tgt = mpos_with_edge.front();
    mpoL_tgt.resize(tenx::array4{1, mpoL_src.dimension(1), mpoL_src.dimension(2), mpoL_src.dimension(3)});
    mpoL_tgt.device(*threads->dev) = Ledge.reshape(tenx::array2{1, Ledge.size()}).contract(mpoL_src, tenx::idx({1}, {0}));

    /* We can append edgeR to the last mpo to reduce the size of subsequent operations.
     * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
     *         2                              1                       2
     *         |                              |                       |
     *    0---[M]---1  0---[R]---1   =  0---[MR]---3  [shuffle]  0---[MR]---1
     *         |                              |                       |
     *         3                              2                       3
     */
    const auto &mpoR_src = mpos.back();
    auto       &mpoR_tgt = mpos_with_edge.back();
    mpoR_tgt.resize(tenx::array4{mpoR_src.dimension(0), 1, mpoR_src.dimension(2), mpoR_src.dimension(3)});
    mpoR_tgt.device(*threads->dev) =
        mpoR_src.contract(Redge.reshape(tenx::array2{Redge.size(), 1}), tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
    return mpos_with_edge;
}

std::vector<Eigen::Tensor<cplx_t, 4>> tools::finite::mpo::get_mpos_with_edges_t(const std::vector<Eigen::Tensor<cplx_t, 4>> &mpos,
                                                                            const Eigen::Tensor<cplx_t, 1> &Ledge, const Eigen::Tensor<cplx_t, 1> &Redge) {
    auto mpos_with_edge = mpos;
    auto &threads = tenx::threads::get();

    /* We can prepend edgeL to the first mpo to reduce the size of subsequent operations.
     * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
     *                        2               2
     *                        |               |
     *    0---[L]---(1)(0)---[M]---1 =  0---[LM]---1
     *                        |               |
     *                        3               3
     */
    const auto &mpoL_src = mpos.front();
    auto       &mpoL_tgt = mpos_with_edge.front();
    mpoL_tgt.resize(tenx::array4{1, mpoL_src.dimension(1), mpoL_src.dimension(2), mpoL_src.dimension(3)});
    mpoL_tgt.device(*threads->dev) = Ledge.reshape(tenx::array2{1, Ledge.size()}).contract(mpoL_src, tenx::idx({1}, {0}));

    /* We can append edgeR to the last mpo to reduce the size of subsequent operations.
     * Start by converting the edge from a rank1 to a rank2 with a dummy index of size 1:
     *         2                              1                       2
     *         |                              |                       |
     *    0---[M]---1  0---[R]---1   =  0---[MR]---3  [shuffle]  0---[MR]---1
     *         |                              |                       |
     *         3                              2                       3
     */
    const auto &mpoR_src = mpos.back();
    auto       &mpoR_tgt = mpos_with_edge.back();
    mpoR_tgt.resize(tenx::array4{mpoR_src.dimension(0), 1, mpoR_src.dimension(2), mpoR_src.dimension(3)});
    mpoR_tgt.device(*threads->dev) =
        mpoR_src.contract(Redge.reshape(tenx::array2{Redge.size(), 1}), tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
    return mpos_with_edge;
}

std::vector<Eigen::Tensor<cplx, 4>> tools::finite::mpo::get_compressed_mpos(std::vector<Eigen::Tensor<cplx, 4>> mpos) {
    tools::log->trace("Compressing MPOs: {} sites", mpos.size());
    //    std::vector<Eigen::Tensor<cplx, 4>> mpos_compressed = mpos;

    // Setup SVD
    // Here we need a lot of precision:
    //  - Use very low svd threshold
    //  - Force the use of JacobiSVD by setting the switchsize_bdc to something large
    //  - Force the use of Lapacke -- it is more precise than Eigen (I don't know why)
    // Eigen Jacobi becomes ?gesvd (i.e. using QR) with the BLAS backend.
    // See here: https://eigen.tuxfamily.org/bz/show_bug.cgi?id=1732
    auto svd_cfg    = svd::config();
    svd_cfg.svd_lib = svd::lib::lapacke;
    svd_cfg.svd_rtn = svd::rtn::gejsv;

    svd::solver svd(svd_cfg);

    // Print the results
    std::vector<std::string> report;
    if(tools::log->level() == spdlog::level::trace)
        for(const auto &mpo : mpos) report.emplace_back(fmt::format("{}", mpo.dimensions()));

    for(size_t iter = 0; iter < 4; iter++) {
        // Next compress from left to right
        Eigen::Tensor<cplx, 2> T_l2r; // Transfer matrix
        Eigen::Tensor<cplx, 4> T_mpo;
        for(auto &&[idx, mpo] : iter::enumerate(mpos)) {
            auto mpo_sq_dim_old = mpo.dimensions();
            if(T_l2r.size() == 0)
                T_mpo = mpo;
            else
                T_mpo = T_l2r.contract(mpo, tenx::idx({1}, {0}));

            if(idx + 1 == mpos.size()) {
                mpo = T_mpo;
            } else {
                std::tie(mpo, T_l2r) = svd.split_mpo_l2r(T_mpo);
                if(idx + 1 == mpos.size())
                    // The remaining transfer matrix T can be multiplied back into the last MPO from the right
                    mpo = Eigen::Tensor<cplx, 4>(mpo.contract(T_l2r, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2}));
            }
            if constexpr(settings::debug) tools::log->trace("iter {} | idx {} | dim {} -> {}", iter, idx, mpo_sq_dim_old, mpo.dimensions());
        }

        // Now we have done left to right. Next we do right to left
        Eigen::Tensor<cplx, 2> T_r2l; // Transfer matrix
        Eigen::Tensor<cplx, 4> mpo_T; // Absorbs transfer matrix
        for(auto &&[idx, mpo] : iter::enumerate_reverse(mpos)) {
            auto mpo_dim_old = mpo.dimensions();
            if(T_r2l.size() == 0)
                mpo_T = mpo;
            else
                mpo_T = mpo.contract(T_r2l, tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
            if(idx == 0) {
                mpo = mpo_T;
            } else {
                std::tie(T_r2l, mpo) = svd.split_mpo_r2l(mpo_T);
                if(idx == 0)
                    // The remaining transfer matrix T can be multiplied back into the first MPO from the left
                    mpo = Eigen::Tensor<cplx, 4>(T_r2l.contract(mpo, tenx::idx({1}, {0})));
            }
            if constexpr(settings::debug) tools::log->trace("iter {} | idx {} | dim {} -> {}", iter, idx, mpo_dim_old, mpo.dimensions());
        }
    }

    // Print the results
    if(tools::log->level() == spdlog::level::trace)
        for(const auto &[idx, msg] : iter::enumerate(report)) tools::log->debug("mpo {}: {} -> {}", idx, msg, mpos[idx].dimensions());

    return mpos;
}
std::vector<Eigen::Tensor<cplx, 4>> tools::finite::mpo::get_compressed_mpos(const std::vector<Eigen::Tensor<cplx, 4>> &mpos,
                                                                            const Eigen::Tensor<cplx, 1> &Ledge, const Eigen::Tensor<cplx, 1> &Redge) {
    return get_compressed_mpos(get_mpos_with_edges(mpos, Ledge, Redge));
}
