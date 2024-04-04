#include "mpo.h"
#include "config/debug.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/float.h"
#include "math/svd.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/site/mpo/MpoSite.h"
#include "tid/tid.h"

//
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include <general/sfinae.h>

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
    auto  mpos_with_edge = mpos;
    auto &threads        = tenx::threads::get();

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
    mpoR_tgt.device(*threads->dev) = mpoR_src.contract(Redge.reshape(tenx::array2{Redge.size(), 1}), tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
    return mpos_with_edge;
}

std::vector<Eigen::Tensor<cplx_t, 4>> tools::finite::mpo::get_mpos_with_edges_t(const std::vector<Eigen::Tensor<cplx_t, 4>> &mpos,
                                                                                const Eigen::Tensor<cplx_t, 1> &Ledge, const Eigen::Tensor<cplx_t, 1> &Redge) {
    auto  mpos_with_edge = mpos;
    auto &threads        = tenx::threads::get();

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
    mpoR_tgt.device(*threads->dev) = mpoR_src.contract(Redge.reshape(tenx::array2{Redge.size(), 1}), tenx::idx({1}, {0})).shuffle(tenx::array4{0, 3, 1, 2});
    return mpos_with_edge;
}

std::vector<Eigen::Tensor<cplx, 4>> tools::finite::mpo::get_compressed_mpos(std::vector<Eigen::Tensor<cplx, 4>> mpos) {
    tools::log->trace("Compressing MPOs: {} sites", mpos.size());
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
    // if(tools::log->level() == spdlog::level::trace)
    for(const auto &mpo : mpos) report.emplace_back(fmt::format("{}", mpo.dimensions()));

    for(size_t iter = 0; iter < 4; iter++) {
        // Next compress from left to right
        Eigen::Tensor<cplx, 2> T_l2r; // Transfer matrix
        Eigen::Tensor<cplx, 4> T_mpo;
        for(auto &&[idx, mpo] : iter::enumerate(mpos)) {
            auto mpo_dim_old = mpo.dimensions();
            if(T_l2r.size() == 0)
                T_mpo = mpo; // First iter
            else
                T_mpo = T_l2r.contract(mpo, tenx::idx({1}, {0})); // Subsequent iters

            if(idx + 1 == mpos.size()) {
                mpo = T_mpo;
            } else {
                std::tie(mpo, T_l2r) = svd.split_mpo_l2r(T_mpo);
            }
            if constexpr(settings::debug) tools::log->trace("iter {} | idx {} | dim {} -> {}", iter, idx, mpo_dim_old, mpo.dimensions());
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
            }
            if constexpr(settings::debug) tools::log->trace("iter {} | idx {} | dim {} -> {}", iter, idx, mpo_dim_old, mpo.dimensions());
        }
    }

    // Print the results
    // if(tools::log->level() == spdlog::level::trace)
    for(const auto &[idx, msg] : iter::enumerate(report)) tools::log->info("mpo {}: {} -> {}", idx, msg, mpos[idx].dimensions());

    return mpos;
}
std::vector<Eigen::Tensor<cplx, 4>> tools::finite::mpo::get_compressed_mpos(const std::vector<Eigen::Tensor<cplx, 4>> &mpos,
                                                                            const Eigen::Tensor<cplx, 1> &Ledge, const Eigen::Tensor<cplx, 1> &Redge) {
    return get_compressed_mpos(get_mpos_with_edges(mpos, Ledge, Redge));
}

template<typename Scalar>
std::vector<Eigen::Tensor<Scalar, 4>> get_inverted_mpos_initial_guess(const std::vector<Eigen::Tensor<Scalar, 4>> &mpos, long extra_dims = 0) {
    auto impos = mpos;
    for(auto &&[idx, impo] : iter::enumerate(impos)) {
        auto dims = impo.dimensions();
        if(idx == 0)
            dims[1] += extra_dims;
        else if(idx + 1 == impos.size())
            dims[0] += extra_dims;
        else {
            dims[0] += extra_dims;
            dims[1] += extra_dims;
        }
        impo.resize(dims);
        impo.setRandom();
    }
    return impos;
}

template<typename Scalar>
Eigen::Tensor<Scalar, 2> get_M(const std::vector<Eigen::Tensor<Scalar, 4>> &mpos, const std::vector<Eigen::Tensor<Scalar, 4>> &impos, size_t pos) {
    assert(mpos.size() == impos.size());
    assert(mpos.front().dimension(0) == 1);
    assert(mpos.back().dimension(1) == 1);
    assert(impos.front().dimension(0) == 1);
    assert(impos.back().dimension(1) == 1);

    Eigen::Tensor<Scalar, 4> LE(1, 1, 1, 1), RE(1, 1, 1, 1); // Left and right environments
    Eigen::Tensor<Scalar, 4> LE_temp, RE_temp;

    LE.setConstant(1.0);
    RE.setConstant(1.0);
    auto &threads = tenx::threads::get();
    for(size_t idx = 0; idx < pos; ++idx) {
        const auto              &mpo         = mpos[idx];
        const auto              &impo        = impos[idx];
        Eigen::Tensor<Scalar, 4> mpo_dagger  = mpo.conjugate().shuffle(tenx::array4{0, 1, 3, 2});
        Eigen::Tensor<Scalar, 4> impo_dagger = impo.conjugate().shuffle(tenx::array4{0, 1, 3, 2});
        LE_temp.resize(impo_dagger.dimension(1), mpo_dagger.dimension(1), mpo.dimension(1), impo.dimension(1));
        LE_temp.device(*threads->dev) = LE.contract(impo_dagger, tenx::idx({0}, {0}))
                                           .contract(mpo_dagger, tenx::idx({0, 5}, {0, 2}))
                                           .contract(mpo, tenx::idx({0, 5}, {0, 2}))
                                           .contract(impo, tenx::idx({0, 5, 2}, {0, 2, 3}));
        LE = std::move(LE_temp);
    }
    for(size_t idx = mpos.size() - 1; idx > pos; --idx) {
        const auto              &mpo         = mpos[idx];
        const auto              &impo        = impos[idx];
        Eigen::Tensor<Scalar, 4> mpo_dagger  = mpo.conjugate().shuffle(tenx::array4{0, 1, 3, 2});
        Eigen::Tensor<Scalar, 4> impo_dagger = impo.conjugate().shuffle(tenx::array4{0, 1, 3, 2});
        RE_temp.resize(impo_dagger.dimension(0), mpo_dagger.dimension(0), mpo.dimension(0), impo.dimension(0));
        RE_temp.device(*threads->dev) = RE.contract(impo_dagger, tenx::idx({0}, {1}))
                                           .contract(mpo_dagger, tenx::idx({0, 5}, {1, 2}))
                                           .contract(mpo, tenx::idx({0, 5}, {1, 2}))
                                           .contract(impo, tenx::idx({0, 5, 2}, {1, 2, 3}));
        RE = std::move(RE_temp);
    }

    const auto              &mpo         = mpos[pos];
    const auto              &impo        = impos[pos];
    Eigen::Tensor<Scalar, 4> mpo_dagger  = mpo.conjugate().shuffle(tenx::array4{0, 1, 3, 2});
    Eigen::Tensor<Scalar, 4> impo_dagger = impo.conjugate().shuffle(tenx::array4{0, 1, 3, 2});
    Eigen::Tensor<Scalar, 2> identity    = tenx::TensorIdentity<Scalar>(impo.dimension(3)); // To set up  2 dummy legs with an identity outer product

    auto rows = LE.dimension(0) * RE.dimension(0) * impo_dagger.dimension(2) * impo_dagger.dimension(3); // Should be the product of impo_dagger dims
    auto cols = LE.dimension(3) * RE.dimension(3) * mpo.dimension(2) * mpo.dimension(3);                 // Should be the product of mpo dims
    Eigen::Tensor<Scalar, 2> M(rows, cols);
    M.device(*threads->dev) = LE.contract(mpo_dagger, tenx::idx({1}, {0}))
                                 .contract(mpo, tenx::idx({1, 5}, {0, 2}))
                                 .contract(RE, tenx::idx({2, 4}, {1, 2}))
                                 .contract(identity, tenx::idx()) // Add two dummy legs
                                 .shuffle(tenx::array8{0, 4, 6, 2, 1, 5, 3, 7})
                                 .reshape(tenx::array2{rows, cols});
    return M;
}
template<typename Scalar>
Eigen::Tensor<Scalar, 1> get_N(const std::vector<Eigen::Tensor<Scalar, 4>> &mpos, const std::vector<Eigen::Tensor<Scalar, 4>> &impos, size_t pos) {
    assert(mpos.size() == impos.size());
    assert(mpos.front().dimension(0) == 1);
    assert(mpos.back().dimension(1) == 1);
    assert(impos.front().dimension(0) == 1);
    assert(impos.back().dimension(1) == 1);

    Eigen::Tensor<Scalar, 2> LE(1, 1), RE(1, 1); // Left and right environments
    Eigen::Tensor<Scalar, 2> LE_temp, RE_temp;

    LE.setConstant(1.0);
    RE.setConstant(1.0);
    auto &threads = tenx::threads::get();

    for(size_t idx = 0; idx < pos; ++idx) {
        Eigen::Tensor<Scalar, 4> mpo_dagger  = mpos[idx].conjugate().shuffle(tenx::array4{0, 1, 3, 2});
        Eigen::Tensor<Scalar, 4> impo_dagger = impos[idx].conjugate().shuffle(tenx::array4{0, 1, 3, 2});
        LE_temp.resize(impo_dagger.dimension(1), mpo_dagger.dimension(1));
        LE_temp.device(*threads->dev) = LE.contract(impo_dagger, tenx::idx({0}, {0})).contract(mpo_dagger, tenx::idx({0, 3, 2}, {0, 2, 3}));
        LE                           = std::move(LE_temp);
    }
    for(size_t idx = mpos.size() - 1; idx > pos; --idx) {
        Eigen::Tensor<Scalar, 4> mpo_dagger  = mpos[idx].conjugate().shuffle(tenx::array4{0, 1, 3, 2});
        Eigen::Tensor<Scalar, 4> impo_dagger = impos[idx].conjugate().shuffle(tenx::array4{0, 1, 3, 2});
        RE_temp.resize(impo_dagger.dimension(0), mpo_dagger.dimension(0));
        RE_temp.device(*threads->dev) = RE.contract(impo_dagger, tenx::idx({0}, {1})).contract(mpo_dagger, tenx::idx({0, 3, 2}, {1, 2, 3}));
        RE                           = std::move(RE_temp);
    }

    const auto              &mpo         = mpos[pos];
    const auto              &impo        = impos[pos];
    Eigen::Tensor<Scalar, 4> mpo_dagger  = mpo.conjugate().shuffle(tenx::array4{0, 1, 3, 2});
    Eigen::Tensor<Scalar, 4> impo_dagger = impo.conjugate().shuffle(tenx::array4{0, 1, 3, 2});

    auto                     rows = LE.dimension(0) * RE.dimension(0) * mpo_dagger.dimension(2) * mpo_dagger.dimension(3);
    Eigen::Tensor<Scalar, 1> N(rows);
    N.device(*threads->dev) =
        LE.contract(mpo_dagger, tenx::idx({1}, {0})).contract(RE, tenx::idx({1}, {1})).shuffle(tenx::array4{0, 3, 1, 2}).reshape(tenx::array1{rows});
    return N;
}
template<typename Scalar>
Eigen::Tensor<Scalar, 1> get_B(const std::vector<Eigen::Tensor<Scalar, 4>> &impos, size_t pos) {
    // Linearize it without shuffling
    return impos[pos].reshape(tenx::array1{impos[pos].size()}).template cast<Scalar>();
}

template<typename Scalar>
void set_B(const Eigen::Tensor<Scalar, 1> &B, std::vector<Eigen::Tensor<Scalar, 4>> &impos, size_t pos) {
    // Assume that B can be interpreted with the correct shape
    impos[pos] = B.reshape(impos[pos].dimensions());
}
template<typename Scalar>
std::vector<Eigen::Tensor<cplx, 4>> get_inverted_mpos_internal(const std::vector<Eigen::Tensor<Scalar, 4>> &mpos) {
    using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    auto   impos     = get_inverted_mpos_initial_guess(mpos, 12);
    auto   lu        = Eigen::LLT<MatrixType>();
    double error     = 1;
    for(size_t iter = 0; iter < 500; ++iter) {
        // Left to right sweep
        for(size_t pos = 0; pos < mpos.size(); ++pos) {
            auto M = get_M(mpos, impos, pos);
            auto N = get_N(mpos, impos, pos);
            auto B = get_B(impos, pos);

            auto M_map = Eigen::Map<MatrixType>(M.data(), M.dimension(0), M.dimension(1));
            auto N_map = Eigen::Map<VectorType>(N.data(), N.size());
            auto B_map = Eigen::Map<VectorType>(B.data(), B.size());

            lu.compute(M_map);
            B_map = lu.solve(N_map);
            set_B(B, impos, pos);
            Scalar BMB = B_map.adjoint() * M_map * B_map;
            Scalar BN  = B_map.dot(N_map);
            Scalar NB  = N_map.dot(B_map);
            error      = std::abs(BMB - BN - NB + std::pow(2, mpos.size()));
            fmt::print("error L2R: {:.3e} | dims {} | {}\n", error, N.dimensions(), sfinae::type_name<Scalar>());
            if(error < 1e-14) break;
        }
        if(error < 1e-14) break;
        // Right to left sweep
        for(size_t pos = mpos.size() - 1; pos < mpos.size(); --pos) {
            auto M = get_M(mpos, impos, pos);
            auto N = get_N(mpos, impos, pos);
            auto B = get_B(impos, pos);

            auto M_map = Eigen::Map<MatrixType>(M.data(), M.dimension(0), M.dimension(1));
            auto N_map = Eigen::Map<VectorType>(N.data(), N.size());
            auto B_map = Eigen::Map<VectorType>(B.data(), B.size());

            lu.compute(M_map);
            B_map = lu.solve(N_map);
            set_B(B, impos, pos);
            Scalar BMB = B_map.adjoint() * M_map * B_map;
            Scalar BN  = B_map.dot(N_map);
            Scalar NB  = N_map.dot(B_map);
            error      = std::abs(BMB - BN - NB + std::pow(2, mpos.size()));
            fmt::print("error R2L: {:.3e} | dims {} | {}\n", error, N.dimensions(), sfinae::type_name<Scalar>());
            if(error < 1e-14) break;
        }
        if(error < 1e-14) break;
    }
    if constexpr(std::is_same_v<Scalar, cplx>)
        return impos;
    else {
        std::vector<Eigen::Tensor<cplx, 4>> impos_cplx;
        impos_cplx.reserve(impos.size());
        for(const auto &impo : impos) impos_cplx.emplace_back(impo.template cast<cplx>());
        return impos_cplx;
    }
}

std::vector<Eigen::Tensor<cplx, 4>> tools::finite::mpo::get_inverted_mpos(const std::vector<Eigen::Tensor<cplx, 4>> &mpos) {
    bool isComplex = false;
    for(const auto &mpo : mpos) {
        if(!tenx::isReal(mpo)) {
            isComplex = true;
            break;
        }
    }
    if(isComplex) {
        return get_inverted_mpos_internal<cplx>(mpos);
    } else {
        std::vector<Eigen::Tensor<real, 4>> mpos_real;
        mpos_real.reserve(mpos.size());
        for(const auto &mpo : mpos) mpos_real.emplace_back(mpo.real());
        return get_inverted_mpos_internal<real>(mpos_real);
    }
}
