#include "split.h"
#include "contraction.h"
#include "log.h"
#include <config/settings.h>
#include <debug/exceptions.h>
#include <deque>
#include <general/iter.h>
#include <h5pp/h5pp.h>
#include <math/linalg/tensor.h>
#include <math/svd.h>
#include <optional>
#include <tensors/site/mps/MpsSite.h>
#include <tid/tid.h>

namespace settings {
    inline constexpr bool debug_split = false;
}

template<typename Scalar>
std::vector<MpsSite> tools::common::split::split_mps(const Eigen::Tensor<Scalar, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                     const std::vector<size_t> &positions, long center_position, long chi_limit,
                                                     std::optional<svd::settings> svd_settings) {
    /*  Here we split an mps containing multiple sites into its consituent sites.
     *  Consider a case with 3 sites:
     *  spin_dims = {2,2,2}
     *  positions = {6,7,8}
     *  center_position = 6
     *
     *  then
     *
     *  chiL ---[mps]--- chiR
     *            |
     *         d^3=2*2*2
     * Becomes
     *
     *  chiL---[A6]---chi  chi---[LC]---chi   chi---[B7]---chi  chi---[L7]---chi   chi---[B8]---chiR
     *          |                                    |                                    |
     *         d=2                                  d=2                                  d=2
     *
     *
     * If instead we had
     *
     * spin_dims = {2,2,2}
     * positions = {3,4,5}
     * center_position = 4
     *
     * then we would get
     * chiL---[A3]---chi   chi---[L]---chi   chi---[A4]---chi  chi---[LC]---chi   chi---[B5]---chiR
     *         |                                    |                                    |
     *        d=2                                  d=2                                  d=2
     *
     * There is a special case where the multisite tensor has a single site. For instance, if going left-to-right:
     *
     * spin_dims = {2}
     * positions = {3}
     * center_position = 4
     *
     * then we have the following one-site tensor
     * chiL---[A3]---chi    chi---[LC]---chi
     *         |
     *        d=2
     *
     * or, right-to-left:
     *
     * spin_dims = {2}
     * positions = {7}
     * center_position = 6
     *
     * then we have the following one-site tensor
     * chi---[LC]---chi  chiL---[B7]---chi
     *                           |
     *                          d=2
     *
     *
     *
     * The split is done in 2 steps. We start with.
     *
     * (1)chiL ---[mps]--- (2)chiR
     *             |
     *        (0) dL*dR
     *
     * where the index order is in parentheses. This is is reinterpreted as
     *
     *
     * (2)chiL---[    mps    ]---(3)chiR
     *            |         |
     *         (0)dL     (1)dR
     *
     * Here dL=d*d*d... as many as there are A-matrices,
     * and dR=d*d*d*... as many as there are B-matrices.
     *
     * Then the mps is reshuffled into
     *
     * (1)chiL---[    mps    ]---(3)chiR
     *            |         |
     *         (0)dL     (2)dR
     *
     * which then is split up into
     *
     * chiL ---[U]----chi chi---[S]---chi chi---[V]---chiR
     *          |                                |
     *         dL                                dR
     *
     * The U and V tensors contain the sites to the left and right of the center bond,
     * and are sent away for further splitting.
     *
     *
     */
    if constexpr(std::is_same_v<Scalar, cplx>) {
        auto chiL   = multisite_mps.dimension(1);
        auto chiR   = multisite_mps.dimension(2);
        auto t_real = tid::tic_scope("isReal");
        if(chiL * chiR > 512 * 512 and tenx::isReal(multisite_mps)) {
            tools::log->info("Converting to real!");
            return split_mps<real>(multisite_mps.real(), spin_dims, positions, center_position, chi_limit, svd_settings);
        }
    }

    if(positions.empty()) throw std::runtime_error("Could not split multisite tensor: positions list is empty");
    if(multisite_mps.size() == 0) throw std::runtime_error("Could not split multisite tensor: tensor is empty");
    if(spin_dims.empty()) throw std::runtime_error("Could not split multisite tensor: spin_dims list is empty");
    if(spin_dims.size() != positions.size())
        throw std::runtime_error(fmt::format("Could not split multisite tensor: size mismatch in given lists: spin_dims {} != positions {} -- sizes not equal",
                                             spin_dims, positions));
    if(chi_limit <= 0) throw std::runtime_error(fmt::format("Invalid bond dimension limit:  chi_limit = {}", chi_limit));

    auto t_split = tid::tic_scope("split", tid::level::detail);
    // Setup the svd settings if not given explicitly
    if(not svd_settings) svd_settings = svd::settings();
    if(not svd_settings->threshold) svd_settings->threshold = settings::precision::svd_threshold;
    if(not svd_settings->switchsize_bdc) svd_settings->switchsize_bdc = settings::precision::svd_switchsize_bdc;

    // Split the multisite tensor at the given center position.
    std::vector<size_t> positions_left, positions_right;
    std::vector<long>   spin_dims_left, spin_dims_right;

    // Define dimensions and positions after the desired split
    long chiL = multisite_mps.dimension(1);
    long chiR = multisite_mps.dimension(2);
    long dL   = 1;
    long dR   = 1;

    auto pos_it = positions.begin();
    auto dim_it = spin_dims.begin();
    while(pos_it != positions.end() and dim_it != spin_dims.end()) {
        if(num::cmp_less_equal(*pos_it, center_position)) {
            positions_left.emplace_back(*pos_it);
            spin_dims_left.emplace_back(*dim_it);
            dL *= *dim_it;
        } else {
            positions_right.emplace_back(*pos_it);
            spin_dims_right.emplace_back(*dim_it);
            dR *= *dim_it;
        }
        pos_it++;
        dim_it++;
    }

    // Initialize the new mps_sites
    std::vector<MpsSite> mps_sites_As;
    std::deque<MpsSite>  mps_sites_Bs;
    // Initialize the left/right unitary matrices
    Eigen::Tensor<Scalar, 3> U, V;
    Eigen::Tensor<Scalar, 1> S;

    // Now there are 5 options
    // 1) All sites will become A-sites
    // 2) All sites will become B-sites
    // 3) One A and one B-site
    // 4) Mixture of A and B-sites, with more (or eq) A's than B's
    // 5) Mixture of A and B-sites, with more B's than A's

    if constexpr(settings::debug_split)
        tools::log->trace("split: positions {} | L {} | R {} | C {} | spins {} | spinL {} | spinR {} ", positions, positions_left, positions_right,
                          center_position, spin_dims, spin_dims_left, spin_dims_right);
    if(positions.size() == positions_left.size()) {
        if constexpr(settings::debug_split) tools::log->trace("split: option 1");
        auto t_o1    = tid::tic_scope("o1");
        mps_sites_As = internal::split_mps_into_As(multisite_mps, spin_dims_left, positions_left, chi_limit, svd_settings);
        // Remnant stash
        auto &mps     = mps_sites_As.back();
        auto &S_stash = mps.get_S_stash();
        auto  pos     = mps.get_position<long>();
        if(pos == center_position and S_stash) { // Steal LC from stash
            mps.set_LC(S_stash->data, S_stash->error);
            S_stash.reset();
        } else {
            auto &V_stash = mps.get_V_stash();
            if(V_stash) {
                auto vdim = V_stash->data.dimensions();
                if(vdim[0] * vdim[1] != vdim[2])
                    tools::log->error(FMT_STRING("V_stash with dimensions {} for pos {} is not a diagonal matrix!"), vdim, V_stash->pos_dst);
            }
        }

    } else if(positions.size() == positions_right.size()) {
        if constexpr(settings::debug_split) tools::log->trace("split: option 2 - sites {} become B's", positions);
        auto t_o2    = tid::tic_scope("o2");
        mps_sites_Bs = internal::split_mps_into_Bs(multisite_mps, spin_dims_right, positions_right, chi_limit, svd_settings);
        // Take care of stash
        auto &mps     = mps_sites_Bs.front();
        auto &U_stash = mps.get_U_stash();
        auto &S_stash = mps.get_S_stash();
        auto  pos     = mps.get_position<long>();
        if constexpr(settings::debug_split) {
            if(S_stash) {
                if(not tenx::isReal(S_stash->data)) throw except::runtime_error("S_stash is not real:\n{}", linalg::tensor::to_string(S_stash->data, 16, 18));
                if(not tenx::isPositive(S_stash->data))
                    throw except::runtime_error("S_stash is not positive:\n{}", linalg::tensor::to_string(S_stash->data, 16, 18));
            }
        }

        if(pos == center_position + 1 and S_stash) { // The stash in S becomes the LC for the site on the left
            mps.stash_C(S_stash->data, S_stash->error, static_cast<size_t>(center_position));
            S_stash.reset();
        } else if(U_stash) {
            auto udim = U_stash->data.dimensions();
            if(udim[1] != udim[0] * udim[2])
                tools::log->error(FMT_STRING("U_stash with dimensions {} for pos {} is not a diagonal matrix!"), udim, U_stash->pos_dst);
        }

    } else if(positions.size() == 2 and positions_left.size() == 1 and positions_right.size() == 1) {
        if constexpr(settings::debug_split) tools::log->trace("split: option 3");
        auto t_o3 = tid::tic_scope("o3");
        // Set up the SVD
        svd::solver svd(svd_settings);
        std::tie(U, S, V) = svd.schmidt_multisite(multisite_mps, dL, dR, chiL, chiR, chi_limit);
        mps_sites_As.emplace_back(U, std::nullopt, positions.front(), 0, "AC");
        mps_sites_As.back().set_LC(S.template cast<cplx>(), svd.truncation_error);
        mps_sites_Bs.emplace_back(V, std::nullopt, positions.back(), 0, "B");

    } else if(positions_left.size() >= positions_right.size()) {
        if constexpr(settings::debug_split) tools::log->trace("split: option 4");
        auto t_o4    = tid::tic_scope("o4");
        mps_sites_As = internal::split_mps_into_As(multisite_mps, spin_dims_left, positions_left, chi_limit, svd_settings);
        // We expect stashed S and V. Merge these and send onward to B's
        auto &V_stash = mps_sites_As.back().get_V_stash();
        auto &S_stash = mps_sites_As.back().get_S_stash();
        if(V_stash and S_stash) {
            if constexpr(std::is_same_v<Scalar, real>)
                V = tools::common::contraction::contract_bnd_mps_temp(S_stash->data, V_stash->data).real();
            else
                V = tools::common::contraction::contract_bnd_mps_temp(S_stash->data, V_stash->data);
            V_stash.reset();
            S_stash.reset();
            mps_sites_Bs = internal::split_mps_into_Bs(V, spin_dims_right, positions_right, chi_limit, svd_settings);
            // Remnant S is now an LC for the last A
            auto &LC_stash = mps_sites_Bs.front().get_S_stash();
            auto &mpsA     = mps_sites_As.back();
            mpsA.set_LC(LC_stash->data, LC_stash->error);
            LC_stash.reset();
            // Remnant U is merged back into A
            mpsA.take_stash(mps_sites_Bs.front());
        }
    } else if(positions_left.size() < positions_right.size()) {
        if constexpr(settings::debug_split) tools::log->trace("split: option 5");
        auto t_o5    = tid::tic_scope("o5");
        mps_sites_Bs = internal::split_mps_into_Bs(multisite_mps, spin_dims_right, positions_right, chi_limit, svd_settings);
        // We expect stashed U and S. Merge these and send onward to A's
        auto &U_stash = mps_sites_Bs.front().get_U_stash();
        auto &S_stash = mps_sites_Bs.front().get_S_stash();
        if(U_stash and S_stash) {
            if constexpr(std::is_same_v<Scalar, real>)
                U = tools::common::contraction::contract_mps_bnd_temp(U_stash->data, S_stash->data).real();
            else
                U = tools::common::contraction::contract_mps_bnd_temp(U_stash->data, S_stash->data);

            U_stash.reset();
            S_stash.reset();
            mps_sites_As = internal::split_mps_into_As(U, spin_dims_left, positions_left, chi_limit, svd_settings);
            // Remnant S is now an LC for the last A
            auto &LC_stash = mps_sites_As.back().get_S_stash();
            auto &mpsA     = mps_sites_As.back();
            mpsA.set_LC(LC_stash->data, LC_stash->error);
            LC_stash.reset();
            mps_sites_Bs.front().take_stash(mpsA);
        }
    }

    if(mps_sites_As.empty() and mps_sites_Bs.empty()) throw std::runtime_error("Got empty mps_sites from both left and right");

    // Sanity checks
    if(mps_sites_As.size() != spin_dims_left.size())
        throw std::runtime_error(
            fmt::format("Could not split multisite tensor: Got mps_sites_As.size() {} != spins_dims_left {}", mps_sites_As.size(), spin_dims_left.size()));
    if(mps_sites_Bs.size() != spin_dims_right.size())
        throw std::runtime_error(
            fmt::format("Could not split multisite tensor: Got mps_sites_Bs.size() {} != spins_dims_right {}", mps_sites_Bs.size(), spin_dims_right.size()));

    // Move the right side onto the left side (This is equivalent to std::list::splice)
    mps_sites_As.insert(mps_sites_As.end(), std::make_move_iterator(mps_sites_Bs.begin()), std::make_move_iterator(mps_sites_Bs.end()));

    for(const auto &&[idx, mps] : iter::enumerate(mps_sites_As)) {
        auto pos = positions[idx];
        if(not mps.is_at_position(pos)) {
            throw std::runtime_error(fmt::format("Could not split multisite tensor: Position mismatch: expected site {} != mps pos {} | positions to merge {}",
                                                 pos, mps.get_position(), positions));
        }
        if(mps.is_at_position(center_position) and not mps.isCenter()) {
            throw std::runtime_error(
                fmt::format("Could not split multisite tensor: Specified center position {} did not become a center MPS", mps.get_position()));
        }
    }

    //    if(chiL > 500 and chiR > 500){
    //        auto file = h5pp::File("../output/svd-benchmark.h5", h5pp::FilePermission::READWRITE);
    //        std::string dsetbase = "multisite_tensor";
    //        if constexpr (std::is_same_v<Scalar,real>) dsetbase.append("_real");
    //        if constexpr (std::is_same_v<Scalar,cplx>) dsetbase.append("_cplx");
    //        static size_t count = 0;
    //        auto dsetname = fmt::format("{}_{}",dsetbase,count++);
    //        file.writeDataset(multisite_tensor, dsetname, H5D_CHUNKED);
    //        file.writeAttribute(tid::get("split").get_last_interval(), "t_split", dsetname);
    //        file.writeAttribute(spin_dims, "spin_dims", dsetname);
    //        file.writeAttribute(positions, "positions", dsetname);
    //        file.writeAttribute(center_position, "center_position", dsetname);
    //        file.writeAttribute(chi_limit, "chi_limit", dsetname);
    //        file.writeAttribute(t_split->get_last_interval(), "t_split", dsetname);
    //    }

    // Return the all the positions including the bond matrix
    return mps_sites_As;
}

template std::vector<MpsSite> tools::common::split::split_mps(const Eigen::Tensor<real, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                              const std::vector<size_t> &positions, long center_position, long chi_limit,
                                                              std::optional<svd::settings> svd_settings);
template std::vector<MpsSite> tools::common::split::split_mps(const Eigen::Tensor<cplx, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                              const std::vector<size_t> &positions, long center_position, long chi_limit,
                                                              std::optional<svd::settings> svd_settings);

template<typename Scalar>
std::vector<MpsSite> tools::common::split::internal::split_mps_into_As(const Eigen::Tensor<Scalar, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                                       const std::vector<size_t> &positions, long chi_limit,
                                                                       std::optional<svd::settings> svd_settings) {
    /*  Here we split an mps containing multiple sites into its consituent sites from the left.
     *  Consider a case with 3 sites and
     *  spin_dims = {2,2,2}
     *  sites = {6,7,8}
     *
     *  then
     *
     *  chiL ---[mps]--- chiR
     *            |
     *         d^3=2*2*2
     * Becomes
     *
     *  chiL---[A6]---chi  chi---[L7]---chi   chi---[A7]---chi  chi---[L8]---chi   chi---[A8]---chiR
     *          |                                    |                                    |
     *         d=2                                  d=2                                  d=2
     *
     *
     * By convention, a multisite_tensor is created by merging sites from left to right, that is
     *
     * 1-2-3-4-5
     * 12-3-4-5
     * 123-4-5
     * 1234-5
     * 12345
     *
     * Here we do the process in reverse order, i.e.
     *
     * 12345
     * 1234-5
     * 123-4-5
     * 12-3-4-5
     * 1-2-3-4-5
     *
     * Note that the spin dimensions are given in left-to-right order, i.e. 12345
     */

    //    if(spin_dims.empty()) throw std::runtime_error("Could not split multisite tensor from the left: spin_dims list is empty");
    if(spin_dims.size() != positions.size())
        throw std::runtime_error(fmt::format(
            "Could not split multisite tensor from the left: size mismatch in given lists: spin_dims {} != sites {} -- sizes not equal", spin_dims, positions));

    // A special case is when we do one-site tensors. Then we expect
    // this function to receive a "U" without sites in it ( U*S will become a center bond).
    if(positions.empty()) return std::vector<MpsSite>();
    auto t_to_a = tid::tic_scope("to_a");
    // Initialize the resulting container of split sites
    std::vector<MpsSite> mps_sites;

    // Now we have multiple spin dimensions

    // Set up the SVD
    svd::solver svd(svd_settings);

    // Declare the tensors that will catch the schmidt (SVD) decompositions
    Eigen::Tensor<Scalar, 3>                U;                           // This will become the first site to be extracted
    Eigen::Tensor<Scalar, 1>                S;                           // The singular values
    Eigen::Tensor<Scalar, 3>                V = multisite_mps;           // This side contains all the remaining sites
    Eigen::Tensor<Scalar, 3>                SV_temp;                     // Temporary for contracting S*V
    std::optional<Eigen::Tensor<Scalar, 1>> S_prev       = std::nullopt; // Starts out empty, carries the schmidt values from the previous iteration
    double                                  S_prev_error = 0;            // Truncation error from the previous iteration
    for(const auto &[idx, spin_dim] : iter::enumerate(spin_dims)) {
        /* The schmidt decomposition gives us the 3 matrices to the right of the line |:
         *                       |
         *                       |
         *  chi---[S_prev]---chi | chiL ---[U]----chi chi---[S]---chi chi---[V]---chiR
         *                       |          |                                |
         *                       |          d                             d^(sites-1)
         *                       |
         * Here U is an "M" matrix of type A = Lambda * Gamma
         * See the notes on svd.schmidt at its definition
         *
         */

        if(S_prev) {
            auto t_sv = tid::tic_token("contract_sv");
            // Let V absorb S from the previous SVD (see note below)
            V = tools::common::contraction::contract_bnd_mps_temp(S_prev.value(), V, SV_temp);
        }

        if(&spin_dim == &spin_dims.back()) {
            // In the last SVD, it's important that we don't truncate by threshold.
            // If we do, we risk making the V non-diagonal, which would truncate whatever site is on the right, later on.
            svd.threshold = std::numeric_limits<double>::epsilon();
        }

        std::tie(U, S, V) = svd.schmidt_into_left_normalized(V, spin_dim, chi_limit);
        if(S.size() == 0) throw std::runtime_error("Could not split multisite tensor: Got 0 singular values from left svd");

        // Note: Now we have the three components
        // U:   A left-unitary "A" matrix which is a "Lambda * Gamma" in Vidal's notation
        // S:   A set of singular values, "L" matrix, belonging to the site on the right of this one (i.e. the next A to pop out of V)
        // V:   Contains one site less than the it did before. That site is now in U.
        //      Before using V as the next multisite mps, it absorbs S so that the next U to pop out is a left-unitary A-type matrix.

        mps_sites.emplace_back(U, S_prev, positions[idx], S_prev_error, "A"); // S_prev is empty in the first iteration.

        // Store the singular values for the next iteration
        S_prev       = S;
        S_prev_error = svd.truncation_error;
    }

    // Now we have a series of A-A-A-A matrices and their corresponding L's
    // At the last step we have residual S and V left over. Stash them!
    auto &mps = mps_sites.back();
    auto  pos = mps.get_position();
    if(S_prev) mps.stash_S(S_prev.value().template cast<cplx>(), S_prev_error, pos + 1);
    mps.stash_V(V.template cast<cplx>(), pos + 1);
    return mps_sites;
}

template std::vector<MpsSite> tools::common::split::internal::split_mps_into_As(const Eigen::Tensor<real, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                                                const std::vector<size_t> &positions, long chi_limit,
                                                                                std::optional<svd::settings> svd_settings);
template std::vector<MpsSite> tools::common::split::internal::split_mps_into_As(const Eigen::Tensor<cplx, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                                                const std::vector<size_t> &positions, long chi_limit,
                                                                                std::optional<svd::settings> svd_settings);
template<typename Scalar>
std::deque<MpsSite> tools::common::split::internal::split_mps_into_Bs(const Eigen::Tensor<Scalar, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                                      const std::vector<size_t> &positions, long chi_limit,
                                                                      std::optional<svd::settings> svd_settings) {
    /*  Here we split an mps containing multiple sites into its consituent sites from the right
     *  Consider a case with 3 sites and
     *  spin_dims = {2,2,2}
     *  sites = {6,7,8}
     *
     *  then
     *
     *  chiL ---[mps]--- chiR
     *            |
     *         d^3=2*2*2
     * Becomes
     *
     *  chiL---[B6]---chi  chi---[L6]---chi   chi---[B7]---chi  chi---[L7]---chi   chi---[B8]---chiR
     *          |                                    |                                    |
     *         d=2                                  d=2                                  d=2
     *
     *
     * By convention, a multisite_tensor is created by merging sites from left to right, that is
     *
     * 1-2-3-4-5
     * 12-3-4-5
     * 123-4-5
     * 1234-5
     * 12345
     *
     * Here we do the process in reverse order, i.e.
     *
     * 12345
     * 1234-5
     * 123-4-5
     * 12-3-4-5
     * 1-2-3-4-5
     *
     * Note that the spin dimensions are given in left-to-right order, i.e. 12345
     */
    //    if(spin_dims.empty()) throw std::runtime_error("Could not split multisite tensor from the right: spin_dims list is empty");
    if(spin_dims.size() != positions.size())
        throw std::runtime_error(
            fmt::format("Could not split multisite tensor from the right: size mismatch in given lists: spin_dims {} != sites {} -- sizes not equal", spin_dims,
                        positions));

    // A special case is when we do one-site tensors. Then we expect
    // this function to receive a "V^dagger" without sites in it ( S*V^dagger will become a center bond).
    if(positions.empty()) return std::deque<MpsSite>();
    auto t_to_b = tid::tic_scope("to_b");

    // Initialize the resulting container of split sites
    std::deque<MpsSite> mps_sites;
    // Now we have multiple spin dimensions

    // Set up the SVD
    svd::solver svd(svd_settings);

    // Declare the the tensors that will catch the schmidt (SVD) decompositions
    Eigen::Tensor<Scalar, 3>                U = multisite_mps;              // This side contains all the sites
    Eigen::Tensor<Scalar, 1>                S;                              // The singular values
    Eigen::Tensor<Scalar, 3>                V;                              // This will become the first site extracted
    Eigen::Tensor<Scalar, 3>                US_temp;                        // Temporary for contracting U*S
    std::optional<Eigen::Tensor<Scalar, 1>> S_prev       = std::nullopt;    // Starts out empty, needs to be checked outside of this split
    double                                  S_prev_error = 0;               // Truncation error from the previous iteration
    for(const auto &[idx, spin_dim] : iter::enumerate_reverse(spin_dims)) { // Iterate in reverse order
        /* The schmidt decomposition gives us the 3 matrices to the left of the line |:
         *                                                      |
         * chiL ---[U]----chi chi---[S]---chi chi---[V]---chiR  |  chi---[S_prev]---chi
         *          |                                |          |
         *      d^(sites-1)                          d          |
         *
         * Here V is an "M" matrix of type B = Gamma * Lambda
         * See the notes on svd.schmidt at its definition
         */

        if(S_prev) {
            auto t_sv = tid::tic_token("contract_us");
            // Let U absorb S from the previous SVD (see note below)
            U = tools::common::contraction::contract_mps_bnd_temp(U, S_prev.value(), US_temp);
        }

        if(&spin_dim == &spin_dims.front()) {
            // In the last SVD, it's important that we don't truncate by threshold.
            // If we do, we risk making the U non-diagonal, which would truncate whatever site is on the left, later on.
            svd.threshold = std::numeric_limits<double>::epsilon();
        }

        std::tie(U, S, V) = svd.schmidt_into_right_normalized(U, spin_dim, chi_limit);
        if(S.size() == 0) throw std::runtime_error("Could not split multisite tensor: Got 0 singular values from right svd");

        // Note: Now we have the three components
        // U:   Contains one site less than the previous iteration.
        //      Before using U as the next multisite mps, it absorbs the S so that the next V to pop out is a right-unitary B-type matrix.
        // S:   A set of singular values, "L" matrix, belonging to the left site  (i.e. the next B to pop out of U)
        // V:   A right-unitary "B" matrix which is a "Gamma*Lambda" in Vidal's notation

        mps_sites.emplace_front(V, S_prev, positions[idx], S_prev_error, "B"); // S_prev is empty the first iteration.

        // Store the singular values for the next iteration
        S_prev       = S;
        S_prev_error = svd.truncation_error;
    }

    // Now we have a series of B-B-B-B matrices and their corresponding L's
    // At the last step we have residual U and S left over. Stash them!
    auto &mps = mps_sites.front();
    auto  pos = mps.get_position();
    if(pos > 0) { // More left-normalized sites left in U
        if(S_prev) mps.stash_S(S_prev.value().template cast<cplx>(), S_prev_error, pos - 1);
        mps.stash_U(U.template cast<cplx>(), pos - 1);
    }
    return mps_sites;
}

template std::deque<MpsSite> tools::common::split::internal::split_mps_into_Bs(const Eigen::Tensor<real, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                                               const std::vector<size_t> &positions, long chi_limit,
                                                                               std::optional<svd::settings> svd_settings);
template std::deque<MpsSite> tools::common::split::internal::split_mps_into_Bs(const Eigen::Tensor<cplx, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                                               const std::vector<size_t> &positions, long chi_limit,
                                                                               std::optional<svd::settings> svd_settings);