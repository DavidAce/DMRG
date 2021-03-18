#include "split.h"
#include "prof.h"
#include <config/nmspc_settings.h>
#include <deque>
#include <general/nmspc_iter.h>
#include <general/nmspc_tensor_omp.h>
#include <math/svd.h>
#include <optional>
#include <tensors/state/class_mps_site.h>
#include <tools/common/fmt.h>

std::vector<class_mps_site> tools::common::split::split_mps(const Eigen::Tensor<Scalar, 3> &multisite_tensor, const std::vector<long> &spin_dims,
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

    if(positions.empty()) throw std::runtime_error("Could not split multisite tensor: positions list is empty");
    if(multisite_tensor.size() == 0) throw std::runtime_error("Could not split multisite tensor: tensor is empty");
    if(spin_dims.empty()) throw std::runtime_error("Could not split multisite tensor: spin_dims list is empty");
    if(spin_dims.size() != positions.size())
        throw std::runtime_error(fmt::format("Could not split multisite tensor: size mismatch in given lists: spin_dims {} != positions {} -- sizes not equal",
                                             spin_dims, positions));
    if(chi_limit <= 0) throw std::runtime_error(fmt::format("Invalid bond dimension limit:  chi_limit = {}", chi_limit));

    // Setup the svd settings if not given explicitly
    if(not svd_settings) svd_settings = svd::settings();
    if(not svd_settings->threshold) svd_settings->threshold = settings::precision::svd_threshold;
    if(not svd_settings->switchsize) svd_settings->switchsize = settings::precision::svd_switchsize;


    // Split the multisite tensor at the given center position.
    std::vector<size_t> positions_left, positions_right;
    std::vector<long>   spin_dims_left, spin_dims_right;

    // Define dimensions and positions after the desired split
    long chiL = multisite_tensor.dimension(1);
    long chiR = multisite_tensor.dimension(2);
    long dL   = 1;
    long dR   = 1;

    auto pos_it = positions.begin();
    auto dim_it = spin_dims.begin();
    while(pos_it != positions.end() and dim_it != spin_dims.end()) {
        if(static_cast<long>(*pos_it) <= center_position) {
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
    std::vector<class_mps_site> mps_sites_As;
    std::deque<class_mps_site>  mps_sites_Bs;
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
        tools::log->trace("Positions {} | L {} | R {} | C {} | spins {} | spinL {} | spinR {} ", positions, positions_left, positions_right, center_position,
                          spin_dims, spin_dims_left, spin_dims_right);
    if(positions.size() == positions_left.size()) {
        if constexpr(settings::debug_split) tools::log->trace("Option 1");
        auto t_split_svda = tools::common::profile::prof[AlgorithmType::ANY]["t_split_svda"]->tic_token();
        mps_sites_As      = internal::split_mps_into_As(multisite_tensor, spin_dims_left, positions_left, chi_limit, svd_settings);
        // Remnant stash
        auto &mps     = mps_sites_As.back();
        auto &S_stash = mps.get_S_stash();
        auto  pos     = mps.get_position<long>();
        if(pos == center_position and S_stash) { // Steal LC from stash
            mps.set_LC(S_stash->data, S_stash->error);
            S_stash.reset();
        }
    } else if(positions.size() == positions_right.size()) {
        if constexpr(settings::debug_split) tools::log->trace("Option 2 - sites {} become B's", positions);
        auto t_split_svdb = tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdb"]->tic_token();
        mps_sites_Bs      = internal::split_mps_into_Bs(multisite_tensor, spin_dims_right, positions_right, chi_limit, svd_settings);
        // Take care of stash
        auto &mps     = mps_sites_Bs.front();
        auto &S_stash = mps.get_S_stash();
        auto  pos     = mps.get_position<long>();
        if(pos == center_position + 1 and S_stash) { // The stash in S becomes the LC for the site on the left
            mps.stash_C(S_stash->data, S_stash->error, static_cast<size_t>(center_position));
            S_stash.reset();
        }
    } else if(positions.size() == 2 and positions_left.size() == 1 and positions_right.size() == 1) {
        if constexpr(settings::debug_split) tools::log->trace("Option 3");
        auto t_split_svdm = tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdm"]->tic_token();
        auto t_svd        = tools::common::profile::get_default_prof()["t_svd"]->tic_token();
        // Set up the SVD
        svd::solver svd(svd_settings);
        std::tie(U, S, V) = svd.schmidt_multisite(multisite_tensor, dL, dR, chiL, chiR, chi_limit);
        mps_sites_As.emplace_back(U, std::nullopt, positions.front(), 0, "AC");
        mps_sites_As.back().set_LC(S, svd.truncation_error);
        mps_sites_Bs.emplace_back(V, std::nullopt, positions.back(), 0, "B");

        *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_wrk"] += *svd.t_wrk;
        *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_adj"] += *svd.t_adj;
        *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_jac"] += *svd.t_jac;
        *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_svd"] += *svd.t_svd;
    } else if(positions_left.size() >= positions_right.size()) {
        if constexpr(settings::debug_split) tools::log->trace("Option 4");
        auto t_split_svda = tools::common::profile::prof[AlgorithmType::ANY]["t_split_svda"]->tic_token();
        mps_sites_As      = internal::split_mps_into_As(multisite_tensor, spin_dims_left, positions_left, chi_limit, svd_settings);
        // We expect stashed S and V. Merge these and send onward to B's
        auto &V_stash = mps_sites_As.back().get_V_stash();
        auto &S_stash = mps_sites_As.back().get_S_stash();
        if(V_stash and S_stash) {
            V = Textra::asDiagonal(S_stash->data).contract(V_stash->data, Textra::idx({1}, {1})).shuffle(Textra::array3{1, 0, 2});
            V_stash.reset();
            S_stash.reset();
            auto t_split_svdb = tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdb"]->tic_token();
            mps_sites_Bs      = internal::split_mps_into_Bs(V, spin_dims_right, positions_right, chi_limit, svd_settings);
            // Remnant S is now an LC for the last A
            auto &LC_stash = mps_sites_Bs.front().get_S_stash();
            auto &mpsA     = mps_sites_As.back();
            mpsA.set_LC(LC_stash->data, LC_stash->error);
            LC_stash.reset();
            // Remnant U is merged back into A
            mpsA.merge_stash(mps_sites_Bs.front());
        }
    } else if(positions_left.size() < positions_right.size()) {
        if constexpr(settings::debug) tools::log->trace("Option 5");
        auto t_split_svdb = tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdb"]->tic_token();
        mps_sites_Bs      = internal::split_mps_into_Bs(multisite_tensor, spin_dims_right, positions_right, chi_limit, svd_settings);
        // We expect stashed U and S. Merge these and send onward to A's
        auto &U_stash = mps_sites_Bs.front().get_U_stash();
        auto &S_stash = mps_sites_Bs.front().get_S_stash();
        if(U_stash and S_stash) {
            U = U_stash->data.contract(Textra::asDiagonal(S_stash->data), Textra::idx({2}, {0}));
            U_stash.reset();
            S_stash.reset();
            auto t_split_svda = tools::common::profile::prof[AlgorithmType::ANY]["t_split_svda"]->tic_token();
            mps_sites_As      = internal::split_mps_into_As(U, spin_dims_left, positions_left, chi_limit, svd_settings);
            // Remnant S is now an LC for the last A
            auto &LC_stash = mps_sites_As.back().get_S_stash();
            auto &mpsA     = mps_sites_As.back();
            mpsA.set_LC(LC_stash->data, LC_stash->error);
            LC_stash.reset();
            mps_sites_Bs.front().merge_stash(mpsA);
        }
    }

    if(mps_sites_As.empty() and mps_sites_Bs.empty()) throw std::runtime_error("Got empty mps_sites from both left and right");

    // Some sanity checks
    if(mps_sites_As.size() != spin_dims_left.size())
        throw std::runtime_error(
            fmt::format("Could not split multisite tensor: Got mps_sites_As.size() {} != spins_dims_left {}", mps_sites_As.size(), spin_dims_left.size()));
    if(mps_sites_Bs.size() != spin_dims_right.size())
        throw std::runtime_error(
            fmt::format("Could not split multisite tensor: Got mps_sites_Bs.size() {} != spins_dims_right {}", mps_sites_Bs.size(), spin_dims_right.size()));

    // Move the right side onto the left side (This is equivalent to std::list::splice)
    mps_sites_As.insert(mps_sites_As.end(), std::make_move_iterator(mps_sites_Bs.begin()), std::make_move_iterator(mps_sites_Bs.end()));

    for(const auto &[idx, mps] : iter::enumerate(mps_sites_As)) {
        auto pos = positions[idx];
        if(pos != mps.get_position())
            throw std::runtime_error(fmt::format("Could not split multisite tensor: Position mismatch: expected site {} != mps pos {} | positions to merge {}",
                                                 pos, mps.get_position(), positions));
        if(mps.get_position<long>() == center_position and not mps.isCenter())
            throw std::runtime_error(
                fmt::format("Could not split multisite tensor: Specified center position {} did not become a center MPS", mps.get_position<long>()));
    }

    // Return the all the positions including the bond matrix
    return mps_sites_As;
}

std::vector<class_mps_site> tools::common::split::internal::split_mps_into_As(const Eigen::Tensor<Scalar, 3> &multisite_mps, const std::vector<long> &spin_dims,
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
    if(positions.empty()) return std::vector<class_mps_site>();

    // Initialize the resulting container of split sites
    std::vector<class_mps_site> mps_sites;

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
         */

        if(S_prev) {
            // Let V absorb S from the previous SVD (see note below)
            SV_temp.resize(V.dimensions());
            SV_temp.device(Textra::omp::getDevice()) = Textra::asDiagonal(S_prev.value()).contract(V, Textra::idx({1}, {1})).shuffle(Textra::array3{1, 0, 2});
            V                                        = SV_temp;
        }
        auto t_splitA_svd = tools::common::profile::prof[AlgorithmType::ANY]["t_splitA_svd"]->tic_token();
        auto t_svd        = tools::common::profile::get_default_prof()["t_svd"]->tic_token();
        std::tie(U, S, V) = svd.schmidt_into_left_normalized(V, spin_dim, chi_limit);
        t_svd.toc();
        t_splitA_svd.toc();
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

    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_wrk"] += *svd.t_wrk;
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_adj"] += *svd.t_adj;
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_jac"] += *svd.t_jac;
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_svd"] += *svd.t_svd;

    // Now we have a series of A-A-A-A matrices and their corresponding L's
    // At the last step we have residual S and V left over. Stash them!
    auto &mps = mps_sites.back();
    auto  pos = mps.get_position();
    if(S_prev) mps.stash_S(S_prev.value(), S_prev_error, pos + 1);
    mps.stash_V(V, pos + 1);
    return mps_sites;
}

std::deque<class_mps_site>
    // std::pair<Eigen::Tensor<class_mps_site::Scalar,3>,std::vector<class_mps_site>>
    tools::common::split::internal::split_mps_into_Bs(const Eigen::Tensor<class_mps_site::Scalar, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                      const std::vector<size_t> &positions, long chi_limit, std::optional<svd::settings> svd_settings) {
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
    if(positions.empty()) return std::deque<class_mps_site>();

    // Initialize the resulting container of split sites
    std::deque<class_mps_site> mps_sites;
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
            // Let V absorb S from the previous SVD (see note below)
            US_temp.resize(U.dimensions());
            US_temp.device(Textra::omp::getDevice()) = U.contract(Textra::asDiagonal(S), Textra::idx({2}, {0}));
            U                                        = US_temp;
        }
        auto t_splitB_svd = tools::common::profile::prof[AlgorithmType::ANY]["t_splitB_svd"]->tic_token();
        auto t_svd        = tools::common::profile::get_default_prof()["t_svd"]->tic_token();
        std::tie(U, S, V) = svd.schmidt_into_right_normalized(U, spin_dim, chi_limit);
        t_svd.toc();
        t_splitB_svd.toc();
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

    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_wrk"] += *svd.t_wrk;
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_adj"] += *svd.t_adj;
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_jac"] += *svd.t_jac;
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_svd"] += *svd.t_svd;

    // Now we have a series of B-B-B-B matrices and their corresponding L's
    // At the last step we have residual U and S left over. Stash them!
    auto &mps = mps_sites.front();
    auto  pos = mps.get_position();
    if(pos > 0) { // More left-normalized sites left in U
        if(S_prev) mps.stash_S(S_prev.value(), S_prev_error, pos - 1);
        mps.stash_U(U, pos - 1);
    }
    return mps_sites;
}
