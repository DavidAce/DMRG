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
                                                            std::optional<double> svd_threshold) {
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
     * S will later become LC after absorbing the remnants of (then split up) U and V.
     *
     *
     *
     * Special cases:
     *
     * 1) center_position is not in positions, and is to the left of all positions:
     *      Then we should only make a split from the right, turning all positions into "B" type positions.
     *      The U_residue matrix is
     *
     * 2) center_position is not in positions, and is to the right of all positions:
     *      Then we should only make a split from the left, turning all sites into "A" type sites.
     *      The V_residue matrix
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

    // Split the multisite tensor at the given center position.
    std::vector<size_t> positions_left, positions_right;
    std::vector<long>   spin_dims_left, spin_dims_right;

    // Define dimensions for the first reshape into a rank-4 tensor
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

    // Set up the SVD
    svd::solver svd;
    svd.use_lapacke = true;
    svd.use_bdc = true;
    svd.setLogLevel(2);
//    svd.enableProfiling();
    svd.setThreshold(settings::precision::svd_threshold, svd_threshold);
    svd.setSwitchSize(settings::precision::svd_switchsize);

    // Perform the first SVD which splits at the center position
    tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdm"]->tic();
    tools::common::profile::get_default_prof()["t_svd"]->tic();
    auto [U, S, V] = svd.schmidt_multisite(multisite_tensor, dL, dR, chiL, chiR, chi_limit);
    tools::common::profile::get_default_prof()["t_svd"]->toc();
    tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdm"]->toc();
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_wrk"] += *svd.t_wrk;
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_adj"] += *svd.t_adj;
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_jac"] += *svd.t_jac;
    *tools::common::profile::prof[AlgorithmType::ANY]["t_svd_svd"] += *svd.t_svd;

    // Sanity checks
    if(S.size() == 0) throw std::runtime_error("Could not split multisite tensor: Got 0 singular values from main svd");
    if(multisite_tensor.dimension(1) != U.dimension(1))
        throw std::runtime_error(fmt::format("Could not split multisite tensor: left bond dimension mismatch: before split {} != {} after split",
                                             multisite_tensor.dimension(1), U.dimension(1)));
    if(multisite_tensor.dimension(2) != V.dimension(2))
        throw std::runtime_error(fmt::format("Could not split multisite tensor: right bond dimension mismatch: before split {} != {} after split",
                                             multisite_tensor.dimension(2), V.dimension(2)));

    // Perform more SVD's to split the left and right parts around the center
    tools::common::profile::prof[AlgorithmType::ANY]["t_split_svda"]->tic();
    auto mps_sites_As = internal::split_mps_into_As(U, spin_dims_left, positions_left, chi_limit, svd_threshold);
    tools::common::profile::prof[AlgorithmType::ANY]["t_split_svda"]->toc();
    tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdb"]->tic();
    auto mps_sites_Bs = internal::split_mps_into_Bs(V, spin_dims_right, positions_right, chi_limit, svd_threshold);
    tools::common::profile::prof[AlgorithmType::ANY]["t_split_svdb"]->toc();

    if(mps_sites_As.empty() and mps_sites_Bs.empty()) throw std::runtime_error("Got empty mps_sites from both left and right");

    if(not mps_sites_As.empty() and mps_sites_Bs.empty()) {
        // We got only left-normalized "A" matrices from the split.
        if(multisite_tensor.dimension(1) != mps_sites_As.front().get_chiL())
            throw std::runtime_error(fmt::format("Could not split multisite tensor: left bond dimension mismatch: "
                                                 "multisite_tensor {} | mps_sites_As.front() {} | U {} | S {} | V {}",
                                                 multisite_tensor.dimensions(), mps_sites_As.front().dimensions(), U.dimensions(), S.dimensions(),
                                                 V.dimensions()));
        mps_sites_As.back().stash_V(V);
        if(mps_sites_As.back().get_position<long>() == center_position)
            mps_sites_As.back().set_LC(S, svd.get_truncation_error());
        else
            mps_sites_As.back().stash_S(S, svd.get_truncation_error());

    } else if(mps_sites_As.empty() and not mps_sites_Bs.empty()) {
        // We got only right-normalized "B" matrices from the split.
        if(multisite_tensor.dimension(2) != mps_sites_Bs.back().get_chiR())
            throw std::runtime_error(fmt::format("Could not split multisite tensor: right bond dimension mismatch: "
                                                 "multisite_tensor {} | mps_sites_Bs.back() {} | U {} | S {} | V {}",
                                                 multisite_tensor.dimensions(), mps_sites_Bs.back().dimensions(), U.dimensions(), S.dimensions(),
                                                 V.dimensions()));
        mps_sites_Bs.front().stash_S(S, svd.get_truncation_error());
        mps_sites_Bs.front().stash_U(U);
    } else {
        // We find this situation when doing SVD on a multisite tensor, i.e. >= 2 positions straddling the center point
        mps_sites_As.back().set_LC(S, svd.get_truncation_error());
    }

    // Some sanity checks
    if(mps_sites_As.size() != spin_dims_left.size()) throw std::runtime_error("Could not split multisite tensor: Got mps_sites_As of of unexpected size");
    if(mps_sites_Bs.size() != spin_dims_right.size()) throw std::runtime_error("Could not split multisite tensor: Got mps_sites_Bs of of unexpected size");

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
                                                                              std::optional<double> svd_threshold) {
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
    //    if(sites.empty()) throw std::runtime_error("Could not split multisite tensor from the left: sites list is empty");

    // Another  special case is when we do two-site tensors. Then we expect
    // this function to receive a single site. We can simply return it
    if(spin_dims.size() == 1) { return {class_mps_site(multisite_mps, std::nullopt, positions.front(), 0, "A")}; }
    // Initialize the resulting container of split sites
    std::vector<class_mps_site> mps_sites;

    // Now we have multiple spin dimensions

    // Set up the SVD
    svd::solver svd;
    svd.use_lapacke = true;
    svd.setThreshold(settings::precision::svd_threshold, svd_threshold);
    svd.setSwitchSize(settings::precision::svd_switchsize);

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
        tools::common::profile::prof[AlgorithmType::ANY]["t_splitA_svd"]->tic();
        tools::common::profile::get_default_prof()["t_svd"]->tic();
        std::tie(U, S, V) = svd.schmidt_from_left(V, spin_dim, chi_limit);
        tools::common::profile::get_default_prof()["t_svd"]->toc();
        tools::common::profile::prof[AlgorithmType::ANY]["t_splitA_svd"]->toc();
        if(S.size() == 0) throw std::runtime_error("Could not split multisite tensor: Got 0 singular values from left svd");

        // Note: Now we have the three components
        // U:   A left-unitary "A" matrix which is a "Lambda * Gamma" in Vidal's notation
        // S:   A set of singular values, "L" matrix, belonging to the site on the right of this one (i.e. the next A to pop out of V)
        // V:   Contains one site less than the it did before. That site is now in U.
        //      Before using V as the next multisite mps, it absorbs S so that the next U to pop out is a left-unitary A-type matrix.

        mps_sites.emplace_back(U, S_prev, positions[idx], S_prev_error, "A"); // S_prev is empty in the first iteration.

        // Store the singular values for the next iteration
        S_prev       = S;
        S_prev_error = svd.get_truncation_error();
    }

    // The spin dimension of V should be 1 since it does not contain sites anymore!
    if(V.dimension(0) != 1) throw std::runtime_error("Could not split mps from the left: Last V dimension(0) is not 1");

    // Now we have a series of A-A-A-A matrices and their corresponding L's
    // At the last step we have residual S and V left over.
    // U is already left-normalized but has scrambled rows. Use V to sort them!
    // S has all entries 1/dim(S) and can be discarded
    // V is a site-less tensor with scrambled unit-vectors as columns.

    // This operation should not change dimensions!
    Eigen::Tensor<Scalar, 3> UV =
        U.contract(V, Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{U.dimension(0), U.dimension(1), V.dimension(2)});
    // Note that after this M M^dagger = Identity
    mps_sites.back().set_M(UV);
    return mps_sites;
}

std::deque<class_mps_site>
    // std::pair<Eigen::Tensor<class_mps_site::Scalar,3>,std::vector<class_mps_site>>
    tools::common::split::internal::split_mps_into_Bs(const Eigen::Tensor<class_mps_site::Scalar, 3> &multisite_mps, const std::vector<long> &spin_dims,
                                                      const std::vector<size_t> &positions, long chi_limit, std::optional<double> svd_threshold) {
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

    // Another special case is when we do two-site tensors. Then we expect
    // this function to receive a single site. We can simply return it
    if(spin_dims.size() == 1) return {class_mps_site(multisite_mps, std::nullopt, positions.front(), 0, "B")};

    // Initialize the resulting container of split sites
    std::deque<class_mps_site> mps_sites;
    // Now we have multiple spin dimensions

    // Set up the SVD
    svd::solver svd;
    svd.use_lapacke = true;
    svd.setThreshold(settings::precision::svd_threshold, svd_threshold);
    svd.setSwitchSize(settings::precision::svd_switchsize);

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
        tools::common::profile::prof[AlgorithmType::ANY]["t_splitB_svd"]->tic();
        tools::common::profile::get_default_prof()["t_svd"]->tic();
        std::tie(U, S, V) = svd.schmidt_from_right(U, spin_dim, chi_limit);
        tools::common::profile::get_default_prof()["t_svd"]->toc();
        tools::common::profile::prof[AlgorithmType::ANY]["t_splitB_svd"]->toc();
        if(S.size() == 0) throw std::runtime_error("Could not split multisite tensor: Got 0 singular values from right svd");

        // Note: Now we have the three components
        // U:   Contains one site less than the previous iteration.
        //      Before using U as the next multisite mps, it absorbs the S so that the next V to pop out is a right-unitary B-type matrix.
        // S:   A set of singular values, "L" matrix, belonging to the left site  (i.e. the next B to pop out of U)
        // V:   A right-unitary "B" matrix which is a "Gamma*Lambda" in Vidal's notation

        mps_sites.emplace_front(V, S_prev, positions[idx], S_prev_error, "B"); // S_prev is empty the first iteration.

        // Store the singular values for the next iteration
        S_prev       = S;
        S_prev_error = svd.get_truncation_error();
    }
    // The spin dimension of U should be 1, since it contains no sites anymore!
    if(U.dimension(0) != 1) throw std::runtime_error("Could not split mps from the right: Last U dimension(0) is not 1");

    // Now we have a series of B-B-B-B matrices and their corresponding L's
    // At the last step we have residual U and S left over.
    // U is a site-less tensor with scrambled unit-vectors as rows.
    // S has all entries 1/dim(S) and can be discarded
    // V is already right-normalized but has scrambled columns. Use U to sort them!
    // This operation should not change dimensions!
    Eigen::Tensor<Scalar, 3> UV =
        U.contract(V, Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(Textra::array3{V.dimension(0), U.dimension(1), V.dimension(2)});
    // Note that after this M M^dagger = Identity
    mps_sites.front().set_M(UV);
    return mps_sites;
}
