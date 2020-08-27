#include "split.h"
#include "prof.h"
#include <config/nmspc_settings.h>
#include <math/svd.h>
#include <optional>
#include <tensors/state/class_mps_site.h>
#include <tools/common/fmt.h>

std::list<class_mps_site> tools::common::split::split_mps(const Eigen::Tensor<Scalar, 3> &multisite_tensor, const std::vector<long> &spin_dims,
                                                          const std::vector<size_t> &sites, size_t center_position, long chi_limit,
                                                          std::optional<double> svd_threshold) {
    /*  Here we split an mps containing multiple sites into its consituent sites.
     *  Consider a case with 3 sites and
     *  spin_dims = {2,2,2}
     *  sites = {6,7,8}
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
     * sites = {3,4,5}
     * center_position = 4
     *
     * then we would get
     * chiL---[A3]---chi   chi---[L]---chi   chi---[A4]---chi  chi---[LC]---chi   chi---[B5]---chiR
     *         |                                    |                                    |
     *        d=2                                  d=2                                  d=2
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
     * 1) center_position is not in sites, and is smaller than all sites (i.e. to the left)
     *      Then we should only make a split from the right, turning all sites into "B" type sites.
     *      The U_residue matrix is
     *
     * 2) center_position is not in sites, and is smaller than all sites (i.e. to the left)
     *      Then we should only make a split from the left, turning all sites into "A" type sites.
     *      The V_residue matrix
     *
     *
     */
    if(sites.empty()) throw std::runtime_error("Could not split multisite tensor: sites list is empty");
    if(multisite_tensor.size() == 0) throw std::runtime_error("Could not split multisite tensor: tensor is empty");
    if(spin_dims.empty()) throw std::runtime_error("Could not split multisite tensor: spin_dims list is empty");
    if(center_position < sites.front() or center_position > sites.back())
        throw std::runtime_error(
            fmt::format("Could not split multisite tensor: center position [{}] is not in the given list of sites {}", center_position, sites));
    if(spin_dims.size() != sites.size())
        throw std::runtime_error(
            fmt::format("Could not split multisite tensor: size mismatch in given lists: spin_dims {} != sites {} -- sizes not equal", spin_dims, sites));
    if(chi_limit <= 0) throw std::runtime_error(fmt::format("Invalid bond dimension limit:  chi_limit = {}", chi_limit));

    // Split the multisite tensor at the given center position.
    std::vector<size_t> sites_left;
    std::vector<size_t> sites_right;
    std::vector<long>   spin_dims_left;
    std::vector<long>   spin_dims_right;

    // Define dimensions for the first reshape into a rank-4 tensor
    long chiL = multisite_tensor.dimension(1);
    long chiR = multisite_tensor.dimension(2);
    long dL   = 1;
    long dR   = 1;

    auto site_it = sites.begin();
    auto dims_it = spin_dims.begin();
    while(site_it != sites.end() and dims_it != spin_dims.end()) {
        if(*site_it <= center_position) {
            sites_left.emplace_back(*site_it);
            spin_dims_left.emplace_back(*dims_it);
            dL *= *dims_it;
        } else {
            sites_right.emplace_back(*site_it);
            spin_dims_right.emplace_back(*dims_it);
            dR *= *dims_it;
        }
        site_it++;
        dims_it++;
    }

    // Set up the SVD
    if(not svd_threshold) svd_threshold = settings::precision::svd_threshold;
    svd::solver svd;
    svd.setThreshold(svd_threshold.value());
    tools::common::profile::t_svd->tic();
    auto [U, S, V] = svd.schmidt_multisite(multisite_tensor, dL, dR, chiL, chiR, chi_limit);
    tools::common::profile::t_svd->toc();
    if(S.size() == 0) throw std::runtime_error("Could not split multisite tensor: Got 0 singular values from main svd");
    auto mps_sites_left = internal::split_mps_from_left(U, spin_dims_left, sites_left, chi_limit, svd_threshold);
    //    tools::log->debug("SV dims: {} {} {}", SV.dimension(0),SV.dimension(1), SV.dimension(2));
    auto mps_sites_right = internal::split_mps_from_right(V, spin_dims_right, sites_right, chi_limit, svd_threshold);
    //    tools::log->debug("US dims: {} {} {}", US.dimension(0),US.dimension(1), US.dimension(2));

    //    Eigen::Tensor<Scalar,2> C = SV.reshape(Textra::array2{SV.dimension(0)*SV.dimension(1),SV.dimension(2)})
    //                                 .contract(Textra::asDiagonal(S), Textra::idx({1},{0}))
    //                                 .contract(US.reshape(Textra::array2{US.dimension(0)*US.dimension(1),US.dimension(2)}), Textra::idx({1},{0}));
    mps_sites_left.back().set_LC(S, svd.get_truncation_error());
    if(sites.size() > 2) {
        Eigen::Tensor<Scalar, 4> theta = mps_sites_left.back()
                                             .get_M_bare()
                                             .contract(Textra::asDiagonal(S), Textra::idx({2}, {0}))
                                             .contract(mps_sites_right.front().get_M_bare(), Textra::idx({2}, {1}));
        auto [u, s, v] = svd.schmidt(theta, chi_limit);
        mps_sites_left.back().set_M(u);
        mps_sites_left.back().set_LC(s, svd.get_truncation_error());
        mps_sites_right.front().set_M(v);
    }

    if constexpr(settings::debug) {
        //    Eigen::Tensor<Scalar,2> C = SV.reshape(Textra::array2{SV.dimension(0)*SV.dimension(1),SV.dimension(2)})
        //                                 .contract(Textra::asDiagonal(S), Textra::idx({1},{0}))
        //                                 .contract(US.reshape(Textra::array2{US.dimension(0)*US.dimension(1),US.dimension(2)}), Textra::idx({1},{0}));
        //        Eigen::Tensor<Scalar,4> theta = mps_sites_left.back().get_M_bare()
        //            .contract(Textra::asDiagonal(S), Textra::idx({2},{0}))
        //            .contract(mps_sites_right.front().get_M_bare(), Textra::idx({2},{1}));
        //    std::cout << "S = \n" << S << std::endl;
        //        std::tie(U,S,V) = svd.schmidt(theta,chi_limit);
        //        mps_sites_left.back().set_M(U);
        //    mps_sites_left.back().set_LC(S,svd.get_truncation_error());
        //        mps_sites_right.front().set_M(V);
        //    std::cout << "C = \n" << C << std::endl;
        //    std::cout << "C = \n" << Textra::extractDiagonal(C) << std::endl;
        //    std::cout << "S = \n" << S << std::endl;
        //    S = Textra::extractDiagonal(C);
    }

    // Some sanity checks
    if(mps_sites_left.size() != spin_dims_left.size()) throw std::runtime_error("Could not split multisite tensor: Got mps_sites_left of of unexpected size");
    if(mps_sites_right.size() != spin_dims_right.size())
        throw std::runtime_error("Could not split multisite tensor: Got mps_sites_right of of unexpected size");

    // Append the center bond to the left side
    //    mps_sites_left.back().set_LC(S, svd.get_truncation_error());
    // Move the right side onto the left side.
    mps_sites_left.splice(mps_sites_left.end(), mps_sites_right);

    site_it = sites.begin();
    for(const auto &mps : mps_sites_left) {
        if(*site_it != mps.get_position())
            throw std::runtime_error(fmt::format("Could not split multisite tensor: Position mismatch: expected site {} != mps pos {} | sites to merge {}",
                                                 *site_it, mps.get_position(), sites));
        if(mps.get_position() == center_position and not mps.isCenter())
            throw std::runtime_error(
                fmt::format("Could not split multisite tensor: Specified center position {} did not become a center MPS", mps.get_position()));
        site_it++;
    }

    // Return the all the sites including the bond matrix
    return mps_sites_left;
}

std::list<class_mps_site>
    // std::pair<std::list<class_mps_site>,Eigen::Tensor<class_mps_site::Scalar,3>>
    tools::common::split::internal::split_mps_from_left(const Eigen::Tensor<Scalar, 3> &multisite_mps, std::vector<long> spin_dims, std::vector<size_t> sites,
                                                        long chi_limit, std::optional<double> svd_threshold) {
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

    if(sites.empty()) throw std::runtime_error("Could not split multisite tensor from the left: sites list is empty");
    if(spin_dims.empty()) throw std::runtime_error("Could not split multisite tensor from the left: spin_dims list is empty");
    if(spin_dims.size() != sites.size())
        throw std::runtime_error(fmt::format(
            "Could not split multisite tensor from the left: size mismatch in given lists: spin_dims {} != sites {} -- sizes not equal", spin_dims, sites));

    // A special case is when we do two-site tensors. Then we expect
    // this function to receive a single site. We can simply return it
    if(spin_dims.size() == 1) {
        return {class_mps_site(multisite_mps, std::nullopt, sites.front())};
        //        mps_sites.emplace_back(multisite_mps,std::nullopt,sites.front());
        //        return mps_sites;
        //        auto diag = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>::Identity(multisite_mps.dimension(2),multisite_mps.dimension(2));
        //        Eigen::Tensor<Scalar,3> diag_tensor = Textra::MatrixTensorMap(diag,1,multisite_mps.dimension(2),multisite_mps.dimension(2));
        //        return {mps_sites,diag_tensor};
    }
    // Initialize the resulting container of split sites
    std::list<class_mps_site> mps_sites;

    // Now we have multiple spin dimensions

    // Set up the SVD
    if(not svd_threshold) svd_threshold = settings::precision::svd_threshold;
    svd::solver svd;
    svd.setThreshold(svd_threshold.value());

    // Declare the tensors that will catch the schmidt (SVD) decompositions
    Eigen::Tensor<Scalar, 3>                U;                      // This will become the first site extracted
    Eigen::Tensor<Scalar, 1>                S;                      // The singular values
    Eigen::Tensor<Scalar, 3>                V      = multisite_mps; // This side contains all the sites
    std::optional<Eigen::Tensor<Scalar, 1>> S_prev = std::nullopt;  // Starts out empty, needs to be checked outside of this split
    for(const auto &spin_dim : spin_dims) {
        /* The schmidt decomposition gives us 3 matrices:
         *
         * chiL ---[U]----chi chi---[S]---chi chi---[V]---chiR
         *          |                                |
         *          d                             d^(sites-1)
         *
         * Here U is an "M" matrix of type A = Lambda * Gamma
         * See the notes on svd.schmidt at its definition
         */
        if(sites.empty()) throw std::logic_error("Could not split multisite tensor from the right: Site list became empty");

        tools::common::profile::t_svd->tic();
        std::tie(U, S, V) = svd.schmidt_from_left(V, spin_dim, chi_limit);
        tools::common::profile::t_svd->toc();
        if(S.size() == 0) throw std::runtime_error("Could not split multisite tensor: Got 0 singular values from left svd");

        // Now we have the three components
        // U:   A left-unitary "A" matrix which already contains its left singular values
        // S:   A set of singular values, "L" matrix, belonging to the right site
        //      (i.e. the next A to pop out of V)
        // V:   Contains one site less than the previous iteration.
        //      Absorbs the S so that the next U to pop out is a left-unitary A-type matrix.
        mps_sites.emplace_back(U, S_prev, sites.front(), svd.get_truncation_error()); // S_prev is empty the first iteration.

        // Let V absorb S
        Eigen::Tensor<Scalar, 3> temp = Textra::asDiagonal(S).contract(V, Textra::idx({1}, {1})).shuffle(Textra::array3{1, 0, 2});
        V                             = temp;

        // Store the singular values for the next iteration
        S_prev = S;
        sites.erase(sites.begin());
    }
    // Now we have a series of A-A-A-A matrices and their corresponding L's
    // At the last step we have a residual S*V left over, where V has already absorbed S.
    // The spin dimension of S*V should be 1, and together S*V may contain important signs that can't be ignored.
    // We just tuck them back into the last  A matrix.
    // This operation should not change dimensions!
    if(V.dimension(0) != 1) throw std::runtime_error("Could not split mps from the left: Last V dimension(0) is not 1");
    Eigen::Tensor<Scalar, 3> USV = U.contract(V, Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(mps_sites.back().dimensions());
    mps_sites.back().set_M(USV);
    return mps_sites;
}

std::list<class_mps_site>
    // std::pair<Eigen::Tensor<class_mps_site::Scalar,3>,std::list<class_mps_site>>
    tools::common::split::internal::split_mps_from_right(const Eigen::Tensor<class_mps_site::Scalar, 3> &multisite_mps, std::vector<long> spin_dims,
                                                         std::vector<size_t> sites, long chi_limit, std::optional<double> svd_threshold) {
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

    if(sites.empty()) throw std::runtime_error("Could not split multisite tensor from the right: sites list is empty");
    if(spin_dims.empty()) throw std::runtime_error("Could not split multisite tensor from the right: spin_dims list is empty");
    if(spin_dims.size() != sites.size())
        throw std::runtime_error(fmt::format(
            "Could not split multisite tensor from the right: size mismatch in given lists: spin_dims {} != sites {} -- sizes not equal", spin_dims, sites));

    // A special case is when we do two-site tensors. Then we expect
    // this function to receive a single site. We can simply return it
    if(spin_dims.size() == 1) { return {class_mps_site(multisite_mps, std::nullopt, sites.front())}; }

    // Initialize the resulting container of split sites
    std::list<class_mps_site> mps_sites;
    // Now we have multiple spin dimensions
    // Reverse the order of spin_dims so we can iterate starting from its right side.
    std::reverse(spin_dims.begin(), spin_dims.end());

    // Set up the SVD
    if(not svd_threshold) svd_threshold = settings::precision::svd_threshold;
    svd::solver svd;
    svd.setThreshold(svd_threshold.value());

    // Declare the the tensors that will catch the schmidt (SVD) decompositions
    Eigen::Tensor<Scalar, 3>                U = multisite_mps;     // This side contains all the sites
    Eigen::Tensor<Scalar, 1>                S;                     // The singular values
    Eigen::Tensor<Scalar, 3>                V;                     // This will become the first site extracted
    std::optional<Eigen::Tensor<Scalar, 1>> S_prev = std::nullopt; // Starts out empty, needs to be checked outside of this split
    for(const auto &spin_dim : spin_dims) {
        /* The schmidt decomposition gives us 3 matrices:
         *
         * chiL ---[U]----chi chi---[S]---chi chi---[V]---chiR
         *          |                                |
         *      d^(sites-1)                          d
         *
         * Here V is an "M" matrix of type B = Gamma * Lambda
         * See the notes on svd.schmidt at its definition
         */

        if(sites.empty()) throw std::logic_error("Could not split multisite tensor from the right: No sites left");

        tools::common::profile::t_svd->tic();
        std::tie(U, S, V) = svd.schmidt_from_right(U, spin_dim, chi_limit);
        tools::common::profile::t_svd->toc();
        if(S.size() == 0) throw std::runtime_error("Could not split multisite tensor: Got 0 singular values from right svd");

        // Now we have the three components
        // U:   Contains one site less than the previous iteration.
        //      Absorbs the S so that the next V to pop out is a right-unitary B-type matrix.
        // S:   A set of singular values, "L" matrix, belonging to the left site
        //      (i.e. the next B to pop out of U)
        // V:   A right-unitary "B" matrix which already contains its right singular values

        mps_sites.emplace_front(V, S_prev, sites.back(), svd.get_truncation_error()); // S_prev is empty the first iteration.

        // Let U absorb S
        Eigen::Tensor<Scalar, 3> temp = U.contract(Textra::asDiagonal(S), Textra::idx({2}, {0}));
        U                             = temp;

        // Store the singular values for the next iteration
        S_prev = S;
        sites.pop_back();
    }
    // Now we have a series of B-B-B-B matrices and their corresponding L's
    // At the last step we have a residual U*S, where U has already absorbed S.
    // The spin dimension of U should be 1, and together U*S may contain important data that can't be ignored.
    // We just tuck them back into the front B matrix, and perform one more SVD on A S B later.
    // This operation should not change dimensions!
    if(U.dimension(0) != 1) throw std::runtime_error("Could not split mps from the right: Last U dimension(0) is not 1");
    Eigen::Tensor<Scalar, 3> USV = U.contract(V, Textra::idx({2}, {1})).shuffle(Textra::array4{0, 2, 1, 3}).reshape(V.dimensions());

    mps_sites.front().set_M(USV);
    return mps_sites;
}
