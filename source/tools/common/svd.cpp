#include "svd.h"
#include "prof.h"
#include <optional>
#include <config/nmspc_settings.h>
#include <math/class_svd_wrapper.h>
#include <math/nmspc_math.h>
#include <tensors/state/class_mps_site.h>

//extern std::list<class_mps_site> split_mps (const Eigen::Tensor<Scalar,3> & multisite_mps,
//                                            long                            chi_limit,
//                                            std::optional<double>           svd_threshold = std::nullopt);
//
//extern std::list<class_mps_site> split_mps (const Eigen::Tensor<Scalar,3> & multisite_mps,
//                                            const std::list<long>         & spin_dims,
//                                            const std::list<size_t>       & positions,
//                                            size_t                          center_position,
//                                            long                            chi_limit,
//                                            std::optional<double>           svd_threshold = std::nullopt);


std::list<class_mps_site> tools::common::svd::split_mps (const Eigen::Tensor<Scalar,3> & multisite_mps,
                                                        long                             chi_limit,
                                                        std::optional<std::list<long>>   spin_dims,
                                                        std::optional<double>            svd_threshold){
    /* This function is used primarily to split 2-site infinite mps */
    // Start by figuring out the spin dimensions d
    if(not spin_dims) {
        long multisite_spin_dim = multisite_mps.dimension(0);
        if(math::mod(multisite_spin_dim, 2) != 0) throw std::runtime_error("Could not split svd: spin dimensions not given and multisite mps spin dim is not divisible by 2 (the default)");
        auto sites = static_cast<size_t>(std::log2(multisite_spin_dim));
        spin_dims = std::list<long>(sites, 2);
    }
    std::list<size_t> positions;
    for(size_t pos = 0; pos < spin_dims->size(); pos++)
        positions.emplace_back(pos);

    return tools::common::svd::split_mps(multisite_mps,spin_dims.value(),positions,0,chi_limit,svd_threshold);
}


std::list<class_mps_site> tools::common::svd::split_mps(const Eigen::Tensor<Scalar,3> & multisite_mps,
                                                        std::list<long>                 spin_dims,
                                                        std::list<size_t>               positions,
                                                        size_t                          center_position,
                                                        long                            chi_limit,
                                                        std::optional<double>           svd_threshold){

    /*  Here we split an mps containing multiple sites into its consituent sites.
     *  Consider a case with 3 sites and
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
     * 1) center_position is not in positions, and is smaller than all positions (i.e. to the left)
     *      Then we should only make a split from the right, turning all sites into "B" type sites.
     *      The U_residue matrix is
     *
     * 2) center_position is not in positions, and is smaller than all positions (i.e. to the left)
     *      Then we should only make a split from the left, turning all sites into "A" type sites.
     *      The V_residue matrix
     *
     *
     */

    if(center_position < positions.front() or center_position > positions.front())
        throw std::runtime_error("Center position {} is not in the given list of positions {}", center_position, positions);

    if(spin_dims.size() != positions.size())
        throw std::runtime_error("Size mismatch in given lists: spin dims "
                                 + std::to_string(spin_dims.size())
                                 + " != positions "
                                 + std::to_string(positions.size()));



    // Split the multisite mps at the given center position.
    std::list<size_t> positions_left;
    std::list<long>   spin_dims_left;
    std::list<size_t> positions_right;
    std::list<long>   spin_dims_right;

    // Define dimensions for the first reshape into a rank-4 tensor
    long chiL = multisite_mps.dimension(1);
    long chiR = multisite_mps.dimension(2);
    long dL  = 1;
    long dR  = 1;

    auto pos_it = positions.begin();
    auto dim_it = spin_dims.begin();
    while (pos_it != positions.end() and dim_it != spin_dims.end()){
        if(*pos_it <= center_position){
            positions_left.emplace_back(*pos_it);
            spin_dims_left.emplace_back(*dim_it);
            dL *= *dim_it;
        }else{
            positions_right.emplace_back(*pos_it);
            spin_dims_right.emplace_back(*dim_it);
            dR *= *dim_it;
        }
        pos_it++;
        dim_it++;
    }


    // Set up the SVD
    if(not svd_threshold) svd_threshold = settings::precision::svd_threshold;
    class_SVD svd;
    svd.setThreshold(svd_threshold.value());

    auto [U,S,V] = svd.schmidt(multisite_mps,dL,dR,chiL,chiR,chi_limit);

    auto [mps_sites_left, V_residue]  = internal::split_mps_from_left (U, spin_dims_left, positions_left,chi_limit,svd_threshold);
    auto [mps_sites_right, U_residue] = internal::split_mps_from_right(V, spin_dims_right, positions_right,chi_limit,svd_threshold);

    // Some sanity checks
    if(mps_sites_left.size() != spin_dims_left.size())
        throw std::runtime_error("Got mps_sites_left of of unexpected size");
    if(mps_sites_right.size() != spin_dims_right.size())
        throw std::runtime_error("Got mps_sites_right of of unexpected size");


    //////////////// DEPRECATED ////////////////////
    // We have 3 possible cases here depending on the center_position
    //  * center_point < all positions:
    //      * dL = 1, dR = multisite_mps.dimension(0) = prod(spin_dims)
    //      * mps_sites_left is empty
    //      * V_residue = U
    //      * "LC" = V_residue * S * U_residue = U * S * U_residue
    //      * Contract LC into the first B (not setting LC!)
    //      * Splice mps_sites_right onto mps_sites_left as usual
    //  * center_point > all positions:
    //      * dL = multisite_mps.dimension(0) = prod(spin_dims) , dR = 1
    //      * mps_sites_right is empty
    //      * U_residue = V
    //      * "LC" = V_residue * S * U_residue = V_residue * S * V
    //      * Contract LC into the last A (not setting LC!)
    //      * Splice mps_sites_right onto mps_sites_left is not needed, but can be done anyway to reuse code
    //  * center_point is somewhere in positions (normal case):
    //      * LC = V_residue * S * U_residue
    //      * Set LC into the last site in mps_sites_left
    //      * Splice mps_sites_right onto mps_sites_left as usual
    //
    //////////////// DEPRECATED ////////////////////



    // Take care of the center bond, which should be V_residue * S * U_residue.
    Eigen::Tensor<Scalar,2> center_matrix = U_residue
                                          .contract(Textra::asDiagonal(S), Textra::idx({2},{0}))
                                          .contract(V_residue            , Textra::idx({2},{1}))
                                          .reshape(Textra::array2{S.size(),S.size()});

    Eigen::Tensor<Scalar,1> LC = Textra::extractDiagonal(center_matrix);

    // Append the center bond to the left side
    mps_sites_left.back().set_LC(LC,svd.get_truncation_error());
    //Move the right side onto the left side.
    mps_sites_left.splice(mps_sites_left.end(), mps_sites_right);
    // Return the all the sites including the bond matrix
    return  mps_sites_left;

//    // Initialize the resulting container of split sites
//    std::list<class_mps_site> mps_sites;
//
//    // There is a special case when there is only one site
//    // which cannot be split further. Simply return it.
//    if(spin_dims.size() == 1){
//        if(multisite_mps.dimension(0) != spin_dims.front())
//            throw std::logic_error("Spin dimension mismatch in given mps and spin dim list");
//        mps_sites.emplace_front(class_mps_site());
//        mps_sites.front().set_position(0);
//        mps_sites.front().set_M(multisite_mps);
//        return mps_sites;
//    }
//
//
//    // Now we have multiple spin dimensions
//    // Reverse the order of spin_dims so we can iterate starting from its right side.
//    std::reverse(spin_dims.begin(), spin_dims.end());
//
//
//    // Declare the the tensors that will catch the schmidt (SVD) decompositions
//    Eigen::Tensor<Scalar, 3>  U = multisite_mps; // This side contains all the sites
//    Eigen::Tensor<Scalar, 1>  S;                 // The singular values
//    Eigen::Tensor<Scalar, 3>  V;                 // This will become the first site extracted
//    Eigen::Tensor<Scalar, 1>  S_prev;            // Starts out empty, needs to be checked outside of this split
//    auto                      position = positions.end();
//    for(const auto &spin_dim : spin_dims.value()) {
//        /* The schmidt decomposition gives us 3 matrices:
//         *
//         * chiL ---[U]----chi chi---[S]---chi chi---[V]---chiR
//         *          |                                |
//         *      d^(sites-1)                          d
//         *
//         * Here V is an "M" matrix of type B = Gamma * Lambda
//         * See the notes on svd.schmidt at its definition
//         */
//
//        tools::common::profile::t_svd->tic();
//        std::tie(U, S, V) = svd.schmidt(U, chi_limit);
//        tools::common::profile::t_svd->toc();
//
//        // Now we have the three components
//        // U:   Contains one site less than the previous iteration.
//        //      Absorbs the S so that the next V to pop out is a right-unitary B-type matrix.
//        // S:   A set of singular values, "L" matrix, belonging to the left site
//        //      (i.e. the next B to pop out of U)
//        // V:   A right-unitary "B" matrix which already contains its right singular values
//
//        mps_sites.emplace_front(V, S_prev, --position); // S_prev is empty the first iteration.
//
//        // Let U absorb S
//        Eigen::Tensor<Scalar, 3> temp = U.contract(Textra::asDiagonal(S), Textra::idx({2}, {0}));
//        U                             = temp;
//
//        // If the spin dimension of U is equal to spin_dim, then we have reached the last site!
//        if(U.dimension(0) != spin_dim){
//            // Store the singular values for the next iteration
//            S_prev = S;
//        }else{
//            // Now U is a left-unitary "A" matrix, and the singular value S is a center bond "LC"
//            // We can distinguish it by setting S in its LC
//            mps_sites.emplace_front(class_mps_site());
//            mps_sites.front().set_position(--position);
//            mps_sites.front().set_M(U);
//            mps_sites.front().set_LC(S);
//            break;
//        }
//    }
//    // Hopefully we have collected all sites now.
//    return mps_sites;
}



std::tuple<
    std::list<class_mps_site>,
    Eigen::Tensor<class_mps_site::Scalar,3>>
tools::common::svd::internal::split_mps_from_left(const Eigen::Tensor<Scalar,3> & multisite_mps,
                                                  std::list<long>                 spin_dims,
                                                  std::list<size_t>               positions,
                                                  long                            chi_limit,
                                                  std::optional<double>           svd_threshold){

    /*  Here we split an mps containing multiple sites into its consituent sites from the left.
     *  Consider a case with 3 sites and
     *  spin_dims = {2,2,2}
     *  positions = {6,7,8}
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
     * By convention, a multisite_mps is created by merging sites from left to right, that is
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

    // Initialize the resulting container of split sites
    std::list<class_mps_site> mps_sites;

    // There is a special case when there is only one site
    // which cannot be split further. Simply return it.
//    if(spin_dims.size() == 1){
//        if(multisite_mps.dimension(0) != spin_dims.front())
//            throw std::logic_error("Spin dimension mismatch in given mps and spin dim list");
//        mps_sites.emplace_front(class_mps_site());
//        mps_sites.front().set_position(0);
//        mps_sites.front().set_M(multisite_mps);
//        return mps_sites;
//    }


    // Now we have multiple spin dimensions

    // Set up the SVD
    if(not svd_threshold) svd_threshold = settings::precision::svd_threshold;
    class_SVD svd;
    svd.setThreshold(svd_threshold.value());

    // Declare the the tensors that will catch the schmidt (SVD) decompositions
    Eigen::Tensor<Scalar, 3>  U;  // This will become the first site extracted
    Eigen::Tensor<Scalar, 1>  S;                 // The singular values
    Eigen::Tensor<Scalar, 3>  V  = multisite_mps; // This side contains all the sites
    Eigen::Tensor<Scalar, 1>  S_prev;            // Starts out empty, needs to be checked outside of this split
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
        if(positions.empty())
            throw std::logic_error("Could not split multisite mps from the right: No positions left");

        tools::common::profile::t_svd->tic();
        std::tie(U, S, V) = svd.schmidt_from_left(V, spin_dim, chi_limit);
        tools::common::profile::t_svd->toc();

        // Now we have the three components
        // U:   A left-unitary "A" matrix which already contains its left singular values
        // S:   A set of singular values, "L" matrix, belonging to the right site
        //      (i.e. the next A to pop out of V)
        // V:   Contains one site less than the previous iteration.
        //      Absorbs the S so that the next U to pop out is a left-unitary A-type matrix.


        mps_sites.emplace_front(U, S_prev, positions.front(),svd.get_truncation_error()); // S_prev is empty the first iteration.

        // Let V absorb S
        Eigen::Tensor<Scalar,3> temp = Textra::asDiagonal(S).contract(V, Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
        V = temp;

        // Store the singular values for the next iteration
        S_prev = S;

        positions.pop_front();
    }
    // Now we have a series of A-A-A-A matrices and their corresponding L's
    // At the last step we have one V and one S left over, where S has already been absorbed into V.
    // These are important for calculating the center bond correctly, so we return them as well.
    if(U.dimension(0) != 1)
        throw std::runtime_error("Could not split mps from the right: Last U dimension is not 1");

    return std::make_tuple(mps_sites,V);
}



std::tuple<
    std::list<class_mps_site>,
    Eigen::Tensor<class_mps_site::Scalar,3>>
tools::common::svd::internal::split_mps_from_right(const Eigen::Tensor<class_mps_site::Scalar,3> & multisite_mps,
                                                         std::list<long>                           spin_dims,
                                                         std::list<size_t>                         positions,
                                                         long                                      chi_limit,
                                                         std::optional<double>                     svd_threshold){

    /*  Here we split an mps containing multiple sites into its consituent sites from the right
     *  Consider a case with 3 sites and
     *  spin_dims = {2,2,2}
     *  positions = {6,7,8}
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
     * By convention, a multisite_mps is created by merging sites from left to right, that is
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

    // Initialize the resulting container of split sites
    std::list<class_mps_site> mps_sites;

    // There is a special case when there is only one site
    // which cannot be split further. Simply return it.
//    if(spin_dims.size() == 1){
//        if(multisite_mps.dimension(0) != spin_dims.front())
//            throw std::logic_error("Spin dimension mismatch in given mps and spin dim list");
//        mps_sites.emplace_front(class_mps_site());
//        mps_sites.front().set_position(0);
//        mps_sites.front().set_M(multisite_mps);
//        return mps_sites;
//    }


    // Now we have multiple spin dimensions
    // Reverse the order of spin_dims so we can iterate starting from its right side.
    std::reverse(spin_dims.begin(), spin_dims.end());

    // Set up the SVD
    if(not svd_threshold) svd_threshold = settings::precision::svd_threshold;
    class_SVD svd;
    svd.setThreshold(svd_threshold.value());

    // Declare the the tensors that will catch the schmidt (SVD) decompositions
    Eigen::Tensor<Scalar, 3>  U = multisite_mps; // This side contains all the sites
    Eigen::Tensor<Scalar, 1>  S;                 // The singular values
    Eigen::Tensor<Scalar, 3>  V;                 // This will become the first site extracted
    Eigen::Tensor<Scalar, 1>  S_prev;            // Starts out empty, needs to be checked outside of this split
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

        if(positions.empty())
            throw std::logic_error("Could not split multisite mps from the right: No positions left");

        tools::common::profile::t_svd->tic();
        std::tie(U, S, V) = svd.schmidt_from_right(U,spin_dim, chi_limit);
        tools::common::profile::t_svd->toc();

        // Now we have the three components
        // U:   Contains one site less than the previous iteration.
        //      Absorbs the S so that the next V to pop out is a right-unitary B-type matrix.
        // S:   A set of singular values, "L" matrix, belonging to the left site
        //      (i.e. the next B to pop out of U)
        // V:   A right-unitary "B" matrix which already contains its right singular values

        mps_sites.emplace_front(V, S_prev, positions.back(),svd.get_truncation_error()); // S_prev is empty the first iteration.

        // Let U absorb S
        Eigen::Tensor<Scalar, 3> temp = U.contract(Textra::asDiagonal(S), Textra::idx({2}, {0}));
        U                             = temp;

        // Store the singular values for the next iteration
        S_prev = S;
        positions.pop_back();


    }
    // Now we have a series of B-B-B-B matrices and their corresponding L's
    // At the last step we have one U and one S left over, where S has already been absorbed into U.
    // These are important for calculating the center bond correctly, so we return them as well.
    if(U.dimension(0) != 1)
        throw std::runtime_error("Could not split mps from the right: Last U dimension is not 1");

    return std::make_tuple(mps_sites,U);

}
