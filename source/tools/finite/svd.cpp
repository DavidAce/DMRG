//
// Created by david on 2019-06-29.
//
#include <tools/finite/svd.h>
#include <tools/finite/opt.h>
#include <tools/finite/measure.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/common/svd.h>
#include <math/class_svd_wrapper.h>
#include <math/nmspc_math.h>
#include <config/nmspc_settings.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/class_tensors_finite.h>





//void tools::finite::svd::truncate_theta(const Eigen::Tensor<Scalar,3> &theta, class_state_finite & state, std::optional<size_t> chi_lim, std::optional<double> svd_threshold){
//    state.clear_measurements();
//    if(not chi_lim.has_value()) chi_lim = state.get_chi_lim();
//    if (state.active_sites.empty()) throw std::runtime_error("truncate_theta: No active sites to truncate");
//    if (theta.size() == 0)          throw std::runtime_error("truncate_theta: Theta is empty");
//    auto theta_map = Eigen::Map<const Eigen::VectorXcd>(theta.data(),theta.size());
//    auto fullnorm  = tools::finite::measure::norm(state);
//    auto thetanorm = theta_map.norm();
//    if(std::abs(fullnorm  - 1.0) > settings::precision::max_norm_error){ tools::log->error("Norm before truncation too far from unity: {:.16f}", fullnorm);}
//    if(std::abs(thetanorm - 1.0) > settings::precision::max_norm_error){ tools::log->error("Norm of given theta too far from unity: {:.16f}", thetanorm);}
//
//    //Start by clearing the environments that will become stale.
//    while(true){
//        if(state.active_sites.back() > state.ENV_R.front().get_position()){
//            state.ENV_R.pop_front();
//            state.ENV2_R.pop_front();
//        }
//        else if(state.active_sites.front() < state.ENV_L.back().get_position()){
//            state.ENV_L.pop_back();
//            state.ENV2_L.pop_back();
//        }
//        else break;
//    }
//    if (state.get_direction() == 1){
//        tools::finite::svd::truncate_left(theta,state, chi_lim,svd_threshold);
//
//    }else{
//        tools::finite::svd::truncate_right(theta,state,chi_lim,svd_threshold);
//    }
//
//    state.clear_measurements();
//    state.clear_cache();
//    if (state.ENV_L.size() + state.ENV_R.size() != state.get_length())
//        throw std::runtime_error(fmt::format("Environment sizes do not add up to system size after truncation: {} + {} != {}", state.ENV_L.size() , state.ENV_R.size(), state.get_length()));
//    if (state.ENV2_L.size() + state.ENV2_R.size() != state.get_length())
//        throw std::runtime_error(fmt::format("Environment sq sizes do not add up to system size after truncation: {} + {} != {}", state.ENV2_L.size() , state.ENV2_R.size(), state.get_length()));
//
//    fullnorm  = tools::finite::measure::norm(state);
//    if(std::abs(fullnorm  - 1.0) > 1e-10){
//        tools::log->warn("Norm after truncation too far from unity: {:.16f} -- Normalizing",fullnorm);
//        tools::finite::mps::normalize(state);
//        fullnorm  = tools::finite::measure::norm(state);
//        tools::log->warn("New norm: {:.16f}",fullnorm);
//    }
//}

//
//void tools::finite::svd::truncate_right(const Eigen::Tensor<std::complex<double>,3> &theta, class_state_finite & state, std::optional<size_t> chi_lim, std::optional<double> svd_threshold){
//    tools::log->trace("Truncating multitheta from left to right");
//    if(not svd_threshold.has_value()){
//        svd_threshold = settings::precision::svd_threshold;
//    }
//    class_SVD SVD;
//    SVD.setThreshold(svd_threshold.value());
//    using Scalar = class_state_finite::Scalar;
//    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
//    if(not chi_lim.has_value()) chi_lim = state.get_chi_lim();
//    Eigen::Tensor<Scalar,4> theta4;
//    Eigen::Tensor<Scalar,3> U;
//    Eigen::Tensor<Scalar,3> V = theta;
//    Eigen::Tensor<Scalar,1> S;
//    auto active_sites = state.active_sites;
//    double norm;
//    tools::common::profile::t_svd->tic();
//    while (active_sites.size() >= 2){
//        size_t site = active_sites.front();
//        long dim0 = state.get_mps(site).spin_dim();
//        long dim1 = V.dimension(0) / state.get_mps(site).spin_dim();
//        long dim2 = V.dimension(1);
//        long dim3 = V.dimension(2);
//        theta4 = V
//                .reshape(Textra::array4{dim0,dim1,dim2,dim3})
//                .shuffle(Textra::array4{0,2,1,3});
//        std::tie(U, S, V,norm) = SVD.schmidt_with_norm(theta4, chi_lim);
//
//        state.set_truncation_error(site, SVD.get_truncation_error());
//        state.get_mps(site).set_M(U);
//        state.get_mps(site).unset_LC();
//
//
//        if (not Eigen::Map<VectorType>(U.data(),U.size()).allFinite() )
//            throw std::runtime_error("L_U has nan's or inf's");
//
//        Eigen::Tensor<Scalar,2> leftID = state.get_mps(site).get_M()
//                .contract(state.get_mps(site).get_M().conjugate(), Textra::idx({0,1},{0,1}) );
//        if(not Textra::TensorMatrixMap(leftID).isIdentity(1e-12)) throw std::runtime_error(fmt::format("Not left normalized at site {} with threshold 1e-12.", site));
//
//
//
//        if(active_sites.size() >= 3){
//            Eigen::Tensor<Scalar,3> temp = Textra::asDiagonal(S).contract(V, Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
//            state.get_mps(site + 1).set_L(S);
//            V = temp;
//            if(state.ENV_L. back().get_position() != site) throw std::runtime_error(fmt::format("Site and postion mismatch in ENV_L while truncating  left to right {} != {}",state.ENV_L. front().get_position() , site));
//            if(state.ENV2_L.back().get_position() != site) throw std::runtime_error(fmt::format("Site and postion mismatch in ENV2_L while truncating  left to right{} != {}",state.ENV2_L.front().get_position() , site));
//
//            state.ENV_L .emplace_back(state.ENV_L .back().enlarge(state.get_mps(site), state.get_MPO(site)));
//            state.ENV2_L.emplace_back(state.ENV2_L.back().enlarge(state.get_mps(site), state.get_MPO(site)));
//        } else{
//            //Always set LC on the last "A" matrix
//            state.get_mps(site).set_LC(S);
//        }
//        tools::log->trace("SVD site {:2} log₁₀ trunc: {:12.8f} χlim: {:4} χ: {:4}", site, std::log10(state.get_truncation_error(site)),state.get_chi_lim(),
//                          state.get_mps(site).get_chiR());
//        active_sites.pop_front();
//    }
//
//
//    size_t site = active_sites.front();
////    Eigen::Tensor<Scalar,3> V_L = V.contract(Textra::asDiagonalInversed(state.get_L(site+1)), Textra::idx({2},{0}));
////    state.get_G(site) = V_L;
//    state.get_mps(site).set_M(V);
//
//    if (not Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1 >>(V.data(),V.size()).allFinite() )
//        throw std::runtime_error("V_L has nan's or inf's");
//
//
//    Eigen::Tensor<Scalar,2> rightID = state.get_mps(site).get_M()
//            .contract(state.get_mps(site).get_M().conjugate(), Textra::idx({0,2},{0,2}) );
//    if(not Textra::TensorMatrixMap(rightID).isIdentity(1e-12)) {
//        std::cout << "L site   : \n" << state.get_mps(site).get_L() << std::endl;
//        std::cout << "L site+1 : \n" << state.get_mps(site + 1).get_L() << std::endl;
//
//        std::cout << "rightID: \n" << rightID << std::endl;
//        throw std::runtime_error(fmt::format("Not right normalized at site {} with threshold 1e-12", site));
//    }
//    tools::log->trace("SVD site {:2} log₁₀ trunc: {:12.8f} χlim: {:4} χ: {:4}", site, std::log10(state.get_truncation_error(site)),state.get_chi_lim(),
//                      state.get_mps(site).get_chiR());
//    tools::common::profile::t_svd->toc();
//
//
//}
//
//
//void tools::finite::svd::truncate_left(const Eigen::Tensor<std::complex<double>,3> &theta, class_state_finite & state, std::optional<size_t> chi_lim, std::optional<double> svd_threshold){
//    tools::log->trace("Truncating multitheta from right to left (used when direction == 1 == ---> )");
//    if(not svd_threshold.has_value()){
//        svd_threshold = settings::precision::svd_threshold;
//    }
//
//    class_SVD SVD;
//    SVD.setThreshold(svd_threshold.value());
//    using Scalar = class_state_finite::Scalar;
//    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
//    if(not chi_lim.has_value()) chi_lim = state.get_chi_lim();
//    Eigen::Tensor<Scalar,4> theta4;
//    Eigen::Tensor<Scalar,3> U = theta;
//    Eigen::Tensor<Scalar,3> V;
//    Eigen::Tensor<Scalar,1> S;
//    auto reverse_active_sites = state.active_sites;
//    std::reverse(reverse_active_sites.begin(),reverse_active_sites.end());
//    double norm;
//    tools::common::profile::t_svd->tic();
//    while (reverse_active_sites.size() >= 2){
//        size_t site = reverse_active_sites.front();
//        long dim0 = U.dimension(0) / state.get_mps(site).spin_dim();
//        long dim1 = state.get_mps(site).spin_dim();
//        long dim2 = U.dimension(1);
//        long dim3 = U.dimension(2);
//        theta4 = U
//                .reshape(Textra::array4{dim0,dim1,dim2,dim3})
//                .shuffle(Textra::array4{0,2,1,3});
//        try {
//            std::tie(U, S, V,norm) = SVD.schmidt_with_norm(theta4, chi_lim);
//        }
//        catch(std::exception &ex){
//            std::cerr << "U :\n" << U << std::endl;
//            std::cerr << "S :\n" << S << std::endl;
//            std::cerr << "V :\n" << V << std::endl;
//            std::cerr << "theta4:\n" << theta4 << std::endl;
//            throw std::runtime_error(fmt::format("Truncation failed at site {}: {}", site, ex.what()));
//        }
//
//
//        state.set_truncation_error(site-1,SVD.get_truncation_error());
//        state.get_mps(site).set_M(V);
//        state.get_mps(site).unset_LC();
//
//
//
//        if(reverse_active_sites.size() >= 3){
//            Eigen::Tensor<Scalar,3> temp =  U.contract(Textra::asDiagonal(S), Textra::idx({2},{0}));
//            U = temp;
//            state.get_mps(site - 1).set_L(S);
//
//            if(state.ENV_R. front().get_position() != site) throw std::runtime_error(fmt::format("Site and postion mismatch in ENV_R while truncating right to left {} != {}", state.ENV_R. front().get_position() , site));
//            if(state.ENV2_R.front().get_position() != site) throw std::runtime_error(fmt::format("Site and postion mismatch in ENV2_R while truncating right to left {} != {}",state.ENV2_R.front().get_position() , site));
////            class_environment     R  = state.ENV_R.front();
////            class_environment_var R2 = state.ENV2_R.front();
////            R.enlarge (state.get_mps(site), state.get_mpo(site));
////            R2.enlarge(state.get_mps(site), state.get_mpo(site));
//            state.ENV_R .emplace_front(state.ENV_R .front().enlarge(state.get_mps(site), state.get_MPO(site)));
//            state.ENV2_R.emplace_front(state.ENV2_R.front().enlarge(state.get_mps(site), state.get_MPO(site)));
//
//        }else{
//            state.get_mps(site - 1).set_LC(S);
//        }
//        tools::log->trace("SVD site {:2} log₁₀ trunc: {:12.8f} χlim: {:4} χ: {:4}", site, std::log10(state.get_truncation_error(site)),state.get_chi_lim(),
//                          state.get_mps(site).get_chiR());
//        reverse_active_sites.pop_front();
//
//        if (not Eigen::Map<VectorType>(V.data(),V.size()).allFinite() )
//            throw std::runtime_error("V_L has nan's or inf's");
//
//        Eigen::Tensor<Scalar,2> rightID = state.get_mps(site).get_M()
//                .contract(state.get_mps(site).get_M().conjugate(), Textra::idx({0,2},{0,2}) );
//        if(not Textra::TensorMatrixMap(rightID).isIdentity(1e-12)) throw std::runtime_error(fmt::format("Not right normalized at site {} with threshold 1e-12", site));
//
//    }
//    size_t site = reverse_active_sites.front();
////    Eigen::Tensor<Scalar,3> L_U = Textra::asDiagonalInversed(state.get_L(site)).contract(U,Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
////    state.get_G(site) = L_U;
//    state.get_mps(site).set_M(U);
//    if (not Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1 >>(U.data(),U.size()).allFinite() )
//        throw std::runtime_error("L_U has nan's or inf's");
//
//    Eigen::Tensor<Scalar,2> leftID = state.get_mps(site).get_M_bare()
//            .contract(state.get_mps(site).get_M_bare().conjugate(), Textra::idx({0,1},{0,1}) );
//    if(not Textra::TensorMatrixMap(leftID).isIdentity(1e-12)) throw std::runtime_error(fmt::format("Not left normalized at site {} with threshold 1e-12", site));
//    tools::log->trace("SVD site {:2} log₁₀ trunc: {:12.8f} χlim: {:4} χ: {:4}", site, std::log10(state.get_truncation_error(site)),state.get_chi_lim(),
//                      state.get_mps(site).get_chiR());
//    tools::common::profile::t_svd->toc();
//
//}
//
//
//void tools::finite::svd::truncate_theta(const Eigen::Tensor<std::complex<double>,4> &theta, class_state_finite & state, std::optional<size_t> chi_lim, std::optional<double> svd_threshold) {
//    tools::common::profile::t_svd->tic();
//    if(not chi_lim.has_value()) chi_lim = state.get_chi_lim();
//    if(not svd_threshold.has_value()) svd_threshold = settings::precision::svd_threshold;
//
//    // Calculate the minimum chi for this position
//    // We always need to make sure we aren't truncating too aggressively.
//    // We can only truncate down by a factor of d = {spin dimension} compared to adjacent sites.
//
//    size_t chiL = theta.dimension(1);
//    size_t chiR = theta.dimension(3);
//    size_t dimL = theta.dimension(0);
//    size_t dimR = theta.dimension(2);
//    size_t chi_min = (size_t) std::ceil(std::max(chiL,chiR)/std::min(dimL,dimR));
//    chi_min = std::max(chi_min,chi_lim.value());
//    if(chi_min == 0 or chi_min > 10000) throw std::logic_error("Chi lim out of bounds for truncation: " + std::to_string(chi_min));
//
//    class_SVD SVD;
//    SVD.setThreshold(svd_threshold.value());
//    auto[U, S, V] = SVD.schmidt(theta, chi_min);
//
//    if(chi_lim.value() <= (size_t) state.get_chi_lim())
//        state.set_truncation_error(SVD.get_truncation_error());
//
//    state.MPS_L.back().set_M(U);
//    state.MPS_L.back().set_LC(S);
//    state.MPS_R.front().set_M(V);
//    state.clear_measurements();
//    tools::common::profile::t_svd->toc();
//}
//
//
//void tools::finite::svd::truncate_all_sites(class_state_finite & state, const size_t & chi_lim, std::optional<double> svd_threshold){
//    tools::log->trace("Truncating all sites to bond dimension {}",chi_lim);
//
//    auto original_position = state.get_position();
//    auto original_direction = state.get_direction();
//    // Start by truncating at the current position.
//    while(true){
//        move_center_point(state,chi_lim,svd_threshold);
//        if(state.get_position() == original_position and state.get_direction() == original_direction){
//            // Check if all bond dimensions less than or equal to below chi_lim
//            auto bond_dimensions = tools::finite::measure::bond_dimensions(state);
//            if (std::all_of(bond_dimensions.begin(), bond_dimensions.end(),
//                           [chi_lim](const size_t &chi){return chi <= chi_lim;}))
//                break;
//        }
//    }
//
//    tools::log->trace("Truncated all sites");
//}
//
//
//void tools::finite::svd::truncate_active_sites(class_state_finite & state, const size_t & chi_lim, std::optional<double> svd_threshold){
//    if(state.active_sites.empty()) return;
//    if(state.active_sites.size() < 2) throw std::runtime_error("Need 2 or more active sites");
//
//    tools::log->trace("Currently active sites {}", state.active_sites);
//
//    auto target_sites = state.active_sites;
//    // When going to the right, we shouldn't truncate the last site since
//    // the corresponding LC will not be part of the update.
//    target_sites.pop_back();
//
//    if(target_sites.size() <= 1){
//        //We only got two active, sites, this is simple.
//        tools::finite::svd::truncate_theta(state.get_theta(),state,chi_lim,svd_threshold);
//        return;
//    }
//    tools::log->trace("Truncating interior sites {}", target_sites);
//    tools::log->trace("Bond dimensions: {}", tools::finite::measure::bond_dimensions(state));
//
//
//    size_t original_position  = state.get_position();
//    int original_direction    = state.get_direction();
//    auto last_bond_dimensions = tools::finite::measure::bond_dimensions(state);
//    while(true){
//        size_t next_position = state.get_position() + state.get_direction();
//        bool next_position_is_interior = std::find(target_sites.begin(), target_sites.end(), next_position) != target_sites.end();
//        if(next_position_is_interior)
//            tools::finite::mps::move_center_point(state,chi_lim);
//        else
//            tools::finite::mps::move_center_point(state,state.get_chi_lim());
//
//        //Check if we are done
//        if(state.get_position() == original_position and state.get_direction() == original_direction){
//            // Check if all bond dimensions less than or equal to below chi_lim
//            auto prev_bond_dimensions = last_bond_dimensions;
//            last_bond_dimensions = tools::finite::measure::bond_dimensions(state);
//            if (last_bond_dimensions == prev_bond_dimensions) break;
//        }
//
//        //We don't need to go too far away from the active sites, you can come back earlier
//        if(state.get_direction() > 0 and state.get_position()  >  state.active_sites.back()) state.flip_direction();
//        if(state.get_direction() < 0 and state.get_position()  <  state.active_sites.front()) state.flip_direction();
//
////        tools::log->info("Pos {} dir {} active_sites {}", state.get_position(), state.get_direction(),state.active_sites);
//    }
//    tools::log->trace("Bond dimensions: {}", tools::finite::measure::bond_dimensions(state));
//
//}
//
//
//void tools::finite::svd::truncate_next_sites(class_state_finite & state, const size_t & chi_lim, size_t num_sites, std::optional<double> svd_threshold){
//    state.clear_cache();
//    state.clear_measurements();
//    if(not state.active_sites.empty()) throw std::logic_error("Make sure the active site list is empty before reaching this point");
//
//    size_t from = state.get_direction() > 0 ? state.get_position() : state.get_position() + 1;
//    size_t upto = from + (static_cast<size_t>(state.get_direction()) * (num_sites - 1));
//    from = std::max(from,0ul);
//    from = std::min(from,state.get_length()-1);
//    upto = std::max(upto,0ul);
//    upto = std::min(upto,state.get_length()-1);
//    tools::log->trace("Truncating starting from {} upto {}, num sites = {}",from,upto,num_sites);
//
//    auto sites  =  math::range(std::min(from,upto),std::max(from,upto),1);
//    state.active_sites = std::list<size_t>(sites.begin(),sites.end());
//    truncate_active_sites(state,chi_lim,svd_threshold);
//    state.active_sites.clear();
//}
//

