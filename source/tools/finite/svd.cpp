//
// Created by david on 2019-06-29.
//
#include <tools/finite/opt.h>
#include <state/class_state_finite.h>
#include <math/class_svd_wrapper.h>
#include <simulation/nmspc_settings.h>


void tools::finite::mps::normalize(class_state_finite & state){

    tools::log->trace("Normalizing state");
    using namespace Textra;
    using Scalar = class_state_finite::Scalar;
    state.unset_measurements();
    tools::common::profile::t_svd.tic();
    size_t num_moves = 2*(state.get_length()-2);
    tools::log->trace("Norm before normalization = {:.16f}", tools::finite::measure::norm(state));
    tools::log->trace("Bond dimensions before normalization: {}", tools::finite::measure::bond_dimensions(state));

    for(size_t moves = 0; moves < num_moves; moves++){
        auto & MPS_L  = state.MPS_L;
        auto & MPS_R  = state.MPS_R;
        if(MPS_L.empty()) throw std::runtime_error("MPS_L is empty");
        if(MPS_R.empty()) throw std::runtime_error("MPS_R is empty");

        //Store the special LC bond in a temporary.
        Eigen::Tensor<Scalar,1> LC = MPS_L.back().get_LC();
        MPS_L.back().unset_LC();
        if (state.get_direction() == 1){
            MPS_L.emplace_back(class_mps_site(MPS_R.front().get_M(), LC, MPS_R.front().get_position()));
            MPS_R.pop_front();
            Eigen::Tensor<Scalar,4> theta =
                    Textra::asDiagonal(LC)
                            .contract(state.MPS_L.back().get_M(), Textra::idx({1},{1}))
                            .contract(state.MPS_R.front().get_M(), Textra::idx({2},{1}))
                            .shuffle(Textra::array4{1,0,2,3});
            tools::finite::opt::truncate_theta(theta, state);
        }else{
            MPS_R.emplace_front(class_mps_site(MPS_L.back().get_M(), LC, MPS_L.back().get_position()));
            MPS_L.pop_back();
            Eigen::Tensor<Scalar,4> theta =
                    state.MPS_L.back().get_M()
                            .contract(state.MPS_R.front().get_M(), Textra::idx({2},{1}))
                            .contract(Textra::asDiagonal(LC), Textra::idx({3},{0}));
            tools::finite::opt::truncate_theta(theta, state);
        }

        if(MPS_L.empty()) throw std::runtime_error("MPS_L became empty");
        if(MPS_R.empty()) throw std::runtime_error("MPS_R became empty");
        if(MPS_L.size() + MPS_R.size() != state.get_length()) throw std::runtime_error("MPS_L + MPS_R sizes do not add up to chain length anymore");
        //    Check edge
        if (state.position_is_any_edge()){
            state.flip_direction();
        }
    }
    tools::common::profile::t_svd.toc();
    state.unset_measurements();
    tools::log->trace("Norm after normalization = {:.16f}", tools::finite::measure::norm(state));
    tools::log->trace("Bond dimensions after  normalization: {}", tools::finite::measure::bond_dimensions(state));

}


void tools::finite::opt::truncate_theta(Eigen::Tensor<Scalar,3> &theta, class_state_finite & state){
    state.unset_measurements();
    if (state.active_sites.empty()) throw std::runtime_error("truncate_theta: No active sites to truncate");
    if (theta.size() == 0)          throw std::runtime_error("truncate_theta: Theta is empty");
    auto theta_map = Eigen::Map<Eigen::VectorXcd>(theta.data(),theta.size());
    auto fullnorm  = tools::finite::measure::norm(state);
    auto thetanorm = theta_map.norm();
    if(std::abs(fullnorm  - 1.0) > settings::precision::maxNormError){ tools::log->error("Norm before truncation too far from unity: {:.16f}", fullnorm);}
    if(std::abs(thetanorm - 1.0) > settings::precision::maxNormError){ tools::log->error("Norm of theta too far from unity: {:.16f}", thetanorm); theta_map.normalize();}

    //Start by clearing the environments that will become stale.
    while(true){
        if(state.active_sites.back() > state.ENV_R.front().get_position()){
            state.ENV_R.pop_front();
            state.ENV2_R.pop_front();
        }
        else if(state.active_sites.front() < state.ENV_L.back().get_position()){
            state.ENV_L.pop_back();
            state.ENV2_L.pop_back();
        }
        else break;
    }

    if (state.get_direction() == 1){
        tools::finite::opt::truncate_left(theta,state);

    }else{
        tools::finite::opt::truncate_right(theta,state);

    }
    state.unset_measurements();
    state.clear_cache();
    if (state.ENV_L.size() + state.ENV_R.size() != state.get_length())
        throw std::runtime_error(fmt::format("Environment sizes do not add up to system size after truncation: {} + {} != {}", state.ENV_L.size() , state.ENV_R.size(), state.get_length()));
    if (state.ENV2_L.size() + state.ENV2_R.size() != state.get_length())
        throw std::runtime_error(fmt::format("Environment sq sizes do not add up to system size after truncation: {} + {} != {}", state.ENV2_L.size() , state.ENV2_R.size(), state.get_length()));

    fullnorm  = tools::finite::measure::norm(state);
    if(std::abs(fullnorm  - 1.0) > 1e-10){
        tools::log->warn("Norm after truncation too far from unity: {:.16f} -- Normalizing",fullnorm);
        tools::finite::mps::normalize(state);
        fullnorm  = tools::finite::measure::norm(state);
        tools::log->warn("New norm: {:.16f}",fullnorm);
    }


}



void tools::finite::opt::truncate_right(Eigen::Tensor<std::complex<double>,3> &theta, class_state_finite & state){
    tools::log->trace("Truncating multitheta from left to right");
    class_SVD SVD;
    SVD.setThreshold(settings::precision::SVDThreshold);
    using Scalar = class_state_finite::Scalar;
    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    Eigen::Tensor<Scalar,4> theta4;
    Eigen::Tensor<Scalar,3> U;
    Eigen::Tensor<Scalar,3> V = theta;
    Eigen::Tensor<Scalar,1> S;
    auto active_sites = state.active_sites;
    double norm;
    tools::common::profile::t_svd.tic();
    while (active_sites.size() >= 2){
        size_t site = active_sites.front();
        long dim0 = state.get_MPS(site).get_spin_dim();
        long dim1 = V.dimension(0) /  state.get_MPS(site).get_spin_dim();
        long dim2 = V.dimension(1);
        long dim3 = V.dimension(2);
        theta4 = V
                .reshape(Textra::array4{dim0,dim1,dim2,dim3})
                .shuffle(Textra::array4{0,2,1,3});
        std::tie(U, S, V,norm) = SVD.schmidt_with_norm(theta4,state.get_chi_lim());
        state.set_truncation_error(site, SVD.get_truncation_error());

        state.get_MPS(site).set_M(U);
        state.get_MPS(site).unset_LC();


        if (not Eigen::Map<VectorType>(U.data(),U.size()).allFinite() )
            throw std::runtime_error("L_U has nan's or inf's");

        Eigen::Tensor<Scalar,2> leftID = state.get_MPS(site).get_M()
                .contract(state.get_MPS(site).get_M().conjugate(), Textra::idx({0,1},{0,1}) );
        if(not Textra::TensorMatrixMap(leftID).isIdentity(1e-12)) throw std::runtime_error(fmt::format("Not left normalized at site {} with threshold 1e-12.", site));



        if(active_sites.size() >= 3){
            Eigen::Tensor<Scalar,3> temp = Textra::asDiagonal(S).contract(V, Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
            state.get_MPS(site+1).set_L(S);
            V = temp;
            if(state.ENV_L. back().get_position() != site) throw std::runtime_error(fmt::format("Site and postion mismatch in ENV_L while truncating  left to right {} != {}",state.ENV_L. front().get_position() , site));
            if(state.ENV2_L.back().get_position() != site) throw std::runtime_error(fmt::format("Site and postion mismatch in ENV2_L while truncating  left to right{} != {}",state.ENV2_L.front().get_position() , site));
            class_environment     L  = state.ENV_L.back();
            class_environment_var L2 = state.ENV2_L.back();
            L.enlarge (state.get_MPS(site), state.get_MPO(site).MPO());
            L2.enlarge(state.get_MPS(site), state.get_MPO(site).MPO());
            state.ENV_L.emplace_back(L);
            state.ENV2_L.emplace_back(L2);
        } else{
            //Always set LC on the last "A" matrix
            state.get_MPS(site).set_LC(S);
        }
        tools::log->trace("Site {:2} log₁₀ trunc: {:12.8f} χlim: {:4} χ: {:4}", site, std::log10(state.get_truncation_error(site)),state.get_chi_lim(), state.get_MPS(site).get_chiR());
        active_sites.pop_front();
    }


    size_t site = active_sites.front();
//    Eigen::Tensor<Scalar,3> V_L = V.contract(Textra::asDiagonalInversed(state.get_L(site+1)), Textra::idx({2},{0}));
//    state.get_G(site) = V_L;
    state.get_MPS(site).set_M(V);

    if (not Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1 >>(V.data(),V.size()).allFinite() )
        throw std::runtime_error("V_L has nan's or inf's");


    Eigen::Tensor<Scalar,2> rightID = state.get_MPS(site).get_M()
            .contract(state.get_MPS(site).get_M().conjugate(), Textra::idx({0,2},{0,2}) );
    if(not Textra::TensorMatrixMap(rightID).isIdentity(1e-12)) {
        std::cout << "L site   : \n" << state.get_MPS(site).get_L() << std::endl;
        std::cout << "L site+1 : \n" << state.get_MPS(site+1).get_L() << std::endl;

        std::cout << "rightID: \n" << rightID << std::endl;
        throw std::runtime_error(fmt::format("Not right normalized at site {} with threshold 1e-12", site));
    }
    tools::log->trace("Site {:2} log₁₀ trunc: {:12.8f} χlim: {:4} χ: {:4}", site, std::log10(state.get_truncation_error(site)),state.get_chi_lim(), state.get_MPS(site).get_chiR());
    tools::common::profile::t_svd.toc();


}


void tools::finite::opt::truncate_left(Eigen::Tensor<std::complex<double>,3> &theta, class_state_finite & state){
    tools::log->trace("Truncating multitheta from right to left");
    class_SVD SVD;
    SVD.setThreshold(settings::precision::SVDThreshold);
    using Scalar = class_state_finite::Scalar;
    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    Eigen::Tensor<Scalar,4> theta4;
    Eigen::Tensor<Scalar,3> U = theta;
    Eigen::Tensor<Scalar,3> V;
    Eigen::Tensor<Scalar,1> S;
    auto reverse_active_sites = state.active_sites;
    std::reverse(reverse_active_sites.begin(),reverse_active_sites.end());
    double norm;
    tools::common::profile::t_svd.tic();
    while (reverse_active_sites.size() >= 2){
        size_t site = reverse_active_sites.front();
        long dim0 = U.dimension(0) /  state.get_MPS(site).get_spin_dim();
        long dim1 = state.get_MPS(site).get_spin_dim();
        long dim2 = U.dimension(1);
        long dim3 = U.dimension(2);
        theta4 = U
                .reshape(Textra::array4{dim0,dim1,dim2,dim3})
                .shuffle(Textra::array4{0,2,1,3});
        try {std::tie(U,S,V,norm) = SVD.schmidt_with_norm(theta4, state.get_chi_lim());}
        catch(std::exception &ex){
            std::cerr << "U :\n" << U << std::endl;
            std::cerr << "S :\n" << S << std::endl;
            std::cerr << "V :\n" << V << std::endl;
            std::cerr << "theta4:\n" << theta4 << std::endl;
            throw std::runtime_error(fmt::format("Truncation failed at site {}: {}", site, ex.what()));
        }


        state.set_truncation_error(site-1,SVD.get_truncation_error());
        state.get_MPS(site).set_M(V);
        state.get_MPS(site).unset_LC();



        if(reverse_active_sites.size() >= 3){
            Eigen::Tensor<Scalar,3> temp =  U.contract(Textra::asDiagonal(S), Textra::idx({2},{0}));
            U = temp;
            state.get_MPS(site-1).set_L(S);

            if(state.ENV_R. front().get_position() != site) throw std::runtime_error(fmt::format("Site and postion mismatch in ENV_R while truncating right to left {} != {}", state.ENV_R. front().get_position() , site));
            if(state.ENV2_R.front().get_position() != site) throw std::runtime_error(fmt::format("Site and postion mismatch in ENV2_R while truncating right to left {} != {}",state.ENV2_R.front().get_position() , site));
            class_environment     R  = state.ENV_R.front();
            class_environment_var R2 = state.ENV2_R.front();
            R.enlarge (state.get_MPS(site), state.get_MPO(site).MPO());
            R2.enlarge(state.get_MPS(site), state.get_MPO(site).MPO());
            state.ENV_R.emplace_front(R);
            state.ENV2_R.emplace_front(R2);

        }else{
            state.get_MPS(site-1).set_LC(S);
        }
        tools::log->trace("Site {:2} log₁₀ trunc: {:12.8f} χlim: {:4} χ: {:4}", site, std::log10(state.get_truncation_error(site)),state.get_chi_lim(), state.get_MPS(site).get_chiR());
        reverse_active_sites.pop_front();

        if (not Eigen::Map<VectorType>(V.data(),V.size()).allFinite() )
            throw std::runtime_error("V_L has nan's or inf's");

        Eigen::Tensor<Scalar,2> rightID = state.get_MPS(site).get_M()
                .contract(state.get_MPS(site).get_M().conjugate(), Textra::idx({0,2},{0,2}) );
        if(not Textra::TensorMatrixMap(rightID).isIdentity(1e-12)) throw std::runtime_error(fmt::format("Not right normalized at site {} with threshold 1e-12", site));

    }
    size_t site = reverse_active_sites.front();
//    Eigen::Tensor<Scalar,3> L_U = Textra::asDiagonalInversed(state.get_L(site)).contract(U,Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
//    state.get_G(site) = L_U;
    state.get_MPS(site).set_M(U);
    if (not Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1 >>(U.data(),U.size()).allFinite() )
        throw std::runtime_error("L_U has nan's or inf's");

    Eigen::Tensor<Scalar,2> leftID = state.get_MPS(site).get_M_bare()
            .contract(state.get_MPS(site).get_M_bare().conjugate(), Textra::idx({0,1},{0,1}) );
    if(not Textra::TensorMatrixMap(leftID).isIdentity(1e-12)) throw std::runtime_error(fmt::format("Not left normalized at site {} with threshold 1e-12", site));
    tools::log->trace("Site {:2} log₁₀ trunc: {:12.8f} χlim: {:4} χ: {:4}", site, std::log10(state.get_truncation_error(site)),state.get_chi_lim(), state.get_MPS(site).get_chiR());
    tools::common::profile::t_svd.toc();

}


void tools::finite::opt::truncate_theta(Eigen::Tensor<std::complex<double>,4> &theta, class_state_finite & state) {
    tools::common::profile::t_svd.tic();
    class_SVD SVD;
    SVD.setThreshold(settings::precision::SVDThreshold);
    auto[U, S, V] = SVD.schmidt(theta, state.get_chi_lim());
    state.set_truncation_error(SVD.get_truncation_error());
    state.MPS_L.back().set_M(U);
    state.MPS_L.back().set_LC(S);
    state.MPS_R.front().set_M(V);
    state.unset_measurements();
    tools::common::profile::t_svd.toc();
}



int tools::finite::mps::move_center_point(class_state_finite & state){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading
//    std::cout << "Current state -- Direction: " << direction << std::endl;
//    std::cout << "HA: " << state.HA->get_position() << " MPO_L back : " << MPO_L.back()->get_position() << std::endl;
//    std::cout << "HB: " << state.HB->get_position() << " MPO_R front: " << MPO_R.front()->get_position() << std::endl;
//

    auto & MPS_L  = state.MPS_L;
    auto & MPS_R  = state.MPS_R;
    auto & MPO_L  = state.MPO_L;
    auto & MPO_R  = state.MPO_R;
    auto & ENV_L  = state.ENV_L;
    auto & ENV_R  = state.ENV_R;
    auto & ENV2_L = state.ENV2_L;
    auto & ENV2_R = state.ENV2_R;
    if(ENV_L.empty()) throw std::runtime_error("ENVL is empty");
    if(ENV_R.empty()) throw std::runtime_error("ENVR is empty");
    if(MPS_L.empty()) throw std::runtime_error("MPSL is empty");
    if(MPS_R.empty()) throw std::runtime_error("MPSR is empty");
    if(MPS_L.back().get_position()  != ENV_L.back().get_position())  throw std::runtime_error("MPSL and ENVL have mismatching positions");
    if(MPS_R.front().get_position() != ENV_R.front().get_position()) throw std::runtime_error("MPSR and ENVR have mismatching positions");
    if(ENV_L.size() + ENV_R.size() != state.get_length()) throw std::runtime_error("ENVL + ENVR sizes do not add up to chain length");
    if(MPS_L.size() + MPS_R.size() != state.get_length()) throw std::runtime_error("MPSL + MPSR sizes do not add up to chain length");
    assert(ENV_L.size() + ENV_R.size() == state.get_length());
    assert(ENV_L.back().sites + ENV_R.front().sites == state.get_length() - 2);
    //Store the special LC bond in a temporary.
    Eigen::Tensor<Scalar,1> LC = MPS_L.back().get_LC();
    MPS_L.back().unset_LC();

    if (state.get_direction() == 1){
        class_environment     L  = ENV_L.back();
        class_environment_var L2 = ENV2_L.back();
        L.enlarge(MPS_L.back(), MPO_L.back()->MPO());
        L2.enlarge(MPS_L.back(), MPO_L.back()->MPO());
        ENV_L.emplace_back(L);
        ENV2_L.emplace_back(L2);
        MPS_L.emplace_back(class_mps_site(MPS_R.front().get_M(), LC, MPS_R.front().get_position()));
        MPO_L.emplace_back(MPO_R.front()->clone());
        MPS_R.pop_front();
        MPO_R.pop_front();
        ENV_R.pop_front();
        ENV2_R.pop_front();
        Eigen::Tensor<Scalar,4> theta =
                Textra::asDiagonal(LC)
                .contract(state.MPS_L.back().get_M(), Textra::idx({1},{1}))
                .contract(state.MPS_R.front().get_M(), Textra::idx({2},{1}))
                .shuffle(Textra::array4{1,0,2,3});
        tools::finite::opt::truncate_theta(theta,state);
    }else{

        class_environment     R  = ENV_R.front();
        class_environment_var R2 = ENV2_R.front();
        R.enlarge(MPS_R.front(), MPO_R.front()->MPO());
        R2.enlarge(MPS_R.front(), MPO_R.front()->MPO());
        ENV_R.emplace_front(R);
        ENV2_R.emplace_front(R2);
        MPS_R.emplace_front(class_mps_site(MPS_L.back().get_M(), LC, MPS_L.back().get_position()));
        MPO_R.emplace_front(MPO_L.back()->clone());
        MPS_L.pop_back();
        MPO_L.pop_back();
        ENV_L.pop_back();
        ENV2_L.pop_back();
        Eigen::Tensor<Scalar,4> theta =
                state.MPS_L.back().get_M()
                .contract(state.MPS_R.front().get_M(), Textra::idx({2},{1}))
                .contract(Textra::asDiagonal(LC), Textra::idx({3},{0}));
        tools::finite::opt::truncate_theta(theta,state);
    }

    assert(MPO_L.size() + MPO_R.size() == state.get_length());
    if(ENV_L.empty()) throw std::runtime_error("ENVL became empty");
    if(ENV_R.empty()) throw std::runtime_error("ENVR became empty");
    if(MPS_L.empty()) throw std::runtime_error("MPSL became empty");
    if(MPS_R.empty()) throw std::runtime_error("MPSR became empty");
    if(MPS_L.back().get_position()  != ENV_L.back().get_position())  throw std::runtime_error("MPSL and ENVL got mismatching positions");
    if(MPS_R.front().get_position() != ENV_R.front().get_position()) throw std::runtime_error("MPSR and ENVR got mismatching positions");
    if(ENV_L.size() + ENV_R.size() != state.get_length()) throw std::runtime_error("ENVL + ENVR sizes do not add up to chain length anymore");
    if(MPS_L.size() + MPS_R.size() != state.get_length()) throw std::runtime_error("MPSL + MPSR sizes do not add up to chain length anymore");
    //    Check edge
    if (state.position_is_any_edge()){
        state.flip_direction();
        state.increment_sweeps();
    }

    state.increment_moves();
    state.clear_cache();
    state.unset_measurements();
    state.active_sites.clear();

    tools::finite::debug::check_integrity(state);

    return state.get_sweeps();
}
