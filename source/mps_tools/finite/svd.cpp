//
// Created by david on 2019-06-29.
//
#include <mps_tools/finite/opt.h>
#include <mps_state/class_finite_chain_state.h>
#include <general/class_svd_wrapper.h>
#include <sim_parameters/nmspc_sim_settings.h>

void mpstools::finite::mps::normalize(class_finite_chain_state & state){

    mpstools::log->trace("Normalizing state");
    using namespace Textra;
    using Scalar = class_finite_chain_state::Scalar;
    state.unset_measurements();
//    std::cout << "Norm              (before normalization): " << mpstools::finite::measure::norm(state)  << std::endl;
//    std::cout << "Spin component sx (before normalization): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sx)  << std::endl;
//    std::cout << "Spin component sy (before normalization): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sy)  << std::endl;
//    std::cout << "Spin component sz (before normalization): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sz)  << std::endl;
    // Sweep back and forth once on the state

    class_SVD svd;
    svd.setThreshold(settings::precision::SVDThreshold);

    size_t pos_LA = 0; // Lambda left of GA
    size_t pos_LC = 1; //Center Lambda
    size_t pos_LB = 2; // Lambda right of GB
    size_t pos_A  = 0; // Lambda * G
    size_t pos_B  = 1; // G * Lambda
//    bool finished = false;
    size_t num_traversals = 0;
    size_t step = 0;
    int direction = 1;
    Eigen::Tensor<Scalar,3> U;
    Eigen::Tensor<Scalar,1> S;
    Eigen::Tensor<Scalar,3> V;
    double norm;
    while(num_traversals < 4){
        Eigen::Tensor<Scalar,4> theta = state.get_theta(pos_A);
        try {std::tie(U,S,V,norm) = svd.schmidt_with_norm(theta);}
        catch(std::exception &ex){
            std::cerr << "A " << pos_A  << ":\n" << state.get_A(pos_A) << std::endl;
            std::cerr << "L " << pos_LC << ":\n" << state.get_L(pos_LC) << std::endl;
            std::cerr << "B " << pos_B  << ":\n" << state.get_B(pos_B) << std::endl;
            std::cerr << "theta:\n" << theta << std::endl;
            throw std::runtime_error(fmt::format("Normalization failed at positions A:{} C:{} B:{} , step {}: {}", pos_A, pos_LC, pos_B, step, ex.what()));
        }
        Eigen::Tensor<Scalar,3> LA_U = Textra::asDiagonalInversed(state.get_L(pos_LA)).contract(U,idx({1},{1})).shuffle(array3{1,0,2});
        Eigen::Tensor<Scalar,3> V_LB = V.contract(asDiagonalInversed(state.get_L(pos_LB)), idx({2},{0}));
        norm = pos_B == state.get_length()-1 or pos_A == 0 ? 1.0 : norm;
        if (direction == 1){
            state.get_G(pos_A)  = LA_U;
            state.get_L(pos_LC) = S;
            state.get_G(pos_B)  = Scalar(norm)*V_LB;

        }
        if (direction == -1){
            state.get_G(pos_A)  = LA_U * Scalar(norm);
            state.get_L(pos_LC) = S;
            state.get_G(pos_B)  = V_LB;
        }


        if (direction ==  1 and pos_B == state.get_length()-1)  {num_traversals++; direction *= -1;}
        if (direction == -1 and pos_A == 0)                     {num_traversals++; direction *= -1;}

        pos_LA += direction;
        pos_LC += direction;
        pos_LB += direction;
        pos_A  += direction;
        pos_B  += direction;
        step   ++;


    }
    state.unset_measurements();
//    std::cout << "Norm              (after normalization): " << mpstools::finite::measure::norm(state)  << std::endl;
//    std::cout << "Spin component sx (after normalization): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sx)  << std::endl;
//    std::cout << "Spin component sy (after normalization): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sy)  << std::endl;
//    std::cout << "Spin component sz (after normalization): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sz)  << std::endl;


}


void mpstools::finite::opt::truncate_theta(Eigen::Tensor<std::complex<double>,3> &theta, class_finite_chain_state & state, long chi_, double SVDThreshold){
    if (state.get_direction() == 1){
        mpstools::finite::opt::truncate_right(theta,state,chi_,SVDThreshold);
    }else{
        mpstools::finite::opt::truncate_left(theta,state,chi_,SVDThreshold);
    }
}



void mpstools::finite::opt::truncate_right(Eigen::Tensor<std::complex<double>,3> &theta, class_finite_chain_state & state, long chi_, double SVDThreshold){
    mpstools::log->trace("Truncating multitheta from left to right");
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);
    state.unset_measurements();
    using Scalar = class_finite_chain_state::Scalar;
    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    auto theta_map = Eigen::Map<VectorType>(theta.data(),theta.size());
    auto fullnorm  = mpstools::finite::measure::norm(state);
    auto thetanorm = std::abs(theta_map.norm());
    if(std::abs(fullnorm - 1.0) > 1e-14) {
        throw std::runtime_error(fmt::format("Norm before truncation too far from unity: {}",fullnorm));
    }
    if( std::abs(thetanorm - 1.0) > 1e-14) {
        throw std::runtime_error(fmt::format("Norm of theta too far from unity: {}",thetanorm));
    }


    Eigen::Tensor<Scalar,4> theta4;
    Eigen::Tensor<Scalar,3> U;
    Eigen::Tensor<Scalar,3> V = theta;
    Eigen::Tensor<Scalar,1> S;
    auto active_sites = state.active_sites;
    double norm;

    while (active_sites.size() > 1){
        size_t site = active_sites.front();
        long dim0 = state.get_MPS(site).get_spin_dim();
        long dim1 = V.dimension(0) /  state.get_MPS(site).get_spin_dim();
        long dim2 = V.dimension(1);
        long dim3 = V.dimension(2);
        theta4 = V
                .reshape(Textra::array4{dim0,dim1,dim2,dim3})
                .shuffle(Textra::array4{0,2,1,3});
        std::tie(U, S, V,norm) = SVD.schmidt_with_norm(theta4,chi_);
        state.truncation_error[site+1] = SVD.get_truncation_error();
        Eigen::Tensor<Scalar,3> L_U = Textra::asDiagonalInversed(state.get_L(site)).contract(U,Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});

        state.get_G(site)   = L_U;
        state.get_L(site+1) = S;
        auto Snorm2_left  = Textra::Tensor1_to_Vector(state.get_L(site)).squaredNorm();
        auto Snorm2_right = Textra::Tensor1_to_Vector(state.get_L(site+1)).squaredNorm();
        mpstools::log->trace("svd multitheta site {} chi {} norm {:.16f} S norm2 left {:.16f} S norm2 right  {:.16f} trunc {:.16f}",
                             site,S.dimension(0), norm,Snorm2_left, Snorm2_right,SVD.get_truncation_error());

        if(active_sites.size() > 2){
            Eigen::Tensor<Scalar,3> temp = Textra::asDiagonal(S).contract(V, Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
            V = temp;
//            V = temp * Scalar(norm);

        }
        active_sites.pop_front();


        if (not Eigen::Map<VectorType>(L_U.data(),L_U.size()).allFinite() )
            throw std::runtime_error("L_U has nan's or inf's");

        Eigen::Tensor<Scalar,2> leftID = state.get_A(site)
                .contract(state.get_A(site).conjugate(), Textra::idx({0,1},{0,1}) );
        auto leftIDmap = Textra::Tensor2_to_Matrix(leftID);
        if(not leftIDmap.isIdentity(1e-14)) throw std::runtime_error(fmt::format("Not left normalized at site {}", site));

    }


    size_t site = active_sites.front();
    Eigen::Tensor<Scalar,3> V_L = V.contract(Textra::asDiagonalInversed(state.get_L(site+1)), Textra::idx({2},{0}));
    state.get_G(site) = V_L;


    if (not Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1 >>(V_L.data(),V_L.size()).allFinite() )
        throw std::runtime_error("V_L has nan's or inf's");


    Eigen::Tensor<Scalar,2> rightID = state.get_B(site)
            .contract(state.get_B(site).conjugate(), Textra::idx({0,2},{0,2}) );
    auto rightIDmap = Textra::Tensor2_to_Matrix(rightID);
    if(not rightIDmap.isIdentity(1e-14)) {
        std::cout << "L site   : \n" << state.get_L(site) << std::endl;
        std::cout << "L site+1 : \n" << state.get_L(site+1) << std::endl;

        std::cout << "rightID: \n" << rightID << std::endl;
        throw std::runtime_error(fmt::format("Not right normalized at site {}", site));
    }



    state.unset_measurements();
    fullnorm = mpstools::finite::measure::norm(state);
    if(std::abs(fullnorm - 1.0) > 1e-14) {
        throw std::runtime_error(fmt::format("Norm before rebuild of env too far from unity: {}",fullnorm));
    }
    state.unset_measurements();
    mpstools::finite::mps::rebuild_environments(state);

}




void mpstools::finite::opt::truncate_left(Eigen::Tensor<std::complex<double>,3> &theta, class_finite_chain_state & state, long chi_, double SVDThreshold){
    mpstools::log->trace("Truncating multitheta from right to left");
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);
    state.unset_measurements();
    using Scalar = class_finite_chain_state::Scalar;
    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    auto theta_map = Eigen::Map<VectorType>(theta.data(),theta.size());
    auto fullnorm  = mpstools::finite::measure::norm(state);
    auto thetanorm = std::abs(theta_map.norm());
    if(std::abs(fullnorm - 1.0) > 1e-14) {
        throw std::runtime_error(fmt::format("Norm before truncation too far from unity: {:.16f}",fullnorm));
    }
    if( std::abs(thetanorm - 1.0) > 1e-14) {
        throw std::runtime_error(fmt::format("Norm of theta too far from unity: {:.16f}",thetanorm));
    }


    Eigen::Tensor<Scalar,4> theta4;
    Eigen::Tensor<Scalar,3> U = theta;
    Eigen::Tensor<Scalar,3> V;
    Eigen::Tensor<Scalar,1> S;
    auto reverse_active_sites = state.active_sites;
    std::reverse(reverse_active_sites.begin(),reverse_active_sites.end());
    double norm;

    while (reverse_active_sites.size() > 1){
        size_t site = reverse_active_sites.front();
        long dim0 = U.dimension(0) /  state.get_MPS(site).get_spin_dim();
        long dim1 = state.get_MPS(site).get_spin_dim();
        long dim2 = U.dimension(1);
        long dim3 = U.dimension(2);
        theta4 = U
                .reshape(Textra::array4{dim0,dim1,dim2,dim3})
                .shuffle(Textra::array4{0,2,1,3});
        try {std::tie(U,S,V,norm) = SVD.schmidt_with_norm(theta4, chi_);}
        catch(std::exception &ex){
            std::cerr << "U :\n" << U << std::endl;
            std::cerr << "S :\n" << S << std::endl;
            std::cerr << "V :\n" << V << std::endl;
            std::cerr << "theta4:\n" << theta4 << std::endl;
            throw std::runtime_error(fmt::format("Truncation failed at site {}: {}", site, ex.what()));
        }


        state.truncation_error[site-1] = SVD.get_truncation_error();
        Eigen::Tensor<Scalar,3> V_L = V.contract(Textra::asDiagonalInversed(state.get_L(site+1)), Textra::idx({2},{0}));

        state.get_G(site) = V_L;
        state.get_L(site) = S;
        auto Snorm2_left  = Textra::Tensor1_to_Vector(state.get_L(site)).squaredNorm();
        auto Snorm2_right = Textra::Tensor1_to_Vector(state.get_L(site+1)).squaredNorm();
//        mpstools::log->trace("svd multitheta site {} chi {} norm {:.16f} S norm2 {:.16f} trunc {:.16f}", site,S.dimension(0), norm, Snorm2,SVD.get_truncation_error());
        mpstools::log->trace("svd multitheta site {} chi {} norm {:.16f} S norm2 left {:.16f} S norm2 right  {:.16f} trunc {:.16f}",
                             site,S.dimension(0), norm,Snorm2_left, Snorm2_right,SVD.get_truncation_error());

        if(reverse_active_sites.size() > 2){
            Eigen::Tensor<Scalar,3> temp =  U.contract(Textra::asDiagonal(S), Textra::idx({2},{0}));
//            U = Scalar(norm) * temp;
            U = temp;
        }
        reverse_active_sites.pop_front();

        if (not Eigen::Map<VectorType>(V_L.data(),V_L.size()).allFinite() )
            throw std::runtime_error("V_L has nan's or inf's");

        Eigen::Tensor<Scalar,2> rightID = state.get_B(site)
                .contract(state.get_B(site).conjugate(), Textra::idx({0,2},{0,2}) );
        auto rightIDmap = Textra::Tensor2_to_Matrix(rightID);
        if(not rightIDmap.isIdentity(1e-14)) throw std::runtime_error(fmt::format("Not right normalized at site {}", site));

    }
    size_t site = reverse_active_sites.front();
    Eigen::Tensor<Scalar,3> L_U = Textra::asDiagonalInversed(state.get_L(site)).contract(U,Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
    state.get_G(site) = L_U;
    auto Snorm2_left  = Textra::Tensor1_to_Vector(state.get_L(site)).squaredNorm();
    auto Snorm2_right = Textra::Tensor1_to_Vector(state.get_L(site+1)).squaredNorm();
//        mpstools::log->trace("svd multitheta site {} chi {} norm {:.16f} S norm2 {:.16f} trunc {:.16f}", site,S.dimension(0), norm, Snorm2,SVD.get_truncation_error());
    mpstools::log->trace("svd multitheta site {} chi {} norm {:.16f} S norm2 left {:.16f} S norm2 right  {:.16f} trunc {:.16f}",
                         site,S.dimension(0), norm,Snorm2_left, Snorm2_right,SVD.get_truncation_error());


    if (not Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1 >>(L_U.data(),L_U.size()).allFinite() )
        throw std::runtime_error("L_U has nan's or inf's");

    Eigen::Tensor<Scalar,2> leftID = state.get_MPS(site).get_A()
            .contract(state.get_MPS(site).get_A().conjugate(), Textra::idx({0,1},{0,1}) );
    auto leftIDmap = Textra::Tensor2_to_Matrix(leftID);
    if(not leftIDmap.isIdentity(1e-14)) throw std::runtime_error(fmt::format("Not left normalized at site {}", site));



    state.unset_measurements();
    fullnorm = mpstools::finite::measure::norm(state);
    if(std::abs(fullnorm - 1.0) > 1e-14) {
        throw std::runtime_error(fmt::format("Norm before rebuild of env too far from unity: {:.16f}",fullnorm));
    }
    state.unset_measurements();
    mpstools::finite::mps::rebuild_environments(state);
}


void mpstools::finite::opt::truncate_theta(Eigen::Tensor<std::complex<double>,4> &theta, class_finite_chain_state & state, long chi_, double SVDThreshold) {
    mpstools::log->trace("Truncating theta");
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);
    auto[U, S, V] = SVD.schmidt(theta, chi_);
    state.truncation_error[state.get_position()+1] = SVD.get_truncation_error();
    state.MPS_C  = S;
    Eigen::Tensor<std::complex<double>,3> L_U = Textra::asDiagonalInversed(state.MPS_L.back().get_L()).contract(U,Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
    Eigen::Tensor<std::complex<double>,3> V_L = V.contract(Textra::asDiagonalInversed(state.MPS_R.front().get_L()), Textra::idx({2},{0}));
    state.MPS_L.back().set_G(L_U);
    state.MPS_R.front().set_G(V_L);
}

