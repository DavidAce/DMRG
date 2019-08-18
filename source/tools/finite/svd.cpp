//
// Created by david on 2019-06-29.
//
#include <tools/finite/opt.h>
#include <state/class_finite_state.h>
#include <math/class_svd_wrapper.h>
#include <simulation/nmspc_settings.h>

void tools::finite::mps::normalize(class_finite_state & state){

    tools::log->trace("Normalizing state");
    using namespace Textra;
    using Scalar = class_finite_state::Scalar;
    state.unset_measurements();
//    std::cout << "Norm              (before normalization): " << tools::finite::measure::norm(state)  << std::endl;
//    std::cout << "Spin component sx (before normalization): " << tools::finite::measure::spin_component(state, qm::spinOneHalf::sx)  << std::endl;
//    std::cout << "Spin component sy (before normalization): " << tools::finite::measure::spin_component(state, qm::spinOneHalf::sy)  << std::endl;
//    std::cout << "Spin component sz (before normalization): " << tools::finite::measure::spin_component(state, qm::spinOneHalf::sz)  << std::endl;
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
    size_t end_traversals = 1;
    size_t max_traversals = 10;
    size_t step = 0;
    int direction = 1;
    Eigen::Tensor<Scalar,3> U;
    Eigen::Tensor<Scalar,1> S;
    Eigen::Tensor<Scalar,3> V;
    double norm;
    bool svd_success = true;

    while(num_traversals <= end_traversals){
        Eigen::Tensor<Scalar,4> theta = state.get_theta(pos_A);
        try {std::tie(U,S,V,norm) = svd.schmidt_with_norm(theta,state.get_chi_max()); svd_success = true;}
        catch(std::exception &ex){
            tools::log->error("Skipping normalization step.\n\t SVD failed at positions A:{} C:{} B:{} , step {}:\n\t{}", pos_A, pos_LC, pos_B, step, ex.what());
            svd_success = false;
            end_traversals = std::min(end_traversals+1, max_traversals);
        }
        if(svd_success){
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
//    std::cout << "Norm              (after normalization): " << tools::finite::measure::norm(state)  << std::endl;
//    std::cout << "Spin component sx (after normalization): " << tools::finite::measure::spin_component(state, qm::spinOneHalf::sx)  << std::endl;
//    std::cout << "Spin component sy (after normalization): " << tools::finite::measure::spin_component(state, qm::spinOneHalf::sy)  << std::endl;
//    std::cout << "Spin component sz (after normalization): " << tools::finite::measure::spin_component(state, qm::spinOneHalf::sz)  << std::endl;


}


void tools::finite::opt::truncate_theta(Eigen::Tensor<Scalar,3> &theta, class_finite_state & state, long chi_, double SVDThreshold){
    state.unset_measurements();
    if (state.active_sites.empty()) throw std::runtime_error("truncate_theta: No active sites to truncate");
    if (theta.size() == 0)          throw std::runtime_error("truncate_theta: Theta is empty");
    auto theta_map = Eigen::Map<Eigen::VectorXcd>(theta.data(),theta.size());
    auto fullnorm  = tools::finite::measure::norm(state);
    auto thetanorm = theta_map.norm();
    if(std::abs(fullnorm  - 1.0) > 1e-10){ tools::log->error("Norm before truncation too far from unity: {:.16f}",fullnorm);}
    if(std::abs(thetanorm - 1.0) > 1e-10){ tools::log->error("Norm of theta too far from unity: {:.16f}",thetanorm); theta_map.normalize();}



    if (state.get_direction() == 1){
        tools::finite::opt::truncate_right(theta,state,chi_,SVDThreshold);
    }else{
        tools::finite::opt::truncate_left(theta,state,chi_,SVDThreshold);
    }
    state.unset_measurements();
    fullnorm  = tools::finite::measure::norm(state);
    if(std::abs(fullnorm  - 1.0) > 1e-10){
        tools::log->error("Norm after truncation too far from unity: {:.16f}",fullnorm);
        tools::finite::mps::normalize(state);
    }


}



void tools::finite::opt::truncate_right(Eigen::Tensor<std::complex<double>,3> &theta, class_finite_state & state, long chi_, double SVDThreshold){
    tools::log->trace("Truncating multitheta from left to right");
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);
    using Scalar = class_finite_state::Scalar;
    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
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
        tools::log->trace("Truncation error site {} = {}, chi = {}", site,SVD.get_truncation_error(),S.dimension(0));
        Eigen::Tensor<Scalar,3> L_U = Textra::asDiagonalInversed(state.get_L(site)).contract(U,Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});

        state.get_G(site)   = L_U;
        state.get_L(site+1) = S;

        if(active_sites.size() > 2){
            Eigen::Tensor<Scalar,3> temp = Textra::asDiagonal(S).contract(V, Textra::idx({1},{1})).shuffle(Textra::array3{1,0,2});
            V = temp;
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


}


void tools::finite::opt::truncate_left(Eigen::Tensor<std::complex<double>,3> &theta, class_finite_state & state, long chi_, double SVDThreshold){
    tools::log->trace("Truncating multitheta from right to left");
    class_SVD SVD;
    SVD.setThreshold(SVDThreshold);
    using Scalar = class_finite_state::Scalar;
    using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
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
        tools::log->trace("Truncation error site {} = {}, chi = {}", site,SVD.get_truncation_error(),S.dimension(0));
        Eigen::Tensor<Scalar,3> V_L = V.contract(Textra::asDiagonalInversed(state.get_L(site+1)), Textra::idx({2},{0}));

        state.get_G(site) = V_L;
        state.get_L(site) = S;

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

    if (not Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1 >>(L_U.data(),L_U.size()).allFinite() )
        throw std::runtime_error("L_U has nan's or inf's");

    Eigen::Tensor<Scalar,2> leftID = state.get_MPS(site).get_A()
            .contract(state.get_MPS(site).get_A().conjugate(), Textra::idx({0,1},{0,1}) );
    auto leftIDmap = Textra::Tensor2_to_Matrix(leftID);
    if(not leftIDmap.isIdentity(1e-14)) throw std::runtime_error(fmt::format("Not left normalized at site {}", site));

}


void tools::finite::opt::truncate_theta(Eigen::Tensor<std::complex<double>,4> &theta, class_finite_state & state, long chi_, double SVDThreshold) {
    tools::log->trace("Truncating theta");
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

