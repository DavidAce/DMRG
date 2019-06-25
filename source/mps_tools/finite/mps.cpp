//
// Created by david on 2019-01-29.
//


#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_state/class_mps_2site.h>
#include <mps_state/class_environment.h>
#include <general/nmspc_random_numbers.h>
#include <general/class_svd_wrapper.h>
#include <general/nmspc_quantum_mechanics.h>
#include <sim_parameters/nmspc_sim_settings.h>

using Scalar = std::complex<double>;

void mpstools::finite::mps::initialize(class_finite_chain_state &state, const size_t length){
    log->info("Initializing mps");
    //Generate MPS
    Eigen::Tensor<Scalar,3> G;
    Eigen::Tensor<Scalar,1> L;
    state.MPS_C = L;
    while(true){
        state.MPS_L.emplace_back(class_vidal_mps(G,L));
        if(state.MPS_L.size() + state.MPS_R.size() >= length){break;}
        state.MPS_R.emplace_front(class_vidal_mps(G,L));
        if(state.MPS_L.size() + state.MPS_R.size() >= length){break;}
    }
}



void mpstools::finite::mps::normalize(class_finite_chain_state & state){

    mpstools::log->trace("Normalizing state");
    using namespace Textra;
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
    int direction = 1;
    Eigen::Tensor<Scalar,3> U;
    Eigen::Tensor<Scalar,1> S;
    Eigen::Tensor<Scalar,3> V;
    double norm;
    while(num_traversals < 4){
        Eigen::Tensor<Scalar,4> theta = state.get_theta(pos_A);
        bool isZero = Eigen::Map <Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).isZero();
        if(isZero){mpstools::log->warn("Theta is all zeros at positions: {},{}", pos_A, pos_B );}
        try {std::tie(U,S,V,norm) = svd.schmidt_with_norm(theta);}
        catch(std::exception &ex){
            std::cerr << "A:\n" << state.get_A(pos_A) << std::endl;
            std::cerr << "L:\n" << state.get_L(pos_LC) << std::endl;
            std::cerr << "B:\n" << state.get_B(pos_B) << std::endl;
            throw std::runtime_error("Normalization failed at step " + std::to_string(pos_A) + ": " + std::string(ex.what()) );
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


    }
    state.unset_measurements();
//    std::cout << "Norm              (after normalization): " << mpstools::finite::measure::norm(state)  << std::endl;
//    std::cout << "Spin component sx (after normalization): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sx)  << std::endl;
//    std::cout << "Spin component sy (after normalization): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sy)  << std::endl;
//    std::cout << "Spin component sz (after normalization): " << mpstools::finite::measure::spin_component(state, qm::spinOneHalf::sz)  << std::endl;


}



void mpstools::finite::mps::randomize(class_finite_chain_state &state){
    mpstools::log->trace("Randomizing mps");
    state.unset_measurements();
    Eigen::Tensor<Scalar,1> L (1);
    L.setConstant(1);

    for (auto &mpsL : state.MPS_L ){
        auto G = Textra::Matrix_to_Tensor(Eigen::VectorXcd::Random(2).normalized(),2,1,1);
//        auto G = Textra::Matrix_to_Tensor(Eigen::VectorXd::Random(2).cast<Scalar>().normalized(),2,1,1);
        mpsL.set_mps(G,L);
    }
    state.MPS_C = L;
    for (auto &mpsR : state.MPS_R ){
        auto G = Textra::Matrix_to_Tensor(Eigen::VectorXcd::Random(2).normalized(),2,1,1);
//        auto G = Textra::Matrix_to_Tensor(Eigen::VectorXd::Random(2).cast<Scalar>().normalized(),2,1,1);
        mpsR.set_mps(G,L);
    }

    mpstools::finite::mps::rebuild_environments(state);
}


void mpstools::finite::mps::rebuild_environments(class_finite_chain_state &state){
    mpstools::log->trace("Rebuilding environments");

    // Generate new environments
    state.ENV_L.clear();
    state.ENV2_L.clear();
    state.ENV_R.clear();
    state.ENV2_R.clear();

    {
        auto ENV_L  = class_environment    ("L",state.MPS_L.front().get_chiL(), state.MPO_L.front()->MPO().dimension(0));
        auto ENV2_L = class_environment_var("L",state.MPS_L.front().get_chiL(), state.MPO_L.front()->MPO().dimension(0));
        auto mpsL_it   = state.MPS_L.begin();
        auto mpoL_it   = state.MPO_L.begin();

        while(mpsL_it != state.MPS_L.end() and mpoL_it != state.MPO_L.end()) {
            state.ENV_L.push_back(ENV_L);
            state.ENV2_L.push_back(ENV2_L);
            if (mpsL_it == state.MPS_L.end()) break;
            ENV_L.enlarge(mpsL_it->get_A(), mpoL_it->get()->MPO());
            ENV2_L.enlarge(mpsL_it->get_A(), mpoL_it->get()->MPO());
//            throw_if_has_imaginary_part(state.ENV_L.back().block);

            mpsL_it++;
            mpoL_it++;
        }
    }

    {
        auto ENV_R  = class_environment("R",state.MPS_R.back().get_chiR(), state.MPO_L.back()->MPO().dimension(1));
        auto ENV2_R = class_environment_var("R",state.MPS_R.back().get_chiR(), state.MPO_L.back()->MPO().dimension(1));
        auto mpsR_it   = state.MPS_R.rbegin();
        auto mpoR_it   = state.MPO_R.rbegin();
        while(mpsR_it != state.MPS_R.rend() and mpoR_it != state.MPO_R.rend()){
            state.ENV_R.push_front(ENV_R);
            state.ENV2_R.push_front(ENV2_R);
            if (mpsR_it == state.MPS_R.rend()) break;
            ENV_R.enlarge(mpsR_it->get_B(), mpoR_it->get()->MPO());
            ENV2_R.enlarge(mpsR_it->get_B(), mpoR_it->get()->MPO());
            mpsR_it++;
            mpoR_it++;
        }
    }
    state.set_positions();
}


int mpstools::finite::mps::move_center_point(class_finite_chain_state &  state){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading
//    std::cout << "Current state -- Direction: " << direction << std::endl;
//    std::cout << "HA: " << superblock.HA->get_position() << " MPO_L back : " << MPO_L.back()->get_position() << std::endl;
//    std::cout << "HB: " << superblock.HB->get_position() << " MPO_R front: " << MPO_R.front()->get_position() << std::endl;
//
    auto & MPS_L  = state.MPS_L;
    auto & MPS_R  = state.MPS_R;
    auto & MPS_C  = state.MPS_C;
    auto & MPO_L  = state.MPO_L;
    auto & MPO_R  = state.MPO_R;
    auto & ENV_L  = state.ENV_L;
    auto & ENV_R  = state.ENV_R;
    auto & ENV2_L = state.ENV2_L;
    auto & ENV2_R = state.ENV2_R;
    assert(!MPS_L.empty() and !MPS_R.empty());
    assert(MPS_L.size() + MPS_R.size() == state.get_length());
    assert(ENV_L.size() + ENV_R.size() == state.get_length());
    assert(ENV_L.back().size + ENV_R.front().size == state.get_length() - 2);

    if (state.get_direction() == 1){
        class_environment     L  = ENV_L.back();
        class_environment_var L2 = ENV2_L.back();
        L.enlarge(MPS_L.back().get_A(), MPO_L.back()->MPO());
        L2.enlarge(MPS_L.back().get_A(), MPO_L.back()->MPO());
        ENV_L.emplace_back(L);
        ENV2_L.emplace_back(L2);

        //Note that Lblock must just have grown!!
        state.MPS_L.emplace_back(class_vidal_mps(MPS_R.front().get_G(),MPS_C));
        state.MPS_C = MPS_R.front().get_L();
        state.MPO_L.emplace_back(MPO_R.front()->clone());
        MPS_R.pop_front();
        MPO_R.pop_front();
        ENV_R.pop_front();
        ENV2_R.pop_front();
    }else{

        class_environment     R  = ENV_R.back();
        class_environment_var R2 = ENV2_R.back();
        R.enlarge(MPS_R.front().get_A(), MPO_R.front()->MPO());
        R2.enlarge(MPS_R.front().get_A(), MPO_R.front()->MPO());
        ENV_R.emplace_front(R);
        ENV2_R.emplace_front(R2);


        //Note that Rblock must just have grown!!

        state.MPS_R.emplace_front(class_vidal_mps(MPS_L.back().get_G(),MPS_C));
        state.MPS_C = MPS_L.back().get_L();
        state.MPO_L.emplace_front(MPO_L.back()->clone());

        MPS_L.pop_back();
        MPO_L.pop_back();
        ENV_L.pop_back();
        ENV2_L.pop_back();
    }

    assert(MPO_L.size() + MPO_R.size() == state.get_length());


    //    Check edge
    if (state.position_is_any_edge()){
        state.flip_direction();
        state.increment_sweeps();
    }

    state.unset_measurements();
    return state.get_sweeps();
}



void mpstools::finite::mps::project_to_closest_parity   (class_finite_chain_state & state, const std::string paulistring){
    state = mpstools::finite::ops::get_closest_parity_state(state,paulistring);
}
