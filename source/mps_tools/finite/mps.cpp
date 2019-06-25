//
// Created by david on 2019-01-29.
//


#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_state/class_mps_2site.h>
#include <mps_state/class_environment.h>
#include <general/nmspc_random_numbers.h>


using Scalar = std::complex<double>;

void mpstools::finite::mps::initialize_mps(class_finite_chain_state &state, std::string parity, const size_t length){
    state.clear();


    //Generate MPS
    auto spin_dim = state.MPO_L.back()->get_spin_dimension();
    Eigen::Tensor<Scalar,3> G(spin_dim,1,1);
    Eigen::Tensor<Scalar,1> L(1);
    G.setValues({{{1.0}},{{0.0}}});
    L.setValues({1});

    state.MPS_C = L;
    while(true){
        state.MPS_L.emplace_back(class_vidal_mps(G,L));
        if(state.MPS_L.size() + state.MPS_R.size() >= length){break;}
        state.MPS_R.emplace_front(class_vidal_mps(G,L));
        if(state.MPS_L.size() + state.MPS_R.size() >= length){break;}
    }


    state = mpstools::finite::ops::set_random_product_state(state,parity);
    mpstools::finite::ops::rebuild_environments(state);
    mpstools::finite::debug::check_integrity_of_mps(state);
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

