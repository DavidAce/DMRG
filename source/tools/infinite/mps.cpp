//
// Created by david on 2019-06-25.
//

#include <tools/nmspc_tools.h>
#include <state/class_infinite_state.h>
#include <state/class_mps_2site.h>
#include <state/class_mps_site.h>


void tools::infinite::mps::initialize  (class_infinite_state & state,std::string model_type_str){
    tools::log->trace("Constructing 2site MPS");
    long spin_dimension = 2;
    if      (model_type_str == "tf_ising")             spin_dimension = 2;
    else if (model_type_str == "tf_nn_ising")          spin_dimension = 2;
    else if (model_type_str == "selfdual_tf_rf_ising") spin_dimension = 2;

    Eigen::Tensor<Scalar,3> M(spin_dimension,1,1);
    Eigen::Tensor<Scalar,1> L(1);

    // Default is a product state, spins pointing up in z.
    M.setZero();
    M(0,0,0) = 1;
    L.setConstant(1.0);
    state.MPS->MPS_A = std::make_unique<class_mps_site>(M,L,0);
    state.MPS->MPS_B = std::make_unique<class_mps_site>(M,L,1);
    state.MPS->MPS_A->set_LC(L);
}

//void tools::infinite::mps::rebuild_environments(class_infinite_state & state){
//
//}




class_infinite_state tools::infinite::mps::set_random_state(const class_infinite_state & state, [[maybe_unused]] std::string parity,  [[maybe_unused]]  int seed_state){
    throw std::runtime_error("You need to implement set random state for infinite state");
    return state;
}
