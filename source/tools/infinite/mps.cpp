//
// Created by david on 2019-06-25.
//

#include <simulation/nmspc_settings.h>
#include <state/class_mps_2site.h>
#include <state/class_mps_site.h>
#include <state/class_state_infinite.h>
#include <tools/common/log.h>
#include <tools/infinite/mps.h>

using Scalar = std::complex<double>;
void tools::infinite::mps::initialize(class_state_infinite &state, const std::string &model_type_str) {
    tools::log->trace("Constructing 2site MPS");
    size_t spin_dim = 2;
    if(model_type_str == "ising_tf_rf_nn")
        spin_dim = settings::model::ising_tf_rf_nn::spin_dim;
    else if(model_type_str == "ising_selfdual_tf_rf_nn")
        spin_dim = settings::model::ising_selfdual_tf_rf_nn::spin_dim;
    else
        throw std::runtime_error(fmt::format("Unexpected model string: {}", model_type_str));
    Eigen::Tensor<Scalar, 3> M(static_cast<long>(spin_dim), 1, 1);
    Eigen::Tensor<Scalar, 1> L(1);

    // Default is a product state, spins pointing up in z.
    M.setZero();
    M(0, 0, 0) = 1;
    L.setConstant(1.0);
    state.MPS->MPS_A = std::make_unique<class_mps_site>(M, L, 0);
    state.MPS->MPS_B = std::make_unique<class_mps_site>(M, L, 1);
    state.MPS->MPS_A->set_LC(L);
}

// void tools::infinite::mps::rebuild_environments(class_state_infinite & state){
//
//}

class_state_infinite tools::infinite::mps::set_random_state(const class_state_infinite &state, [[maybe_unused]] const std::string &parity) {
    throw std::runtime_error("You need to implement set random state for infinite state");
    return state;
}
