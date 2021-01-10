#include "svd.h"
#include <config/nmspc_settings.h>
#include <math/svd.h>
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/split.h>
void tools::infinite::svd::truncate_theta(Eigen::Tensor<Scalar, 3> &mps, class_state_infinite &state) {
    auto mps_list = tools::common::split::split_mps(mps);
    // Presumably we get only two sites out of this
    if(mps_list.size() != 2) throw std::runtime_error("Got more than two sites from svd split");
    auto mps_it = mps_list.begin();

    for(auto &site : mps_list) { site. }
    state.set_mps(mps_list)
    //    state.MPS->truncation_error = SVD.get_truncation_error();
    //    state.MPS->MPS_A->set_LC(S);
    //    state.MPS->MPS_A->set_M(U);
    //    state.MPS->MPS_B->set_M(V);
}
