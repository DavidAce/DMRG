//
// Created by david on 2019-06-24.
//


#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_finite_chain_state.h>


void mpstools::finite::multisite::compute_best_jump(class_finite_chain_state &state){
    using namespace Textra;
    size_t position  = state.get_position();
    size_t direction = state.get_direction();
    size_t length    = state.get_length();
    std::vector<size_t> costs;
    std::vector<Eigen::DSizes<long,3>> dims;

    while(position >= 0 and position < length){
        dims.emplace_back(state.get_G(position).dimensions());
        size_t cost = direction > 1 ? dims.back()[1]*dims.front()[2] : dims.front()[1]*dims.back()[2]  ;
        for (auto &d : dims ){cost *= d[0];}
        costs.push_back(cost);
        position += direction;
    }

//    auto mapcosts = Eigen::Map<Eigen::Array<size_t, Eigen::Dynamic,1>>(costs.data(),costs.size());
    std::cout << "JUMP COSTS \n" << costs << std::endl;






}
