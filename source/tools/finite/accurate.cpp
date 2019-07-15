//
// Created by david on 2019-07-07.
//
#include <tools/nmspc_tools.h>
#include <state/class_finite_state.h>
#include <general/nmspc_tensor_extra.h>

using namespace Textra;

double tools::finite::measure::accurate::energy(const class_finite_state &state){
    class_finite_state state_tmp = state;
    auto mpo_list = tools::finite::ops::make_mpo_list(state.MPO_L, state.MPO_R);
    auto &ENV_L = state.ENV_L.front().block;
    auto &ENV_R = state.ENV_R.back().block;
    tools::finite::ops::apply_mpos(state_tmp, mpo_list,ENV_L,ENV_R);
    return tools::finite::ops::overlap(state_tmp,state);
}


double tools::finite::measure::accurate::energy_per_site(const class_finite_state &state){
    return tools::finite::measure::accurate::energy(state)/state.get_length();
}


double tools::finite::measure::accurate::energy_variance(const class_finite_state &state){
    double energy_per_site = tools::finite::measure::accurate::energy_per_site(state);
    auto mpo_list = tools::finite::ops::make_mpo_list(state.MPO_L, state.MPO_R);
    size_t pos = 0;
    for (auto & mpo: mpo_list){
        mpo = state.get_MPO(pos++).MPO_reduced_view(energy_per_site);
    }

    auto &ENV_L = state.ENV_L.front().block;
    auto &ENV_R = state.ENV_R.back().block;
    class_finite_state state_tmp = state;
    tools::finite::ops::apply_mpos(state_tmp, mpo_list,ENV_L,ENV_R);
    return tools::finite::ops::overlap(state_tmp,state_tmp);
}


double tools::finite::measure::accurate::energy_variance_per_site(const class_finite_state &state){
    return tools::finite::measure::accurate::energy_variance(state)/state.get_length();
}

