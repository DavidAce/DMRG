//
// Created by david on 2020-05-12.
//

#include "class_env_var.h"
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>

#include <utility>

class_env_var::class_env_var(std::string side_, const class_mps_site &MPS, const class_mpo_site &MPO) : class_env_base(std::move(side_), MPS, MPO) {
    set_edge_dims(MPS, MPO);
}

class_env_var class_env_var::enlarge(const class_mps_site &MPS, const class_mpo_site &MPO) {
    if(MPS.get_position() != MPO.get_position()) throw std::logic_error("MPS and MPO not at the same position!");
    class_env_var env = *this;
    if(env.sites == 0 and not env.edge_has_been_set) {
        env.set_edge_dims(MPS, MPO);
        env.position = MPS.get_position();
        return env;
    }
    env.enlarge(MPS.get_M_bare(), MPO.MPO2().real().cast<Scalar>());
    if(env.side == "L") {
        env.position = MPS.get_position() + 1;
    } else if(env.side == "R") {
        env.position = MPS.get_position() - 1;
    } else {
        throw std::logic_error("Expected environment side L or R, got: " + side);
    }

    return env;
}


void class_env_var::set_edge_dims(const class_mps_site &MPS, const class_mpo_site &MPO) {
    if(edge_has_been_set) return;
    if(side == "L")
        set_edge_dims(MPS.get_M_bare(),MPO.MPO2(),MPO.get_MPO2_edge_left());
    else if (side == "R")
        set_edge_dims(MPS.get_M_bare(),MPO.MPO2(),MPO.get_MPO2_edge_right());
    else throw std::runtime_error(fmt::format("Wrong side: {}", side));
}
