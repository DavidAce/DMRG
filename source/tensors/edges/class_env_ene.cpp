//
// Created by david on 2020-05-12.
//

#include "class_env_ene.h"
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>

#include <utility>

class_env_ene::class_env_ene(std::string side_, const class_mps_site &MPS, const class_mpo_site &MPO) : class_env_base(std::move(side_), MPS, MPO) {
    tag = "ene";
    set_edge_dims(MPS, MPO);
}

class_env_ene class_env_ene::enlarge(const class_mps_site &MPS, const class_mpo_site &MPO) {
    if(MPS.get_position() != MPO.get_position())
        throw std::logic_error(fmt::format("MPS and MPO have different positions: {} != {}", MPS.get_position(), MPO.get_position()));

    if(not edge_has_been_set) throw std::logic_error("Have to set edge dimensions first!");

    class_env_ene env = *this;

    env.enlarge(MPS.get_M_bare(), MPO.MPO());
    // Update positions assuming this is a finite chain.
    // This needs to be corrected (on the right side) on infinite chains
    if(env.side == "L") {
        env.position = MPS.get_position() + 1;
    } else if(env.side == "R") {
        env.position = MPS.get_position() - 1;
    } else {
        throw std::logic_error("Expected environment side L or R, got: " + side);
    }

    return env;
}

void class_env_ene::set_edge_dims(const class_mps_site &MPS, const class_mpo_site &MPO) {
    if(edge_has_been_set) return;
    if(side == "L")
        set_edge_dims(MPS.get_M_bare(),MPO.MPO(),MPO.get_MPO_edge_left());
    else if (side == "R")
        set_edge_dims(MPS.get_M_bare(),MPO.MPO(),MPO.get_MPO_edge_right());
    else throw std::runtime_error(fmt::format("Wrong side: {}", side));

}


