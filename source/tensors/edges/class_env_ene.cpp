//
// Created by david on 2020-05-12.
//

#include "class_env_ene.h"
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <math/num.h>
#include <config/debug.h>
#include <utility>

class_env_ene::class_env_ene(std::string side_, const class_mps_site &MPS, const class_mpo_site &MPO) : class_env_base(std::move(side_), MPS, MPO) {
    tag = "ene";
    set_edge_dims(MPS, MPO);
}

class_env_ene class_env_ene::enlarge(const class_mps_site &mps, const class_mpo_site &mpo) const {
    if(mps.get_position() != mpo.get_position())
        throw std::logic_error(fmt::format("mps and mpo have different positions: {} != {}", mps.get_position(), mpo.get_position()));

    if(not edge_has_been_set) throw std::logic_error("Have to set edge dimensions first!");

    class_env_ene env = *this;

    env.enlarge(mps.get_M_bare(), mpo.MPO());
    // Update positions assuming this is a finite chain.
    // This needs to be corrected (on the right side) on infinite chains
    if(env.side == "L") {
        env.position = mps.get_position() + 1;
    } else if(env.side == "R") {
        env.position = mps.get_position() - 1;
    } else {
        throw std::logic_error("Expected environment side L or R, got: " + side);
    }
    env.tag           = "ene";
    env.unique_id_env = get_unique_id();
    env.unique_id_mps = mps.get_unique_id();
    env.unique_id_mpo = mpo.get_unique_id();
    return env;
}

void class_env_ene::refresh(const class_env_ene & env, const class_mps_site &mps, const class_mpo_site &mpo) {
    // If side == L, env,mps and mpo are all corresponding to the neighbor on the left
    // If side == R, env,mps and mpo are all corresponding to the neighbor on the right
    if(not num::all_equal(env.get_position(), mps.get_position(), mpo.get_position()))
        throw std::logic_error(fmt::format("class_env_ene::refresh(pos == {}): All positions are not equal: env {} | mps {} | mpo {}",
                                           get_position(),env.get_position(), mps.get_position(), mpo.get_position()));

    if(side == "L" and get_position() != mps.get_position()+1)
        throw std::logic_error(fmt::format("class_env_ene::refresh(pos == {}): This env{} needs env, mps and mpo at position {}", get_position(),side, get_position()-1));
    if(side == "R" and get_position()+1 != mps.get_position())
        throw std::logic_error(fmt::format("class_env_ene::refresh(pos == {}): This env{} needs env, mps and mpo at position {}", get_position(),side, get_position()+1));

    // We refresh this block if any of these conditions hold:
    //   not has_block()
    //   unique_id_env != env.unique_id;
    //   unique_id_mps != mps.unique_id;
    //   unique_id_mpo != mpo.unique_id;

    if(not has_block()){
        tools::log->trace("Refreshing {} env{}({}): missing block",tag,side,get_position());
        *this = env.enlarge(mps,mpo);
        return;
    }
    bool refresh = false;
    std::string reason;
    if(env.get_unique_id() != unique_id_env){
        refresh = true;
        if constexpr (settings::debug){
            reason.append(fmt::format("| env({}) {} !=",env.get_position(),env.get_unique_id()));
            if(unique_id_env) reason.append(fmt::format(" {} ",unique_id_env.value()));
        }
    }
    if(mps.get_unique_id() != unique_id_mps){
        refresh = true;
        if constexpr (settings::debug){
            reason.append(fmt::format("| mps({}) {} !=",mps.get_position(),mps.get_unique_id()));
            if(unique_id_mps) reason.append(fmt::format(" {} ", unique_id_mps.value()));
        }
    }
    auto mpo_unique_id = tag == "ene" ? mpo.get_unique_id() : mpo.get_unique_id_sq();
    if(mpo_unique_id != unique_id_mpo){
        refresh = true;
        if constexpr(settings::debug){
            reason.append(fmt::format("| mpo({}) {} !=",mpo.get_position(),mpo_unique_id));
            if(unique_id_mpo) reason.append(fmt::format(" {} ", unique_id_mpo.value()));
        }

    }
    if(refresh){
        if constexpr(settings::debug) tools::log->trace("Refreshing {} env{}({}): modified {}",tag,side,get_position(),reason);
        *this = env.enlarge(mps,mpo);
    }
}



void class_env_ene::set_edge_dims(const class_mps_site &MPS, const class_mpo_site &MPO) {
    if(edge_has_been_set) return;
    if(side == "L")
        set_edge_dims(MPS.get_M_bare(),MPO.MPO(),MPO.get_MPO_edge_left());
    else if (side == "R")
        set_edge_dims(MPS.get_M_bare(),MPO.MPO(),MPO.get_MPO_edge_right());
    else throw std::runtime_error(fmt::format("Wrong side: {}", side));

}


