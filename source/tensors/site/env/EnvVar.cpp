#include "EnvVar.h"
#include <config/debug.h>
#include <math/num.h>
#include <tensors/site/mpo/MpoSite.h>
#include <tensors/site/mps/MpsSite.h>
#include <utility>

EnvVar::EnvVar(std::string side_, const MpsSite &mps, const MpoSite &mpo) : EnvBase(std::move(side_), "var", mps, mpo) { set_edge_dims(mps, mpo); }

EnvVar EnvVar::enlarge(const MpsSite &mps, const MpoSite &mpo) const {
    tools::log->trace("class_env_{}::enlarge(mps,mpo): side({}), pos({})", tag, side, get_position());
    // enlarge() uses "this" block together with mps and mpo to generate a new environment block corresponding to a neighboring site
    if constexpr(settings::debug)
        if(not num::all_equal(get_position(), mps.get_position(), mpo.get_position()))
            throw std::logic_error(fmt::format("class_env_{}::enlarge(): side({}), pos({}): All positions are not equal: env {} | mps {} | mpo {}", tag, side,
                                               get_position(), get_position(), mps.get_position(), mpo.get_position()));

    EnvVar env = *this;

    if(env.sites == 0 and (not block or block->size() == 0)) {
        env.set_edge_dims(mps, mpo);
        env.position = mps.get_position();
        return env;
    }

    env.enlarge(mps.get_M_bare(), mpo.MPO2());
    // Update positions assuming this is a finite chain.
    // This needs to be corrected (on the right side) on infinite chains
    if(env.side == "L")
        env.position = mps.get_position() + 1;
    else if(env.side == "R")
        env.position = mps.get_position() - 1;
    else
        throw std::logic_error("Expected environment side L or R, got: " + side);

    env.tag = "var";
    // Save the hash id's used to create the new block in env
    env.unique_id_env = get_unique_id();
    env.unique_id_mps = mps.get_unique_id();
    env.unique_id_mpo = mpo.get_unique_id_sq();
    if constexpr(settings::debug) {
        tools::log->trace("class_env_{}::enlarge(mps,mpo): side({}), pos({}): unique_id_env: {}", tag, side, get_position(), env.unique_id_env.value());
        tools::log->trace("class_env_{}::enlarge(mps,mpo): side({}), pos({}): unique_id_mps: {}", tag, side, get_position(), env.unique_id_mps.value());
        tools::log->trace("class_env_{}::enlarge(mps,mpo): side({}), pos({}): unique_id_mpo: {}", tag, side, get_position(), env.unique_id_mpo.value());
    }
    return env;
}

void EnvVar::refresh(const EnvVar &env, const MpsSite &mps, const MpoSite &mpo) {
    // If side == L, env,mps and mpo are all corresponding to the neighbor on the left
    // If side == R, env,mps and mpo are all corresponding to the neighbor on the right
    if constexpr(settings::debug)
        if(not num::all_equal(env.get_position(), mps.get_position(), mpo.get_position()))
            throw std::logic_error(fmt::format("class_env_{}::enlarge(): side({}), pos({}),: All positions are not equal: env {} | mps {} | mpo {}", tag, side,
                                               get_position(), get_position(), mps.get_position(), mpo.get_position()));

    if(side == "L" and get_position() != mps.get_position() + 1)
        throw std::logic_error(
            fmt::format("EnvVar::refresh(pos == {}): This env{} needs env, mps and mpo at position {}", get_position(), side, get_position() - 1));
    if(side == "R" and get_position() + 1 != mps.get_position())
        throw std::logic_error(
            fmt::format("EnvVar::refresh(pos == {}): This env{} needs env, mps and mpo at position {}", get_position(), side, get_position() + 1));

    // We refresh this block if any of these conditions hold:
    //   not has_block()
    //   unique_id_env != env.unique_id;
    //   unique_id_mps != mps.unique_id;
    //   unique_id_mpo != mpo.unique_id_sq;

    if(not has_block()) {
        if constexpr(settings::debug) tools::log->trace("Refreshing {} env{}({}): missing block", tag, side, get_position());
        *this = env.enlarge(mps, mpo);
        return;
    }
    bool        refresh = false;
    std::string reason;
    if(env.get_unique_id() != unique_id_env) {
        refresh = true;
        if constexpr(settings::debug) {
            reason.append(fmt::format("| env({}) new {} ", env.get_position(), env.get_unique_id()));
            if(unique_id_env) reason.append(fmt::format("!= old {} ", unique_id_env.value()));
        }
    }
    if(mps.get_unique_id() != unique_id_mps) {
        refresh = true;
        if constexpr(settings::debug) {
            reason.append(fmt::format("| mps({}) new {} ", mps.get_position(), mps.get_unique_id()));
            if(unique_id_mps) reason.append(fmt::format("!= old {} ", unique_id_mps.value()));
        }
    }
    auto mpo_unique_id = tag == "ene" ? mpo.get_unique_id() : mpo.get_unique_id_sq();
    if(mpo_unique_id != unique_id_mpo) {
        refresh = true;
        if constexpr(settings::debug) {
            reason.append(fmt::format("| mpo({}) new {} ", mpo.get_position(), mpo_unique_id));
            if(unique_id_mpo) reason.append(fmt::format("!= old {} ", unique_id_mpo.value()));
        }
    }

    if(refresh) {
        [[maybe_unused]] size_t unique_id_bef;
        if constexpr(settings::debug) {
            unique_id_bef = get_unique_id();
            tools::log->trace("Refreshing {} env{}({}): modified {}", tag, side, get_position(), reason);
        }

        build_block(*env.block, mps.get_M_bare(), mpo.MPO2());
        unique_id_env = env.get_unique_id();
        unique_id_mps = mps.get_unique_id();
        unique_id_mpo = mpo.get_unique_id_sq();

        if constexpr(settings::debug) {
            if(unique_id_bef == get_unique_id()) tools::log->warn("Refreshing {} env{}({}): id did not change: {}", tag, side, get_position(), unique_id_bef);
            //                throw std::logic_error(fmt::format("Refreshing {} env{}({}): failed: id did not change: {}", tag, side, get_position(),
            //                unique_id_bef));
        }
    }
}

void EnvVar::set_edge_dims(const MpsSite &MPS, const MpoSite &MPO) {
    if(block and block->size() != 0) return;
    if constexpr(settings::debug) tools::log->trace("Setting edge dims on env{}({}) {}", side, get_position(), tag);
    if(side == "L")
        set_edge_dims(MPS.get_M_bare(), MPO.MPO2(), MPO.get_MPO2_edge_left());
    else if(side == "R")
        set_edge_dims(MPS.get_M_bare(), MPO.MPO2(), MPO.get_MPO2_edge_right());
    else
        throw std::runtime_error(fmt::format("Wrong side: {}", side));
}
