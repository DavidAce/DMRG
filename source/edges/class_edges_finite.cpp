//
// Created by david on 2020-05-12.
//

#include "class_edges_finite.h"
#include "class_env_ene.h"
#include "class_env_var.h"
#include <tools/common/log.h>
#include <math/nmspc_math.h>

// template<typename env_type>
// size_t class_edges_finite::env_pair<env_type>::get_position()const {
//    if(L.get_position() != R.get_position()) throw std::runtime_error(fmt::format("Position mismatch of edge pair: {} != {}",L.get_position() !=
//    R.get_position() )); return L.get_position();
//}

template<typename env_type>
void class_edges_finite::env_pair<env_type>::assertValidity() const {
    L.assertValidity();
    R.assertValidity();
}
// Explicit instantiations
template struct class_edges_finite::env_pair<class_env_ene>;
template struct class_edges_finite::env_pair<class_env_var>;

size_t class_edges_finite::get_length() const {
    if(not math::all_equal(env_ene.size(), env_var.size()))
        throw std::runtime_error(fmt::format("Size mismatch in environments: ene {} != var {}", env_ene.size(), env_var.size()));
    return env_ene.size();
}
//
//size_t class_edges_finite::get_position() const {
//    // We scan along the environment list. While
//    size_t position = 0;
//    for(size_t pos = 0; pos < get_length(); pos++) {
//        auto pos_env_ene = get_env_ene(pos).get_position();
//        auto pos_env_var = get_env_var(pos).get_position();
//        if(not all_equal(pos_env_ene, pos_env_var))
//            throw std::runtime_error(fmt::format("Position mismatch in environments: ene {} != var {}", pos_env_ene, pos_env_var));
//        if(not all_equal(get_env_ene(pos).side, get_env_var(pos).side))
//            throw std::runtime_error(fmt::format("Side mismatch in environments: ene {} != var {}", get_env_ene(pos).side, get_env_var(pos).side));
//
//        if(get_env_ene(pos).side == "L")
//            position = pos_env_ene;
//        else
//            break;
//    }
//    for(size_t pos = position + 1; position < get_length(); pos++) {
//        if(get_env_ene(pos).side == "L") throw std::logic_error(fmt::format("Left side environment env_ene found on right side at position {}", pos));
//        if(get_env_var(pos).side == "L") throw std::logic_error(fmt::format("Left side environment env_var found on right side at position {}", pos));
//    }
//    return position;
//}

bool class_edges_finite::isReal() const {
    for(const auto &env : env_ene)
        if(not env.L.isReal() or not env.R.isReal()) return false;
    for(const auto &env : env_var)
        if(not env.L.isReal() or not env.R.isReal()) return false;
    return true;
}

bool class_edges_finite::hasNaN() const {
    for(const auto &env : env_ene)
        if(env.L.hasNaN() or env.R.hasNaN()) return true;
    for(const auto &env : env_var)
        if(env.L.hasNaN() or env.R.hasNaN()) return true;
    return false;
}



void class_edges_finite::assertValidity() const {
    for(auto &env : env_ene) env.assertValidity();
    for(auto &env : env_var) env.assertValidity();
//    if(get_position() >= get_length()) throw std::runtime_error("Position out of bounds");
}

class_edges_finite::env_pair_ref<class_env_ene> class_edges_finite::get_env_ene(size_t posL, size_t posR) const {
    if(posL > posR) throw std::range_error(fmt::format("get_env_ene(posL,posR): posL is out of range posL {} > posR {}", posL, posR));
    if(posR >= get_length()) throw std::range_error(fmt::format("get_env_ene(posL,posR): posR is out of range {} | system size {}", posR, get_length()));
    const auto &L = std::next(env_ene.begin(), static_cast<long>(posL))->L;
    const auto &R = std::next(env_ene.begin(), static_cast<long>(posR))->R;
    return {L, R};
}

class_edges_finite::env_pair_ref<class_env_var> class_edges_finite::get_env_var(size_t posL, size_t posR) const {
    if(posL > posR) throw std::range_error(fmt::format("get_env_var(posL,posR): posL is out of range posL {} > posR {}", posL, posR));
    if(posR >= get_length()) throw std::range_error(fmt::format("get_env_var(posL,posR): posR is out of range {} | system size {}", posR, get_length()));
    const auto &L = std::next(env_var.begin(), static_cast<long>(posL))->L;
    const auto &R = std::next(env_var.begin(), static_cast<long>(posR))->R;
    return {L, R};
}

class_edges_finite::env_pair_ref<class_env_ene> class_edges_finite::get_env_ene(std::optional<std::list<size_t>> sites) const {
    if(not sites) sites = active_sites;
    if(sites.value().empty()) throw std::runtime_error("Could not get edges: active site list is empty");
    return get_env_ene(sites.value().front(), sites.value().back());
}

class_edges_finite::env_pair_ref<class_env_var> class_edges_finite::get_env_var(std::optional<std::list<size_t>> sites) const {
    if(not sites) sites = active_sites;
    if(sites.value().empty()) throw std::runtime_error("Could not get edges: active site list is empty");
    return get_env_var(sites.value().front(), sites.value().back());
}


std::pair<class_edges_finite::ene_blk_ref,class_edges_finite::ene_blk_ref> class_edges_finite::get_multisite_ene(std::optional <std::list<size_t>> sites) const{
    const auto & envs = get_env_ene(std::move(sites));
    return std::make_pair(envs.L.get().block,envs.R.get().block);
}

std::pair<class_edges_finite::var_blk_ref,class_edges_finite::var_blk_ref> class_edges_finite::get_multisite_var(std::optional <std::list<size_t>> sites) const{
    const auto & envs = get_env_var(std::move(sites));
    return std::make_pair(envs.L.get().block,envs.R.get().block);
}
