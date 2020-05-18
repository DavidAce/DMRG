//
// Created by david on 2019-08-12.
//
#include "mps.h"
#include <bitset>
#include <general/nmspc_quantum_mechanics.h>
#include <math/nmspc_random.h>
#include <config/nmspc_settings.h>
#include <tensors/state/class_environment.h>
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/measure.h>

using Scalar = tools::finite::mps::Scalar;
int get_sign(const std::string & parity_sector){
    if      (parity_sector.at(0) == '+') return 1;
    else if (parity_sector.at(0) == '-') return -1;
    else return 0;
}

std::string get_axis(const std::string & parity_sector){
    int sign = get_sign(parity_sector);
    if (sign == 0){
        return parity_sector.substr(0,1);
    }
    else{
        return parity_sector.substr(1,1);
    }
}

int get_elem(const std::string & parity_sector){
    int sign = get_sign(parity_sector);
    if (sign == 1 ) return 0;
    if (sign == -1) return 1;
    if (sign == 0 ) return rn::uniform_integer_01();
    throw std::runtime_error("Invalid sign in get_elem.");
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


Eigen::Vector2cd get_eigvec(const std::string & parity, const int sign){

    if (parity == "x" and sign >  0 ) return qm::spinOneHalf::sx_eigvecs[0];
    if (parity == "x" and sign <= 0 ) return qm::spinOneHalf::sx_eigvecs[1];
    if (parity == "y" and sign >  0 ) return qm::spinOneHalf::sy_eigvecs[0];
    if (parity == "y" and sign <= 0 ) return qm::spinOneHalf::sy_eigvecs[1];
    if (parity == "z" and sign >  0 ) return qm::spinOneHalf::sz_eigvecs[0];
    if (parity == "z" and sign <= 0 ) return qm::spinOneHalf::sz_eigvecs[1];
    throw std::runtime_error(fmt::format("get_eigvec given invalid parity sector: {} in {}", parity,sign));
}

void tools::finite::mps::internals::set_product_state_in_parity_sector_from_bitset(class_state_finite & state, const std::string &parity_sector, const long state_number){
    tools::log->trace("Setting product state from bitset of number {} in sector {}", state_number,parity_sector);
    if (state_number < 0){
        throw std::runtime_error(fmt::format("Can't set sector from bitset with negative state number: {}", state_number));
    }
    std::vector<std::string> ok_parity_sectors = {"x","+x","-x","y","+y","-y", "z","+z","-z"};
    bool parity_sector_is_defined = std::find(ok_parity_sectors.begin(), ok_parity_sectors.end(), parity_sector) != ok_parity_sectors.end();
    if (not parity_sector_is_defined)
        throw std::logic_error(fmt::format("Can't use bitfield of state_number {} to set product state in sector [{}]. Choose one of (+-) x,y or z.",state_number, parity_sector));


    constexpr long maxbits = 128;
    if (maxbits < state.get_length()) throw std::range_error("Max supported state length for bitset is 128");
    std::bitset<maxbits> bs (state_number);
    std::vector<int>  bs_vec;
    std::vector<std::string> ud_vec;
    for(size_t i = 0; i < state.get_length() ; i++ ) bs_vec.emplace_back(bs[i]) ;
    std::string axis = get_axis(parity_sector);
    int sector       = get_sign(parity_sector);
    Eigen::Tensor<Scalar,1> L (1);
    L.setConstant(1.0);
    tools::log->info("Initializing state from bitfield of state number {} with eigvecs of σ{} in sector {}: {}", state_number, axis,sector, bs_vec);
    int carry_sign = 1;
    for (auto &mps : state.mps_sites){
        int sign = 2*bs[mps.get_position()] - 1;
        carry_sign *= sign;
        mps.set_mps(Textra::MatrixTensorMap(get_eigvec(axis,sign).normalized(),2,1,1), L);
        std::string arrow = sign < 0 ? "↓" : "↑";
        ud_vec.emplace_back(arrow);
        if(mps.isCenter())
            mps.set_LC(L);
    }

    if(sector * carry_sign == -1){
        //Flip the last spin to get the correct total sign.
        auto &mps = state.mps_sites.back();
        int sign = 2*bs[mps.get_position()] - 1;
        sign *= -1;
        mps.set_mps(Textra::MatrixTensorMap(get_eigvec(axis,sign).normalized(),2,1,1), L);
        std::string arrow = sign < 0 ? "↓" : "↑";
        ud_vec.back() = arrow;
    }
    tools::log->info("Initialized  state from bitfield of state number {} with eigvecs of σ{} in sector {}: {}", state_number, axis,carry_sign, ud_vec);

}



void tools::finite::mps::internals::set_product_state_in_parity_sector_randomly(class_state_finite & state, const std::string &parity_sector){
    tools::log->debug("Setting product state randomly in parity sector {}", parity_sector);
    state.clear_measurements();
    state.clear_cache();
    Eigen::Tensor<Scalar,1> L (1);
    std::string axis = get_axis(parity_sector);
    int sector       = get_sign(parity_sector);
    int sign  = 1;

    L.setConstant(1.0);
    for (auto &mps : state.mps_sites){
        sign = 2* rn::uniform_integer_01() - 1;
        mps.set_mps(Textra::MatrixTensorMap(get_eigvec(axis,sign).normalized(), 2, 1, 1), L);
        if(mps.isCenter())
            state.mps_sites.back().set_LC(L);
    }
    auto spin_component   = tools::finite::measure::spin_component(state,axis);
    if(axis == "x" and spin_component * sector < 0 ){
        state.mps_sites.back().set_mps(Textra::MatrixTensorMap(get_eigvec(axis,-sign).normalized(), 2, 1, 1), L);
        spin_component   = tools::finite::measure::spin_component(state,axis);
    }
    state.clear_measurements();
    log->info("Finished reset to product state with global spin component {} = {}",axis, spin_component);
    if (spin_component * sector < 0) throw std::logic_error("Could not initialize_state in the correct parity sector");

}



void tools::finite::mps::internals::set_product_state_randomly(class_state_finite & state, const std::string &parity_sector, bool use_pauli_eigenstates){
    std::vector<std::string> ok_parity_sectors = {"x","+x","-x","y","+y","-y", "z","+z","-z"};
    bool parity_sector_is_defined = std::find(ok_parity_sectors.begin(), ok_parity_sectors.end(), parity_sector) != ok_parity_sectors.end();
    if (parity_sector_is_defined and use_pauli_eigenstates){
        // Case a)
        set_product_state_in_parity_sector_randomly(state,parity_sector);
    }
    else if (parity_sector_is_defined and not use_pauli_eigenstates){
        set_product_state_randomly(state,"random",false);
        state = tools::finite::ops::get_projection_to_closest_parity_sector(state,parity_sector);
    }
    else if (parity_sector == "randomAxis") {
        std::vector<std::string> possibilities = {"x", "y", "z"};
        std::string chosen_axis = possibilities[static_cast<unsigned long>(rn::uniform_integer_box(0, 2))];
        set_product_state_in_parity_sector_randomly(state, chosen_axis);
    }else if (parity_sector == "random") {
        Eigen::Tensor<Scalar,1> L (1);
        L.setConstant(1.0);
        for (auto &mps : state.mps_sites){
            mps.set_mps(Textra::MatrixTensorMap(Eigen::VectorXcd::Random(2).normalized(), 2, 1, 1), L);
            if(mps.isCenter())
                mps.set_LC(L);
        }
    }else if (parity_sector == "none"){
        return;
    }else{
        throw std::runtime_error(fmt::format(R"(Wrong pauli string. Expected one of (+-) "x","y","z", "randomAxis", "random" or "none". Got: )" + parity_sector));
    }

}
