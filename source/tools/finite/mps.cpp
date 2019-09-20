//
// Created by david on 2019-01-29.
//


#include <tools/nmspc_tools.h>
#include <state/class_finite_state.h>
#include <state/class_mps_2site.h>
#include <state/class_environment.h>
#include <general/nmspc_random_numbers.h>
#include <general/nmspc_quantum_mechanics.h>
#include <simulation/nmspc_settings.h>


void tools::finite::mps::initialize(class_finite_state &state, const size_t length){
    log->info("Initializing mps");
    using Scalar = class_finite_state::Scalar;
    //Generate MPS
    Eigen::Tensor<Scalar,3> G;
    Eigen::Tensor<Scalar,1> L;
    state.MPS_C = L;
    size_t pos = 0;
    state.MPS_L.emplace_back(class_vidal_site(G,L,pos++));
    while(true){
        state.MPS_R.emplace_back(class_vidal_site(G,L,pos++));
        if(state.MPS_L.size() + state.MPS_R.size() >= length){break;}
    }
    state.truncation_error.resize(length+1);
    state.site_update_tags = std::vector<bool>(length,false);
}



void tools::finite::mps::randomize(class_finite_state &state,const std::string &parity_sector, int seed_state, bool use_pauli_eigenstates, bool enumeration)
/*!
 * There are many ways to randomize an initial product state state, based on the
 * arguments (parity_sector,seed_state,use_pauli_eigenstates, enumeration) = (string,int,true/false,true/false).
 * Let "+-sector" mean one of {"x","+x","-x","y","+y","-y", "z","+z","-z"}.

        a) ("+-sector"  ,+- ,t,f)   Set spinors to a random sequence of eigenvectors (up/down) of either
                                    sx, sy or sz pauli matrices (same pauli for all sites). If the global
                                    sign (+-) is omitted, a random sign is chosen with equal probabilities.
                                    In the x and z cases the full state will turn out to be entirely real,
                                    which improves performance.

        b) ("random"    ,+- ,f,f)   Set each spinor randomly on C2


        c) ("+-sector"  ,+- ,f,f)   Set each spinor randomly on C2 (i.e. case b) and then project the  full state
                                    to the given parity sector. If the global sign (+-) is omitted,  a random
                                    sign is chosen with equal probabilities. As a consequence of this, the
                                    full state will have always have nonzero imaginary part.

        d) ("randomAxis",+- ,f,f)   Randomly select one of {"x","y","z"} and go to case a).
        e) ("none"      ,+- ,f,f)   Does not randomize
        f) ("+-sector"  ,>=0,?,t)   Interpret seed_state as bitfield "01100010110..." and interpret these as
                                    up(0)/down(1) of either sx, sy or sz pauli matrices (same pauli for all sites)
 * Note: seed_state is only used if >= 0.
 * Note: we "use" the seed_state only once. Subsequent calls do not keep resetting the seed.
*/
{
    tools::log->trace("Randomizing mps");
    using Scalar = class_finite_state::Scalar;
    state.unset_measurements();
    state.clear_cache();
    state.tag_all_sites_have_been_updated(false);

    if(seed_state >= 0 and internals::seed_state_unused){
        rn::seed(seed_state);
        internals::seed_state_unused = false;
    }

    if   (enumeration)
         internals::set_product_state_in_parity_sector_from_bitset(state, parity_sector, seed_state);
    else internals::set_product_state_randomly(state, parity_sector, use_pauli_eigenstates);
    tools::finite::mps::rebuild_environments(state);
}


void tools::finite::mps::rebuild_environments(class_finite_state &state){
    tools::log->trace("Rebuilding environments");
    if (state.MPS_L.size() != state.MPO_L.size())
        throw std::runtime_error(fmt::format("Size mismatch in MPSL and MPOL: {} != {}",state.MPS_L.size(), state.MPO_L.size()));
    if (state.MPS_R.size() != state.MPO_R.size())
        throw std::runtime_error(fmt::format("Size mismatch in MPSR and MPOR: {} != {}",state.MPS_R.size(), state.MPO_R.size()));
    // Generate new environments
    state.ENV_L.clear();
    state.ENV_R.clear();

    state.ENV2_L.clear();
    state.ENV2_R.clear();
    size_t position = 0;
    {
        long dimMPS = state.MPS_L.front().get_chiL();
        long dimMPO = state.MPO_L.front()->MPO().dimension(0);

        auto ENV_L  = class_environment    ("L",dimMPS,dimMPO);
        auto ENV2_L = class_environment_var("L",dimMPS,dimMPO);
        ENV_L.set_position(state.MPS_L.front().get_position());
        ENV2_L.set_position(state.MPS_L.front().get_position());
        auto mpsL_it   = state.MPS_L.begin();
        auto mpoL_it   = state.MPO_L.begin();
        while(mpsL_it != state.MPS_L.end() and mpoL_it != state.MPO_L.end()) {
            state.ENV_L.emplace_back(ENV_L);
            state.ENV2_L.emplace_back(ENV2_L);
            ENV_L.enlarge(*mpsL_it, mpoL_it->get()->MPO());
            ENV2_L.enlarge(*mpsL_it, mpoL_it->get()->MPO());
            if (mpsL_it->get_position() != state.ENV_L.back().get_position())
                throw std::runtime_error(fmt::format("Size mismatch in MPSL and ENVL: {} != {}",mpsL_it->get_position(), state.ENV_L.back().get_position()));

            position ++;
            mpsL_it++;
            mpoL_it++;
        }
    }

    {
        position = state.MPS_R.back().get_position();
        long dimMPS = state.MPS_R.back().get_chiR();
        long dimMPO = state.MPO_R.back()->MPO().dimension(1);
        auto ENV_R  = class_environment    ("R",dimMPS,dimMPO);
        auto ENV2_R = class_environment_var("R",dimMPS,dimMPO);
        ENV_R.set_position(state.MPS_R.back().get_position());
        ENV2_R.set_position(state.MPS_R.back().get_position());
        auto mpsR_it   = state.MPS_R.rbegin();
        auto mpoR_it   = state.MPO_R.rbegin();
        while(mpsR_it != state.MPS_R.rend() and mpoR_it != state.MPO_R.rend()){
            state.ENV_R .emplace_front(ENV_R);
            state.ENV2_R.emplace_front(ENV2_R);
            ENV_R.enlarge(*mpsR_it, mpoR_it->get()->MPO());
            ENV2_R.enlarge(*mpsR_it, mpoR_it->get()->MPO());
            if (mpsR_it->get_position() != state.ENV_R.front().get_position())
                throw std::runtime_error(fmt::format("Size mismatch in MPSR and ENVR: {} != {}",mpsR_it->get_position(), state.ENV_R.front().get_position()));
            position --;
            mpsR_it++;
            mpoR_it++;


        }
    }
}




int tools::finite::mps::move_center_point(class_finite_state &  state){
    //Take current MPS and generate an Lblock one larger and store it in list for later loading
//    std::cout << "Current state -- Direction: " << direction << std::endl;
//    std::cout << "HA: " << state.HA->get_position() << " MPO_L back : " << MPO_L.back()->get_position() << std::endl;
//    std::cout << "HB: " << state.HB->get_position() << " MPO_R front: " << MPO_R.front()->get_position() << std::endl;
//
    auto & MPS_L  = state.MPS_L;
    auto & MPS_R  = state.MPS_R;
    auto & MPS_C  = state.MPS_C;
    auto & MPO_L  = state.MPO_L;
    auto & MPO_R  = state.MPO_R;
    auto & ENV_L  = state.ENV_L;
    auto & ENV_R  = state.ENV_R;
    auto & ENV2_L = state.ENV2_L;
    auto & ENV2_R = state.ENV2_R;
    if(ENV_L.empty()) throw std::runtime_error("ENVL is empty");
    if(ENV_R.empty()) throw std::runtime_error("ENVR is empty");
    if(MPS_L.empty()) throw std::runtime_error("MPSL is empty");
    if(MPS_R.empty()) throw std::runtime_error("MPSR is empty");
    if(MPS_L.back().get_position()  != ENV_L.back().get_position())  throw std::runtime_error("MPSL and ENVL have mismatching positions");
    if(MPS_R.front().get_position() != ENV_R.front().get_position()) throw std::runtime_error("MPSR and ENVR have mismatching positions");
    if(ENV_L.size() + ENV_R.size() != state.get_length()) throw std::runtime_error("ENVL + ENVR sizes do not add up to chain length");
    if(MPS_L.size() + MPS_R.size() != state.get_length()) throw std::runtime_error("MPSL + MPSR sizes do not add up to chain length");
    assert(ENV_L.size() + ENV_R.size() == state.get_length());
    assert(ENV_L.back().sites + ENV_R.front().sites == state.get_length() - 2);

    if (state.get_direction() == 1){
        class_environment     L  = ENV_L.back();
        class_environment_var L2 = ENV2_L.back();
        L.enlarge(MPS_L.back(), MPO_L.back()->MPO());
        L2.enlarge(MPS_L.back(), MPO_L.back()->MPO());
        ENV_L.emplace_back(L);
        ENV2_L.emplace_back(L2);

        //Note that Lblock must just have grown!!
        state.MPS_L.emplace_back(class_vidal_site(MPS_R.front().get_G(),MPS_C, MPS_R.front().get_position()));
        state.MPS_C = MPS_R.front().get_L();
        state.MPO_L.emplace_back(MPO_R.front()->clone());
        MPS_R.pop_front();
        MPO_R.pop_front();
        ENV_R.pop_front();
        ENV2_R.pop_front();
    }else{

        class_environment     R  = ENV_R.front();
        class_environment_var R2 = ENV2_R.front();
        R.enlarge(MPS_R.front(), MPO_R.front()->MPO());
        R2.enlarge(MPS_R.front(), MPO_R.front()->MPO());
        ENV_R.emplace_front(R);
        ENV2_R.emplace_front(R2);

        //Note that Rblock must just have grown!!

        state.MPS_R.emplace_front(class_vidal_site(MPS_L.back().get_G(),MPS_C, MPS_L.back().get_position()));
        state.MPS_C = MPS_L.back().get_L();
        state.MPO_R.emplace_front(MPO_L.back()->clone());

        MPS_L.pop_back();
        MPO_L.pop_back();
        ENV_L.pop_back();
        ENV2_L.pop_back();
    }

    assert(MPO_L.size() + MPO_R.size() == state.get_length());
    if(ENV_L.empty()) throw std::runtime_error("ENVL became empty");
    if(ENV_R.empty()) throw std::runtime_error("ENVR became empty");
    if(MPS_L.empty()) throw std::runtime_error("MPSL became empty");
    if(MPS_R.empty()) throw std::runtime_error("MPSR became empty");
    if(MPS_L.back().get_position()  != ENV_L.back().get_position())  throw std::runtime_error("MPSL and ENVL got mismatching positions");
    if(MPS_R.front().get_position() != ENV_R.front().get_position()) throw std::runtime_error("MPSR and ENVR got mismatching positions");
    if(ENV_L.size() + ENV_R.size() != state.get_length()) throw std::runtime_error("ENVL + ENVR sizes do not add up to chain length anymore");
    if(MPS_L.size() + MPS_R.size() != state.get_length()) throw std::runtime_error("MPSL + MPSR sizes do not add up to chain length anymore");
    //    Check edge
    if (state.position_is_any_edge()){
        state.flip_direction();
        state.increment_sweeps();
    }
    state.increment_moves();
    state.clear_cache();
    state.unset_measurements();
    state.active_sites.clear();
    return state.get_sweeps();
}



void tools::finite::mps::project_to_closest_parity_sector   (class_finite_state & state, std::string parity_sector, bool keep_bond_dimensions){
    state = tools::finite::ops::get_projection_to_closest_parity_sector(state, parity_sector, keep_bond_dimensions);
}
