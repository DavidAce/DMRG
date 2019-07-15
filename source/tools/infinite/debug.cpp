//
// Created by david on 2019-06-24.
//
#include <iostream>
#include <iomanip>
#include <simulation/class_simulation_status.h>
#include <tools/nmspc_tools.h>
#include <state/class_infinite_state.h>
#include <general/nmspc_quantum_mechanics.h>
#include <spdlog/spdlog.h>

void tools::infinite::debug::check_integrity(const class_infinite_state & state, const class_simulation_status &sim_status)
{
    tools::log->info("Checking integrity...");
   try{
        check_integrity_of_mps(state);
    }catch(std::exception & ex){
        tools::infinite::print::print_state(state) ;
        throw std::runtime_error("Integrity check of MPS failed: " + std::string(ex.what()));
    }
    try {
        tools::log->info("Checking norms");
        auto norm_block = tools::infinite::measure::norm(state);
        if(std::abs(norm_block - 1.0) > 1e-10) {
            tools::log->warn("Norm of state far from unity: {}", norm_block);
            throw std::runtime_error("Norm of state too far from unity: " + std::to_string(norm_block));
        }
    }
    catch(std::exception & ex){
        throw std::runtime_error("Integrity check of norm failed: " + std::string(ex.what()));
    }
}




void tools::infinite::debug::check_integrity_of_mps(const class_infinite_state &state){
//    {
//        tools::log->trace("Checking integrity of MPS");
//        tools::log->trace("\tChecking system sizes");
//        if(state.MPS_L.size() + state.MPS_R.size() != state.get_length() )
//            throw std::runtime_error("Mismatch in MPS size: " + std::to_string(state.MPS_L.size() + state.MPS_R.size()) + " " + std::to_string(state.get_length()));
//
//        if(state.ENV_L.size() + state.ENV_R.size() != state.get_length())
//            throw std::runtime_error("Mismatch in ENV size: " + std::to_string(state.ENV_L.size() + state.ENV_R.size()) + " " + std::to_string(state.get_length()));
//
//        if(state.ENV2_L.size() + state.ENV2_R.size() != state.get_length())
//            throw std::runtime_error("Mismatch in ENV size: " + std::to_string(state.ENV_L.size() + state.ENV_R.size()) + " " + std::to_string(state.get_length()));
//
//        if( state.ENV_L.back().size != state.get_position())
//            throw std::runtime_error("Mismatch in ENV_L size and position: " + std::to_string(state.ENV_L.back().size) + " " + std::to_string(state.get_position()));
//
//        if(state.ENV_R.front().size != state.get_length() - state.get_position() - 2)
//            throw std::runtime_error("Mismatch in ENV_R size+1 and length-position: " + std::to_string(state.ENV_R.front().size) + " " + std::to_string(state.get_length() - state.get_position()-2));
//
//        tools::log->trace("\tChecking matrix sizes on the left side");
//        //Check left side of the state
//        auto mps_it  = state.MPS_L.begin();
//        auto mps_nx  = state.MPS_L.begin();
//        auto env_it  = state.ENV_L.begin();
//        auto env2_it = state.ENV2_L.begin();
//        auto mpo_it  = state.MPO_L.begin();
//        std::advance(mps_nx,1);
//        int i = 0;
//        while(
//                mps_it  != state.MPS_L.end() and
//                mps_nx  != state.MPS_L.end() and
//                env_it  != state.ENV_L.end() and
//                env2_it != state.ENV2_L.end() and
//                mpo_it  != state.MPO_L.end()
//                )
//        {
//            if(mps_it->get_chiR() != mps_nx->get_chiL())
//                throw std::runtime_error("Mismatch in adjacent MPS dimensions (left side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_chiR()) + " " + std::to_string(mps_nx->get_chiL()));
//
//            if(mps_it->get_position() != env_it->get_position())
//                throw std::runtime_error("Mismatch in MPS and ENV positions (left side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(env_it->get_position()));
//
//            if(mps_it->get_position() != env_it->size)
//                throw std::runtime_error("Mismatch in MPS position and ENV size (left side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(env_it->size));
//
//            if(mps_it->get_chiL() != env_it->block.dimension(0))
//                throw std::runtime_error("Mismatch in MPS and ENV dimensions (left side) @ site " + std::to_string(i) + ": " +  std::to_string(mps_it->get_chiL()) + " " + std::to_string(env_it->block.dimension(0)));
//
//            if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
//                throw std::runtime_error("Mismatch in ENV and MPO dimensions (left side) @ site " + std::to_string(i) + ": " + std::to_string(env_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO().dimension(0)));
//
//            if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
//                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (left side) @ site " + std::to_string(i) + ": "+ std::to_string(env2_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO().dimension(0)));
//
//            if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(0))
//                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (left side) @ site " + std::to_string(i) + ": " + std::to_string(env2_it->block.dimension(3)) + " " + std::to_string(mpo_it->get()->MPO().dimension(0)));
//
//
//            mps_it++;
//            mps_nx++;
//            env_it++;
//            env2_it++;
//            mpo_it++;
//            i++;
//        }
//    }

//    {
//        tools::log->trace("\tChecking matrix sizes on the center");
//        //Check center
//        if(state.MPS_C.dimension(0) != state.MPS_L.back().get_chiR())
//            throw std::runtime_error("Mismatch in center bond matrix dimension: " + std::to_string(state.MPS_C.dimension(0)) + " " + std::to_string(state.MPS_L.back().get_chiR()));
//        if(state.MPS_C.dimension(0) != state.MPS_R.front().get_chiL())
//            throw std::runtime_error("Mismatch in center bond matrix dimension: " + std::to_string(state.MPS_C.dimension(0)) + " " + std::to_string(state.MPS_R.front().get_chiL()));

//    }

//    {
//        tools::log->trace("\tChecking matrix sizes on the right side");
//        auto mps_it  = state.MPS_R.rbegin();
//        auto mps_nx  = state.MPS_R.rbegin();
//        auto env_it  = state.ENV_R.rbegin();
//        auto env2_it = state.ENV2_R.rbegin();
//        auto mpo_it  = state.MPO_R.rbegin();
//        std::advance(mps_nx,1);
//        auto i = state.get_length()-1;
//        while(
//                mps_it  != state.MPS_R.rend() and
//                mps_nx  != state.MPS_R.rend() and
//                env_it  != state.ENV_R.rend() and
//                env2_it != state.ENV2_R.rend() and
//                mpo_it  != state.MPO_R.rend())
//        {
//            if(mps_it->get_chiL() != mps_nx->get_chiR())
//                throw std::runtime_error("Mismatch in adjacent MPS dimensions (right side) @ site " + std::to_string(i) + ": "+ std::to_string(mps_nx->get_chiR()) + " " + std::to_string(mps_it->get_chiL()));
//
//            if(mps_it->get_position() != env_it->get_position())
//                throw std::runtime_error("Mismatch in MPS and ENV positions (right side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(env_it->get_position()));
//
//            if(mps_it->get_position() != state.get_length() - (env_it->size + 1))
//                throw std::runtime_error("Mismatch in MPS position and ENV size + 1 (right side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(state.get_length() - (env_it->size + 1)));
//
//            if(mps_it->get_chiR() != env_it->block.dimension(0))
//                throw std::runtime_error("Mismatch in MPS and ENV dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_chiR()) + " " + std::to_string(env_it->block.dimension(0)));
//
//            if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
//                throw std::runtime_error("Mismatch in ENV and MPO dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(env_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO().dimension(1)));
//
//            if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
//                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(env2_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO().dimension(1)));
//
//            if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(1))
//                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(env2_it->block.dimension(3)) + " " + std::to_string(mpo_it->get()->MPO().dimension(1)));
//
//
//            mps_it++;
//            mps_nx++;
//            env_it++;
//            env2_it++;
//            mpo_it++;
//            i--;
//        }
//    }
    tools::log->trace("MPS OK");
}



