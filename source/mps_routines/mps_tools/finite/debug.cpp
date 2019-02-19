//
// Created by david on 2019-02-17.
//
#include <iostream>

#include <algorithms/class_simulation_state.h>
#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/class_superblock.h>
#include <spdlog/spdlog.h>

void MPS_Tools::Finite::Debug::check_integrity(const class_finite_chain_state &state,
                                               const class_superblock &superblock, class_simulation_state &sim_state)
{
    try{
        check_integrity_of_sim(state, superblock, sim_state);
    }catch(std::exception & ex){
        std::cout << sim_state << std::endl;
        throw std::runtime_error("Integrity check of constants failed: " + std::string(ex.what()));
    }
    try{
        check_integrity_of_mps(state);
    }catch(std::exception & ex){
        MPS_Tools::Finite::Print::print_state(state) ;
        throw std::runtime_error("Integrity check of MPS failed: " + std::string(ex.what()));
    }
    try {
        MPS_Tools::Finite::Measure::norm(state);
        MPS_Tools::Common::Measure::norm(superblock);
    }catch(std::exception & ex){
        throw std::runtime_error("Integrity check of norm failed: " + std::string(ex.what()));
    }
}



void MPS_Tools::Finite::Debug::check_integrity_of_sim(const class_finite_chain_state &state,
                                                      const class_superblock &superblock,
                                                      class_simulation_state &sim_state)
{
    spdlog::debug("Checking integrity of SIM");

    if(state.get_length() != superblock.get_length())
        throw std::runtime_error("Length mismatch in state and superblock: " + std::to_string(state.get_length()) + " " + std::to_string(superblock.get_length()));

    if(state.get_length()  != superblock.environment_size + 2 )
        throw std::runtime_error("Length mismatch in state and superblock env + 2: " + std::to_string(state.get_length()) + " " + std::to_string(superblock.environment_size + 2 ));

    if(state.get_position()  != superblock.Lblock->size )
        throw std::runtime_error("Mismatch in state position and sites in superblock env: " + std::to_string(state.get_position()) + " " + std::to_string(superblock.Lblock->size));


    if(state.get_position() != sim_state.position  )
        throw std::runtime_error("Mismatch in state and sim_state positions: " + std::to_string(state.get_position()) + " " + std::to_string(sim_state.position ));


    auto state_MPO_A      = Eigen::Map<Eigen::VectorXcd>(state.get_MPO_L().back()->MPO.data(),state.get_MPO_L().back()->MPO.size() );
    auto superblock_MPO_A = Eigen::Map<Eigen::VectorXcd>(superblock.HA->MPO.data(),superblock.HA->MPO.size() );
    auto state_MPO_B      = Eigen::Map<Eigen::VectorXcd>(state.get_MPO_R().front()->MPO.data(),state.get_MPO_R().front()->MPO.size() );
    auto superblock_MPO_B = Eigen::Map<Eigen::VectorXcd>(superblock.HB->MPO.data(),superblock.HB->MPO.size() );
    if(state_MPO_A != superblock_MPO_A )
        throw std::runtime_error("Mismatch in state and superblock MPOS (left side)");
    if(state_MPO_B != superblock_MPO_B )
        throw std::runtime_error("Mismatch in state and superblock MPOS (right side)");

}



void MPS_Tools::Finite::Debug::check_integrity_of_mps(const class_finite_chain_state &state){
    {
        spdlog::debug("Checking integrity of MPS");

        if(state.get_MPS_L().size() + state.get_MPS_R().size() != state.get_length() )
            throw std::runtime_error("Mismatch in MPS size: " + std::to_string(state.get_MPS_L().size() + state.get_MPS_R().size()) + " " + std::to_string(state.get_length()));

        if(state.get_ENV_L().size() + state.get_ENV_R().size() != state.get_length())
            throw std::runtime_error("Mismatch in ENV size: " + std::to_string(state.get_ENV_L().size() + state.get_ENV_R().size()) + " " + std::to_string(state.get_length()));

        if(state.get_ENV2_L().size() + state.get_ENV2_R().size() != state.get_length())
            throw std::runtime_error("Mismatch in ENV size: " + std::to_string(state.get_ENV_L().size() + state.get_ENV_R().size()) + " " + std::to_string(state.get_length()));

        if(state.get_ENV_L().back().size != state.get_position())
            throw std::runtime_error("Mismatch in ENV_L size and position: " + std::to_string(state.get_ENV_L().back().size) + " " + std::to_string(state.get_position()));

        if(state.get_ENV_R().front().size != state.get_length() - state.get_position() - 2)
            throw std::runtime_error("Mismatch in ENV_R size+1 and length-position: " + std::to_string(state.get_ENV_R().front().size) + " " + std::to_string(state.get_length() - state.get_position()-2));


        //Check left side of the chain
        auto mps_it  = state.get_MPS_L().begin();
        auto mps_nx  = state.get_MPS_L().begin();
        auto env_it  = state.get_ENV_L().begin();
        auto env2_it = state.get_ENV2_L().begin();
        auto mpo_it  = state.get_MPO_L().begin();
        std::advance(mps_nx,1);
        int i = 0;
        while(
            mps_it  != state.get_MPS_L().end() and
            mps_nx  != state.get_MPS_L().end() and
            env_it  != state.get_ENV_L().end() and
            env2_it != state.get_ENV2_L().end() and
            mpo_it  != state.get_MPO_L().end()
            )
        {
            if(mps_it->get_chiR() != mps_nx->get_chiL())
                throw std::runtime_error("Mismatch in adjacent MPS dimensions (left side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_chiR()) + " " + std::to_string(mps_nx->get_chiL()));

            if(mps_it->get_position() != env_it->get_position())
                throw std::runtime_error("Mismatch in MPS and ENV positions (left side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(env_it->get_position()));

            if(mps_it->get_position() != env_it->size)
                throw std::runtime_error("Mismatch in MPS position and ENV size (left side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(env_it->size));

            if(mps_it->get_chiL() != env_it->block.dimension(0))
                throw std::runtime_error("Mismatch in MPS and ENV dimensions (left side) @ site " + std::to_string(i) + ": " +  std::to_string(mps_it->get_chiL()) + " " + std::to_string(env_it->block.dimension(0)));


            if(env_it->block.dimension(2) != mpo_it->get()->MPO.dimension(0))
                throw std::runtime_error("Mismatch in ENV and MPO dimensions (left side) @ site " + std::to_string(i) + ": " + std::to_string(env_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO.dimension(0)));

            if(env2_it->block.dimension(2) != mpo_it->get()->MPO.dimension(0))
                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (left side) @ site " + std::to_string(i) + ": "+ std::to_string(env2_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO.dimension(0)));

            if(env2_it->block.dimension(3) != mpo_it->get()->MPO.dimension(0))
                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (left side) @ site " + std::to_string(i) + ": " + std::to_string(env2_it->block.dimension(3)) + " " + std::to_string(mpo_it->get()->MPO.dimension(0)));


            mps_it++;
            mps_nx++;
            env_it++;
            env2_it++;
            mpo_it++;
            i++;
        }
    }

    {
        //Check center
        if(state.get_MPS_C().dimension(0) != state.get_MPS_L().back().get_chiR())
            throw std::runtime_error("Mismatch in center bond matrix dimension: " + std::to_string(state.get_MPS_C().dimension(0)) + " " + std::to_string(state.get_MPS_L().back().get_chiR()));
        if(state.get_MPS_C().dimension(0) != state.get_MPS_R().front().get_chiL())
            throw std::runtime_error("Mismatch in center bond matrix dimension: " + std::to_string(state.get_MPS_C().dimension(0)) + " " + std::to_string(state.get_MPS_R().front().get_chiL()));

    }

    {
        auto mps_it  = state.get_MPS_R().rbegin();
        auto mps_nx  = state.get_MPS_R().rbegin();
        auto env_it  = state.get_ENV_R().rbegin();
        auto env2_it = state.get_ENV2_R().rbegin();
        auto mpo_it  = state.get_MPO_R().rbegin();
        std::advance(mps_nx,1);
        int i = state.get_length()-1;
        while(
            mps_it  != state.get_MPS_R().rend() and
            mps_nx  != state.get_MPS_R().rend() and
            env_it  != state.get_ENV_R().rend() and
            env2_it != state.get_ENV2_R().rend() and
            mpo_it  != state.get_MPO_R().rend())
        {
            if(mps_it->get_chiL() != mps_nx->get_chiR())
                throw std::runtime_error("Mismatch in adjacent MPS dimensions (right side) @ site " + std::to_string(i) + ": "+ std::to_string(mps_nx->get_chiR()) + " " + std::to_string(mps_it->get_chiL()));

            if(mps_it->get_position() != env_it->get_position())
                throw std::runtime_error("Mismatch in MPS and ENV positions (right side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(env_it->get_position()));

            if(mps_it->get_position() != state.get_length() - (env_it->size + 1))
                throw std::runtime_error("Mismatch in MPS position and ENV size + 1 (right side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(state.get_length() - (env_it->size + 1)));

            if(mps_it->get_chiR() != env_it->block.dimension(0))
                throw std::runtime_error("Mismatch in MPS and ENV dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_chiR()) + " " + std::to_string(env_it->block.dimension(0)));

            if(env_it->block.dimension(2) != mpo_it->get()->MPO.dimension(1))
                throw std::runtime_error("Mismatch in ENV and MPO dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(env_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO.dimension(1)));

            if(env2_it->block.dimension(2) != mpo_it->get()->MPO.dimension(1))
                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(env2_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO.dimension(1)));

            if(env2_it->block.dimension(3) != mpo_it->get()->MPO.dimension(1))
                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(env2_it->block.dimension(3)) + " " + std::to_string(mpo_it->get()->MPO.dimension(1)));


            mps_it++;
            mps_nx++;
            env_it++;
            env2_it++;
            mpo_it++;
            i--;
        }
    }





}
