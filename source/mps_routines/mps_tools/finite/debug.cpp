//
// Created by david on 2019-02-17.
//
#include <iostream>
#include <iomanip>
#include <algorithms/class_simulation_state.h>
#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/class_superblock.h>
#include <general/nmspc_quantum_mechanics.h>
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
        spdlog::debug("\tChecking system size");

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

        spdlog::debug("\tChecking matrix sizes on the left side");
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
        spdlog::debug("\tChecking matrix sizes on the center");
        //Check center
        if(state.get_MPS_C().dimension(0) != state.get_MPS_L().back().get_chiR())
            throw std::runtime_error("Mismatch in center bond matrix dimension: " + std::to_string(state.get_MPS_C().dimension(0)) + " " + std::to_string(state.get_MPS_L().back().get_chiR()));
        if(state.get_MPS_C().dimension(0) != state.get_MPS_R().front().get_chiL())
            throw std::runtime_error("Mismatch in center bond matrix dimension: " + std::to_string(state.get_MPS_C().dimension(0)) + " " + std::to_string(state.get_MPS_R().front().get_chiL()));

    }

    {
        spdlog::debug("\tChecking matrix sizes on the right side");
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

    spdlog::debug("MPS OK!");


}





void MPS_Tools::Finite::Debug::print_parity_properties(const class_finite_chain_state &state) {
    spdlog::info("Printing parity properties");

    spdlog::info("\tComputing spin components");
    const auto sx = MPS_Tools::Finite::Measure::spin_component(state,qm::spinOneHalf::sx);
    const auto sy = MPS_Tools::Finite::Measure::spin_component(state,qm::spinOneHalf::sy);
    const auto sz = MPS_Tools::Finite::Measure::spin_component(state,qm::spinOneHalf::sz);
    spdlog::info("\t<psi | sx | psi>                = {:0.16f}", sx);
    spdlog::info("\t<psi | sy | psi>                = {:0.16f}", sy);
    spdlog::info("\t<psi | sz | psi>                = {:0.16f}", sz);

    spdlog::info("\tComputing parity projected states");
    auto state_up_x = MPS_Tools::Finite::Ops::get_parity_projected_state(state,qm::spinOneHalf::sx, 1);
    auto state_dn_x = MPS_Tools::Finite::Ops::get_parity_projected_state(state,qm::spinOneHalf::sx, -1);
    auto state_up_y = MPS_Tools::Finite::Ops::get_parity_projected_state(state,qm::spinOneHalf::sy, 1);
    auto state_dn_y = MPS_Tools::Finite::Ops::get_parity_projected_state(state,qm::spinOneHalf::sy, -1);
    auto state_up_z = MPS_Tools::Finite::Ops::get_parity_projected_state(state,qm::spinOneHalf::sz, 1);
    auto state_dn_z = MPS_Tools::Finite::Ops::get_parity_projected_state(state,qm::spinOneHalf::sz, -1);


    spdlog::info("\tMore spin components");
    spdlog::info("\t<psi_up_x | sx | psi_up_x>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_up_x,qm::spinOneHalf::sx));
    spdlog::info("\t<psi_up_x | sy | psi_up_x>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_up_x,qm::spinOneHalf::sy));
    spdlog::info("\t<psi_up_x | sz | psi_up_x>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_up_x,qm::spinOneHalf::sz));
    spdlog::info("\t<psi_dn_x | sx | psi_dn_x>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_dn_x,qm::spinOneHalf::sx));
    spdlog::info("\t<psi_dn_x | sy | psi_dn_x>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_dn_x,qm::spinOneHalf::sy));
    spdlog::info("\t<psi_dn_x | sz | psi_dn_x>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_dn_x,qm::spinOneHalf::sz));
    spdlog::info("\t<psi_up_y | sx | psi_up_y>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_up_y,qm::spinOneHalf::sx));
    spdlog::info("\t<psi_up_y | sy | psi_up_y>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_up_y,qm::spinOneHalf::sy));
    spdlog::info("\t<psi_up_y | sz | psi_up_y>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_up_y,qm::spinOneHalf::sz));
    spdlog::info("\t<psi_dn_y | sx | psi_dn_y>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_dn_y,qm::spinOneHalf::sx));
    spdlog::info("\t<psi_dn_y | sy | psi_dn_y>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_dn_y,qm::spinOneHalf::sy));
    spdlog::info("\t<psi_dn_y | sz | psi_dn_y>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_dn_y,qm::spinOneHalf::sz));
    spdlog::info("\t<psi_up_z | sx | psi_up_z>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_up_z,qm::spinOneHalf::sx));
    spdlog::info("\t<psi_up_z | sy | psi_up_z>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_up_z,qm::spinOneHalf::sy));
    spdlog::info("\t<psi_up_z | sz | psi_up_z>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_up_z,qm::spinOneHalf::sz));
    spdlog::info("\t<psi_dn_z | sx | psi_dn_z>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_dn_z,qm::spinOneHalf::sx));
    spdlog::info("\t<psi_dn_z | sy | psi_dn_z>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_dn_z,qm::spinOneHalf::sy));
    spdlog::info("\t<psi_dn_z | sz | psi_dn_z>      = {:0.16f}", MPS_Tools::Finite::Measure::spin_component(state_dn_z,qm::spinOneHalf::sz));


    spdlog::info("\tNormalization check");
    spdlog::info("\tComputing overlaps");
    spdlog::info("\t<psi_up_x|psi_up_x>             = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state_up_x,state_up_x));
    spdlog::info("\t<psi_dn_x|psi_dn_x>             = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state_dn_x,state_dn_x));
    spdlog::info("\t<psi_up_y|psi_up_y>             = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state_up_y,state_up_y));
    spdlog::info("\t<psi_dn_y|psi_dn_y>             = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state_dn_y,state_dn_y));
    spdlog::info("\t<psi_up_z|psi_up_z>             = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state_up_z,state_up_z));
    spdlog::info("\t<psi_dn_z|psi_dn_z>             = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state_dn_z,state_dn_z));

    spdlog::info("\tOverlaps with original state");
    spdlog::info("\t<psi|psi_up_x>                  = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state,state_up_x));
    spdlog::info("\t<psi|psi_dn_x>                  = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state,state_dn_x));
    spdlog::info("\t<psi|psi_up_y>                  = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state,state_up_y));
    spdlog::info("\t<psi|psi_dn_y>                  = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state,state_dn_y));
    spdlog::info("\t<psi|psi_up_z>                  = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state,state_up_z));
    spdlog::info("\t<psi|psi_dn_z>                  = {:0.16f}", MPS_Tools::Finite::Ops::overlap(state,state_dn_z));

    spdlog::info("\tOverlaps between up/down sectors");
    spdlog::info("\t<psi_dn_x|psi_up_x>             ={:0.16f}", MPS_Tools::Finite::Ops::overlap(state_dn_x,state_up_x));
    spdlog::info("\t<psi_dn_y|psi_up_y>             ={:0.16f}", MPS_Tools::Finite::Ops::overlap(state_dn_y,state_up_y));
    spdlog::info("\t<psi_dn_z|psi_up_z>             ={:0.16f}", MPS_Tools::Finite::Ops::overlap(state_dn_z,state_up_z));
    spdlog::info("\tOverlaps between different direction sectors");
    spdlog::info("\t<psi_up_x|psi_up_y>             = {:0.16f}" ,MPS_Tools::Finite::Ops::overlap(state_up_x,state_up_y));
    spdlog::info("\t<psi_up_x|psi_up_z>             = {:0.16f}" ,MPS_Tools::Finite::Ops::overlap(state_up_x,state_up_z));
    spdlog::info("\tGenerating single hamiltonian MPO list");

    auto & Ledge = state.get_ENV_L().front().block;
    auto & Redge = state.get_ENV_R().back().block ;
    auto & Ledge2 = state.get_ENV2_L().front().block;
    auto & Redge2 = state.get_ENV2_R().back().block ;

    auto hamiltonian_mpos = MPS_Tools::Finite::Ops::make_mpo_list(state.get_MPO_L(),state.get_MPO_R());
    spdlog::info("Computing expectation values");
    double energy      = MPS_Tools::Finite::Ops::expectation_value(state     ,state     ,hamiltonian_mpos, Ledge,Redge);
    double energy_up_x = MPS_Tools::Finite::Ops::expectation_value(state_up_x,state_up_x,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_x = MPS_Tools::Finite::Ops::expectation_value(state_dn_x,state_dn_x,hamiltonian_mpos, Ledge,Redge);
    double energy_up_y = MPS_Tools::Finite::Ops::expectation_value(state_up_y,state_up_y,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_y = MPS_Tools::Finite::Ops::expectation_value(state_dn_y,state_dn_y,hamiltonian_mpos, Ledge,Redge);
    double energy_up_z = MPS_Tools::Finite::Ops::expectation_value(state_up_z,state_up_z,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_z = MPS_Tools::Finite::Ops::expectation_value(state_dn_z,state_dn_z,hamiltonian_mpos, Ledge,Redge);
    spdlog::info("Computing variances");
    double variance      = MPS_Tools::Finite::Ops::exp_sq_value(state     ,state     ,hamiltonian_mpos, Ledge2,Redge2) - energy      * energy ;
    double variance_up_x = MPS_Tools::Finite::Ops::exp_sq_value(state_up_x,state_up_x,hamiltonian_mpos, Ledge2,Redge2) - energy_up_x * energy_up_x;
    double variance_dn_x = MPS_Tools::Finite::Ops::exp_sq_value(state_dn_x,state_dn_x,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_x * energy_dn_x;
    double variance_up_y = MPS_Tools::Finite::Ops::exp_sq_value(state_up_y,state_up_y,hamiltonian_mpos, Ledge2,Redge2) - energy_up_y * energy_up_y;
    double variance_dn_y = MPS_Tools::Finite::Ops::exp_sq_value(state_dn_y,state_dn_y,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_y * energy_dn_y;
    double variance_up_z = MPS_Tools::Finite::Ops::exp_sq_value(state_up_z,state_up_z,hamiltonian_mpos, Ledge2,Redge2) - energy_up_z * energy_up_z;
    double variance_dn_z = MPS_Tools::Finite::Ops::exp_sq_value(state_dn_z,state_dn_z,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_z * energy_dn_z;

    auto state_copy2 = state;
    MPS_Tools::Finite::Ops::apply_mpos(state_copy2,hamiltonian_mpos,Ledge,Redge);
    double overlap_H = MPS_Tools::Finite::Ops::overlap(state_copy2,state);
    spdlog::info("\tEnergy per site");
    spdlog::info("\t<psi     | H/L |psi     >       = {:0.16f}" ,  energy    /state.get_length());
    spdlog::info("\t<psi     | H/L  psi     >       = {:0.16f}" ,  overlap_H /state.get_length() );
    spdlog::info("\t<psi_up_x| H/L |psi_up_x>       = {:0.16f}" ,  energy_up_x / state.get_length());
    spdlog::info("\t<psi_dn_x| H/L |psi_dn_x>       = {:0.16f}" ,  energy_dn_x / state.get_length());
    spdlog::info("\t<psi_up_y| H/L |psi_up_y>       = {:0.16f}" ,  energy_up_y / state.get_length());
    spdlog::info("\t<psi_dn_y| H/L |psi_dn_y>       = {:0.16f}" ,  energy_dn_y / state.get_length());
    spdlog::info("\t<psi_up_z| H/L |psi_up_z>       = {:0.16f}" ,  energy_up_z / state.get_length());
    spdlog::info("\t<psi_dn_z| H/L |psi_dn_z>       = {:0.16f}" ,  energy_dn_z / state.get_length());

    spdlog::info("\tVariance per site");
    spdlog::info("\t<psi     | (H2-E2)/L |psi     > = {:0.16f} | log10 = {:0.16f}", variance / state.get_length()     , std::log10(variance / state.get_length())     );
    spdlog::info("\t<psi_up_x| (H2-E2)/L |psi_up_x> = {:0.16f} | log10 = {:0.16f}", variance_up_x / state.get_length(), std::log10(variance_up_x / state.get_length()));
    spdlog::info("\t<psi_dn_x| (H2-E2)/L |psi_dn_x> = {:0.16f} | log10 = {:0.16f}", variance_dn_x / state.get_length(), std::log10(variance_dn_x / state.get_length()));
    spdlog::info("\t<psi_up_y| (H2-E2)/L |psi_up_y> = {:0.16f} | log10 = {:0.16f}", variance_up_y / state.get_length(), std::log10(variance_up_y / state.get_length()));
    spdlog::info("\t<psi_dn_y| (H2-E2)/L |psi_dn_y> = {:0.16f} | log10 = {:0.16f}", variance_dn_y / state.get_length(), std::log10(variance_dn_y / state.get_length()));
    spdlog::info("\t<psi_up_z| (H2-E2)/L |psi_up_z> = {:0.16f} | log10 = {:0.16f}", variance_up_z / state.get_length(), std::log10(variance_up_z / state.get_length()));
    spdlog::info("\t<psi_dn_z| (H2-E2)/L |psi_dn_z> = {:0.16f} | log10 = {:0.16f}", variance_dn_z / state.get_length(), std::log10(variance_dn_z / state.get_length()));

    spdlog::info("\tMidchain entanglement entropies");
    spdlog::info("\tS(L/2) psi                      = {:0.16f}" ,  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state     ));
    spdlog::info("\tS(L/2) psi_up_x                 = {:0.16f}" ,  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_up_x));
    spdlog::info("\tS(L/2) psi_dn_x                 = {:0.16f}" ,  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_dn_x));
    spdlog::info("\tS(L/2) psi_up_y                 = {:0.16f}" ,  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_up_y));
    spdlog::info("\tS(L/2) psi_dn_y                 = {:0.16f}" ,  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_dn_y));
    spdlog::info("\tS(L/2) psi_up_z                 = {:0.16f}" ,  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_up_z));
    spdlog::info("\tS(L/2) psi_dn_z                 = {:0.16f}" ,  MPS_Tools::Finite::Measure::midchain_entanglement_entropy(state_dn_z));



//    rebuild_environments(state_copy);
//
//    rebuild_superblock(state_copy, superblock);


}



void MPS_Tools::Finite::Debug::check_normalization_routine(const class_finite_chain_state &state){
    spdlog::info("Checking normalization routine");
    spdlog::info("\t Generating Pauli Identity mpo");

    auto [mpo,L,R] = class_mpo::pauli_mpo(3*qm::spinOneHalf::Id);
    auto state_3ID = state;
    spdlog::info("\t Measuring original norm");
    auto norm_3ID     = MPS_Tools::Finite::Measure::norm(state_3ID);
    spdlog::info("\t Measuring original overlap");
    auto overlap_3ID  = MPS_Tools::Finite::Ops::overlap(state,state_3ID);
    std::cout << std::setprecision(16) << std::endl;
    std::cout << "Norm 3ID    = " << norm_3ID << std::endl;
    std::cout << "Overlap 3ID = " << overlap_3ID << std::endl;

    spdlog::info("\t Applying Pauli Identity mpo");
    MPS_Tools::Finite::Ops::apply_mpo(state_3ID,mpo,L,R);
    spdlog::info("\t Measuring new norm");
    norm_3ID     = MPS_Tools::Finite::Measure::norm(state_3ID);
    spdlog::info("\t Measuring new overlap");

    overlap_3ID  = MPS_Tools::Finite::Ops::overlap(state,state_3ID);
    std::cout << std::setprecision(16) << std::endl;
    std::cout << "Norm 3ID    = " << norm_3ID << std::endl;
    std::cout << "Overlap 3ID = " << overlap_3ID << std::endl;
    spdlog::info("\t Normalizing state");

    MPS_Tools::Finite::Ops::normalize_chain(state_3ID);
    spdlog::info("\t Measuring new norm");
    norm_3ID     = MPS_Tools::Finite::Measure::norm(state_3ID);
    spdlog::info("\t Measuring new overlap");
    overlap_3ID  = MPS_Tools::Finite::Ops::overlap(state,state_3ID);
    std::cout << "Norm 3ID    = " << norm_3ID << std::endl;
    std::cout << "Overlap 3ID = " << overlap_3ID << std::endl;



    spdlog::info("\t Generating Pauli  sx up/dn mpos");
    auto [mpo_up,L_up,R_up] = class_mpo::parity_projector_mpos(qm::spinOneHalf::sx, state.get_length() ,1);
    auto [mpo_dn,L_dn,R_dn] = class_mpo::parity_projector_mpos(qm::spinOneHalf::sx, state.get_length() ,-1);
    auto state_sx_up    = state;
    auto state_sx_dn    = state;
    auto state_sx_dn_up = state;
    auto state_sx_up_up = state;
    auto state_sx_dn_dn = state;
    spdlog::info("\t Applying Pauli sx up/dn mpos");
    MPS_Tools::Finite::Ops::apply_mpos(state_sx_up,mpo_up,L_up,R_up);
    MPS_Tools::Finite::Ops::apply_mpos(state_sx_dn,mpo_dn,L_dn,R_dn);

    MPS_Tools::Finite::Ops::apply_mpos(state_sx_dn_up,mpo_up,L_dn,R_dn);
    MPS_Tools::Finite::Ops::apply_mpos(state_sx_dn_up,mpo_dn,L_dn,R_dn);

    MPS_Tools::Finite::Ops::apply_mpos(state_sx_up_up,mpo_up,L_dn,R_dn);
    MPS_Tools::Finite::Ops::apply_mpos(state_sx_up_up,mpo_up,L_dn,R_dn);

    MPS_Tools::Finite::Ops::apply_mpos(state_sx_dn_dn,mpo_dn,L_dn,R_dn);
    MPS_Tools::Finite::Ops::apply_mpos(state_sx_dn_dn,mpo_dn,L_dn,R_dn);

    spdlog::info("\t Measuring new norms");
    auto norm_sx_up     = MPS_Tools::Finite::Measure::norm(state_sx_up);
    auto norm_sx_dn     = MPS_Tools::Finite::Measure::norm(state_sx_dn);
    auto norm_sx_dn_up  = MPS_Tools::Finite::Measure::norm(state_sx_dn_up);
    auto norm_sx_up_up  = MPS_Tools::Finite::Measure::norm(state_sx_up_up);
    auto norm_sx_dn_dn  = MPS_Tools::Finite::Measure::norm(state_sx_dn_dn);
    auto overlap_sx_up_up_up = MPS_Tools::Finite::Ops::overlap(state_sx_up,state_sx_up_up);
    auto overlap_sx_dn_dn_dn = MPS_Tools::Finite::Ops::overlap(state_sx_dn,state_sx_dn_dn);
    std::cout << "<P+   psi | P+   psi>      = " << norm_sx_up << std::endl;
    std::cout << "<P-   psi | P-   psi>      = " << norm_sx_dn << std::endl;
    std::cout << "<P-P+ psi | P-P+ psi>      = " << norm_sx_dn_up << std::endl;
    std::cout << "<P+ psi   | P+P+ psi>      = " << overlap_sx_up_up_up << std::endl;
    std::cout << "<P- psi   | P-P- psi>      = " << overlap_sx_dn_dn_dn << std::endl;

    spdlog::info("\t Normalizing states");
    MPS_Tools::Finite::Ops::normalize_chain(state_sx_up);
    MPS_Tools::Finite::Ops::normalize_chain(state_sx_dn);
    MPS_Tools::Finite::Ops::normalize_chain(state_sx_dn_up);
    MPS_Tools::Finite::Ops::normalize_chain(state_sx_up_up);
    MPS_Tools::Finite::Ops::normalize_chain(state_sx_dn_dn);
    spdlog::info("\t Measuring new norms");
    norm_sx_up     = MPS_Tools::Finite::Measure::norm(state_sx_up);
    norm_sx_dn     = MPS_Tools::Finite::Measure::norm(state_sx_dn);
    norm_sx_dn_up  = MPS_Tools::Finite::Measure::norm(state_sx_dn_up);
    norm_sx_up_up  = MPS_Tools::Finite::Measure::norm(state_sx_up_up);
    norm_sx_dn_dn  = MPS_Tools::Finite::Measure::norm(state_sx_dn_dn);
    std::cout << "<N(P+   psi) psi | N(P+   psi)>      = " << norm_sx_up << std::endl;
    std::cout << "<N(P-   psi) psi | N(P-   psi)>      = " << norm_sx_dn << std::endl;
    std::cout << "<N(P-P+ psi) psi | N(P-P+ psi)>      = " << norm_sx_dn_up << std::endl;
    std::cout << "<N(P+P+ psi) psi | N(P+P+ psi)>      = " << norm_sx_up_up << std::endl;
    std::cout << "<N(P-P- psi) psi | N(P-P- psi)>      = " << norm_sx_dn_dn << std::endl;

    spdlog::info("\t Measuring new overlap");
    auto overlap_sx_up  = MPS_Tools::Finite::Ops::overlap(state,state_sx_up);
    auto overlap_sx_dn  = MPS_Tools::Finite::Ops::overlap(state,state_sx_dn);
    auto overlap_sx_dn_up = MPS_Tools::Finite::Ops::overlap(state,state_sx_dn_up);
    auto overlap_sx_up_up = MPS_Tools::Finite::Ops::overlap(state,state_sx_up_up);
    auto overlap_sx_dn_dn = MPS_Tools::Finite::Ops::overlap(state,state_sx_dn_dn);
    overlap_sx_up_up_up = MPS_Tools::Finite::Ops::overlap(state_sx_up,state_sx_up_up);
    overlap_sx_dn_dn_dn = MPS_Tools::Finite::Ops::overlap(state_sx_dn,state_sx_dn_dn);
    std::cout << "<psi | N(P+   psi)>      = " << overlap_sx_up << std::endl;
    std::cout << "<psi | N(P-   psi)>      = " << overlap_sx_dn << std::endl;
    std::cout << "<psi | N(P-P+ psi)>      = " << overlap_sx_dn_up << std::endl;
    std::cout << "<psi | N(P+P+ psi)>      = " << overlap_sx_up_up << std::endl;
    std::cout << "<psi | N(P-P- psi)>      = " << overlap_sx_dn_dn << std::endl;
    std::cout << "<N(P+ psi) | N(P+P+ psi)>      = " << overlap_sx_up_up_up << std::endl;
    std::cout << "<N(P- psi) | N(P-P- psi)>      = " << overlap_sx_dn_dn_dn << std::endl;
}
