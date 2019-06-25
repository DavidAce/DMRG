//
// Created by david on 2019-02-17.
//
#include <iostream>
#include <iomanip>
#include <algorithms/class_simulation_state.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_mpo.h>
#include <mps_state/class_finite_chain_state.h>
#include <general/nmspc_quantum_mechanics.h>
#include <spdlog/spdlog.h>

void mpstools::finite::debug::check_integrity(const class_finite_chain_state &state,
                                              const class_simulation_state &sim_state)
{
    mpstools::log->info("Checking integrity...");
    try{
        check_integrity_of_mps(state);
    }catch(std::exception & ex){
        mpstools::finite::print::print_state(state) ;
        throw std::runtime_error("Integrity check of MPS failed: " + std::string(ex.what()));
    }
    try {
        mpstools::log->info("Checking norms");
        auto norm_chain = mpstools::finite::measure::norm(state);
        if(std::abs(norm_chain - 1.0) > 1e-10) {
            mpstools::log->warn("Norm of state far from unity: {}", norm_chain);
            throw std::runtime_error("Norm of state too far from unity: " + std::to_string(norm_chain));
        }
    }
    catch(std::exception & ex){
        throw std::runtime_error("Integrity check of norm failed: " + std::string(ex.what()));
    }
}



void mpstools::finite::debug::check_integrity_of_mps(const class_finite_chain_state &state){
    {
        mpstools::log->trace("Checking integrity of MPS");
        mpstools::log->trace("\tChecking system sizes");
        if(state.MPS_L.size() + state.MPS_R.size() != state.get_length() )
            throw std::runtime_error("Mismatch in MPS size: " + std::to_string(state.MPS_L.size() + state.MPS_R.size()) + " " + std::to_string(state.get_length()));

        if(state.ENV_L.size() + state.ENV_R.size() != state.get_length())
            throw std::runtime_error("Mismatch in ENV size: " + std::to_string(state.ENV_L.size() + state.ENV_R.size()) + " " + std::to_string(state.get_length()));

        if(state.ENV2_L.size() + state.ENV2_R.size() != state.get_length())
            throw std::runtime_error("Mismatch in ENV size: " + std::to_string(state.ENV_L.size() + state.ENV_R.size()) + " " + std::to_string(state.get_length()));

        if( state.ENV_L.back().size != state.get_position())
            throw std::runtime_error("Mismatch in ENV_L size and position: " + std::to_string(state.ENV_L.back().size) + " " + std::to_string(state.get_position()));

        if(state.ENV_R.front().size != state.get_length() - state.get_position() - 2)
            throw std::runtime_error("Mismatch in ENV_R size+1 and length-position: " + std::to_string(state.ENV_R.front().size) + " " + std::to_string(state.get_length() - state.get_position()-2));

        mpstools::log->trace("\tChecking matrix sizes on the left side");
        //Check left side of the state
        auto mps_it  = state.MPS_L.begin();
        auto mps_nx  = state.MPS_L.begin();
        auto env_it  = state.ENV_L.begin();
        auto env2_it = state.ENV2_L.begin();
        auto mpo_it  = state.MPO_L.begin();
        std::advance(mps_nx,1);
        int i = 0;
        while(
            mps_it  != state.MPS_L.end() and
            mps_nx  != state.MPS_L.end() and
            env_it  != state.ENV_L.end() and
            env2_it != state.ENV2_L.end() and
            mpo_it  != state.MPO_L.end()
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

            if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
                throw std::runtime_error("Mismatch in ENV and MPO dimensions (left side) @ site " + std::to_string(i) + ": " + std::to_string(env_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO().dimension(0)));

            if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (left side) @ site " + std::to_string(i) + ": "+ std::to_string(env2_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO().dimension(0)));

            if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(0))
                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (left side) @ site " + std::to_string(i) + ": " + std::to_string(env2_it->block.dimension(3)) + " " + std::to_string(mpo_it->get()->MPO().dimension(0)));


            mps_it++;
            mps_nx++;
            env_it++;
            env2_it++;
            mpo_it++;
            i++;
        }
    }

    {
        mpstools::log->trace("\tChecking matrix sizes on the center");
        //Check center
        if(state.MPS_C.dimension(0) != state.MPS_L.back().get_chiR())
            throw std::runtime_error("Mismatch in center bond matrix dimension: " + std::to_string(state.MPS_C.dimension(0)) + " " + std::to_string(state.MPS_L.back().get_chiR()));
        if(state.MPS_C.dimension(0) != state.MPS_R.front().get_chiL())
            throw std::runtime_error("Mismatch in center bond matrix dimension: " + std::to_string(state.MPS_C.dimension(0)) + " " + std::to_string(state.MPS_R.front().get_chiL()));

    }

    {
        mpstools::log->trace("\tChecking matrix sizes on the right side");
        auto mps_it  = state.MPS_R.rbegin();
        auto mps_nx  = state.MPS_R.rbegin();
        auto env_it  = state.ENV_R.rbegin();
        auto env2_it = state.ENV2_R.rbegin();
        auto mpo_it  = state.MPO_R.rbegin();
        std::advance(mps_nx,1);
        auto i = state.get_length()-1;
        while(
            mps_it  != state.MPS_R.rend() and
            mps_nx  != state.MPS_R.rend() and
            env_it  != state.ENV_R.rend() and
            env2_it != state.ENV2_R.rend() and
            mpo_it  != state.MPO_R.rend())
        {
            if(mps_it->get_chiL() != mps_nx->get_chiR())
                throw std::runtime_error("Mismatch in adjacent MPS dimensions (right side) @ site " + std::to_string(i) + ": "+ std::to_string(mps_nx->get_chiR()) + " " + std::to_string(mps_it->get_chiL()));

            if(mps_it->get_position() != env_it->get_position())
                throw std::runtime_error("Mismatch in MPS and ENV positions (right side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(env_it->get_position()));

            if(mps_it->get_position() != state.get_length() - (env_it->size + 1))
                throw std::runtime_error("Mismatch in MPS position and ENV size + 1 (right side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_position()) + " " + std::to_string(state.get_length() - (env_it->size + 1)));

            if(mps_it->get_chiR() != env_it->block.dimension(0))
                throw std::runtime_error("Mismatch in MPS and ENV dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(mps_it->get_chiR()) + " " + std::to_string(env_it->block.dimension(0)));

            if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
                throw std::runtime_error("Mismatch in ENV and MPO dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(env_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO().dimension(1)));

            if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(env2_it->block.dimension(2)) + " " + std::to_string(mpo_it->get()->MPO().dimension(1)));

            if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(1))
                throw std::runtime_error("Mismatch in ENV2 and MPO dimensions (right side) @ site " + std::to_string(i) + ": " + std::to_string(env2_it->block.dimension(3)) + " " + std::to_string(mpo_it->get()->MPO().dimension(1)));


            mps_it++;
            mps_nx++;
            env_it++;
            env2_it++;
            mpo_it++;
            i--;
        }
    }
    mpstools::log->trace("MPS OK");
}





void mpstools::finite::debug::print_parity_properties(const class_finite_chain_state &state) {
    mpstools::log->info("Printing parity properties");

    mpstools::log->info("\tComputing spin components");
    const auto sx = mpstools::finite::measure::spin_component(state,qm::spinOneHalf::sx);
    const auto sy = mpstools::finite::measure::spin_component(state,qm::spinOneHalf::sy);
    const auto sz = mpstools::finite::measure::spin_component(state,qm::spinOneHalf::sz);
    mpstools::log->info("\t<psi | sx | psi>                = {:0.16f}", sx);
    mpstools::log->info("\t<psi | sy | psi>                = {:0.16f}", sy);
    mpstools::log->info("\t<psi | sz | psi>                = {:0.16f}", sz);

    mpstools::log->info("\tComputing parity projected states");
    auto state_up_x = mpstools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sx, 1);
    auto state_dn_x = mpstools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sx, -1);
    auto state_up_y = mpstools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sy, 1);
    auto state_dn_y = mpstools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sy, -1);
    auto state_up_z = mpstools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sz, 1);
    auto state_dn_z = mpstools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sz, -1);


    mpstools::log->info("\tMore spin components");
    mpstools::log->info("\t<psi_up_x | sx | psi_up_x>      = {:0.16f}", mpstools::finite::measure::spin_component(state_up_x,qm::spinOneHalf::sx));
    mpstools::log->info("\t<psi_up_x | sy | psi_up_x>      = {:0.16f}", mpstools::finite::measure::spin_component(state_up_x,qm::spinOneHalf::sy));
    mpstools::log->info("\t<psi_up_x | sz | psi_up_x>      = {:0.16f}", mpstools::finite::measure::spin_component(state_up_x,qm::spinOneHalf::sz));
    mpstools::log->info("\t<psi_dn_x | sx | psi_dn_x>      = {:0.16f}", mpstools::finite::measure::spin_component(state_dn_x,qm::spinOneHalf::sx));
    mpstools::log->info("\t<psi_dn_x | sy | psi_dn_x>      = {:0.16f}", mpstools::finite::measure::spin_component(state_dn_x,qm::spinOneHalf::sy));
    mpstools::log->info("\t<psi_dn_x | sz | psi_dn_x>      = {:0.16f}", mpstools::finite::measure::spin_component(state_dn_x,qm::spinOneHalf::sz));
    mpstools::log->info("\t<psi_up_y | sx | psi_up_y>      = {:0.16f}", mpstools::finite::measure::spin_component(state_up_y,qm::spinOneHalf::sx));
    mpstools::log->info("\t<psi_up_y | sy | psi_up_y>      = {:0.16f}", mpstools::finite::measure::spin_component(state_up_y,qm::spinOneHalf::sy));
    mpstools::log->info("\t<psi_up_y | sz | psi_up_y>      = {:0.16f}", mpstools::finite::measure::spin_component(state_up_y,qm::spinOneHalf::sz));
    mpstools::log->info("\t<psi_dn_y | sx | psi_dn_y>      = {:0.16f}", mpstools::finite::measure::spin_component(state_dn_y,qm::spinOneHalf::sx));
    mpstools::log->info("\t<psi_dn_y | sy | psi_dn_y>      = {:0.16f}", mpstools::finite::measure::spin_component(state_dn_y,qm::spinOneHalf::sy));
    mpstools::log->info("\t<psi_dn_y | sz | psi_dn_y>      = {:0.16f}", mpstools::finite::measure::spin_component(state_dn_y,qm::spinOneHalf::sz));
    mpstools::log->info("\t<psi_up_z | sx | psi_up_z>      = {:0.16f}", mpstools::finite::measure::spin_component(state_up_z,qm::spinOneHalf::sx));
    mpstools::log->info("\t<psi_up_z | sy | psi_up_z>      = {:0.16f}", mpstools::finite::measure::spin_component(state_up_z,qm::spinOneHalf::sy));
    mpstools::log->info("\t<psi_up_z | sz | psi_up_z>      = {:0.16f}", mpstools::finite::measure::spin_component(state_up_z,qm::spinOneHalf::sz));
    mpstools::log->info("\t<psi_dn_z | sx | psi_dn_z>      = {:0.16f}", mpstools::finite::measure::spin_component(state_dn_z,qm::spinOneHalf::sx));
    mpstools::log->info("\t<psi_dn_z | sy | psi_dn_z>      = {:0.16f}", mpstools::finite::measure::spin_component(state_dn_z,qm::spinOneHalf::sy));
    mpstools::log->info("\t<psi_dn_z | sz | psi_dn_z>      = {:0.16f}", mpstools::finite::measure::spin_component(state_dn_z,qm::spinOneHalf::sz));


    mpstools::log->info("\tNormalization check");
    mpstools::log->info("\tComputing overlaps");
    mpstools::log->info("\t<psi_up_x|psi_up_x>             = {:0.16f}", mpstools::finite::ops::overlap(state_up_x,state_up_x));
    mpstools::log->info("\t<psi_dn_x|psi_dn_x>             = {:0.16f}", mpstools::finite::ops::overlap(state_dn_x,state_dn_x));
    mpstools::log->info("\t<psi_up_y|psi_up_y>             = {:0.16f}", mpstools::finite::ops::overlap(state_up_y,state_up_y));
    mpstools::log->info("\t<psi_dn_y|psi_dn_y>             = {:0.16f}", mpstools::finite::ops::overlap(state_dn_y,state_dn_y));
    mpstools::log->info("\t<psi_up_z|psi_up_z>             = {:0.16f}", mpstools::finite::ops::overlap(state_up_z,state_up_z));
    mpstools::log->info("\t<psi_dn_z|psi_dn_z>             = {:0.16f}", mpstools::finite::ops::overlap(state_dn_z,state_dn_z));

    mpstools::log->info("\tOverlaps with original state");
    mpstools::log->info("\t<psi|psi_up_x>                  = {:0.16f}", mpstools::finite::ops::overlap(state,state_up_x));
    mpstools::log->info("\t<psi|psi_dn_x>                  = {:0.16f}", mpstools::finite::ops::overlap(state,state_dn_x));
    mpstools::log->info("\t<psi|psi_up_y>                  = {:0.16f}", mpstools::finite::ops::overlap(state,state_up_y));
    mpstools::log->info("\t<psi|psi_dn_y>                  = {:0.16f}", mpstools::finite::ops::overlap(state,state_dn_y));
    mpstools::log->info("\t<psi|psi_up_z>                  = {:0.16f}", mpstools::finite::ops::overlap(state,state_up_z));
    mpstools::log->info("\t<psi|psi_dn_z>                  = {:0.16f}", mpstools::finite::ops::overlap(state,state_dn_z));

    mpstools::log->info("\tOverlaps between up/down sectors");
    mpstools::log->info("\t<psi_dn_x|psi_up_x>             ={:0.16f}", mpstools::finite::ops::overlap(state_dn_x,state_up_x));
    mpstools::log->info("\t<psi_dn_y|psi_up_y>             ={:0.16f}", mpstools::finite::ops::overlap(state_dn_y,state_up_y));
    mpstools::log->info("\t<psi_dn_z|psi_up_z>             ={:0.16f}", mpstools::finite::ops::overlap(state_dn_z,state_up_z));
    mpstools::log->info("\tOverlaps between different direction sectors");
    mpstools::log->info("\t<psi_up_x|psi_up_y>             = {:0.16f}" ,mpstools::finite::ops::overlap(state_up_x,state_up_y));
    mpstools::log->info("\t<psi_up_x|psi_up_z>             = {:0.16f}" ,mpstools::finite::ops::overlap(state_up_x,state_up_z));
    mpstools::log->info("\tGenerating single hamiltonian MPO list");

    auto & Ledge = state.ENV_L.front().block;
    auto & Redge = state.ENV_R.back().block ;
    auto & Ledge2 = state.ENV2_L.front().block;
    auto & Redge2 = state.ENV2_R.back().block ;

    auto hamiltonian_mpos = mpstools::finite::ops::make_mpo_list(state.MPO_L,state.MPO_R);
    mpstools::log->info("Computing expectation values");
    double energy      = mpstools::finite::ops::expectation_value(state     ,state     ,hamiltonian_mpos, Ledge,Redge);
    double energy_up_x = mpstools::finite::ops::expectation_value(state_up_x,state_up_x,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_x = mpstools::finite::ops::expectation_value(state_dn_x,state_dn_x,hamiltonian_mpos, Ledge,Redge);
    double energy_up_y = mpstools::finite::ops::expectation_value(state_up_y,state_up_y,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_y = mpstools::finite::ops::expectation_value(state_dn_y,state_dn_y,hamiltonian_mpos, Ledge,Redge);
    double energy_up_z = mpstools::finite::ops::expectation_value(state_up_z,state_up_z,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_z = mpstools::finite::ops::expectation_value(state_dn_z,state_dn_z,hamiltonian_mpos, Ledge,Redge);
    mpstools::log->info("Computing variances");
    double variance      = mpstools::finite::ops::exp_sq_value(state     ,state     ,hamiltonian_mpos, Ledge2,Redge2) - energy      * energy ;
    double variance_up_x = mpstools::finite::ops::exp_sq_value(state_up_x,state_up_x,hamiltonian_mpos, Ledge2,Redge2) - energy_up_x * energy_up_x;
    double variance_dn_x = mpstools::finite::ops::exp_sq_value(state_dn_x,state_dn_x,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_x * energy_dn_x;
    double variance_up_y = mpstools::finite::ops::exp_sq_value(state_up_y,state_up_y,hamiltonian_mpos, Ledge2,Redge2) - energy_up_y * energy_up_y;
    double variance_dn_y = mpstools::finite::ops::exp_sq_value(state_dn_y,state_dn_y,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_y * energy_dn_y;
    double variance_up_z = mpstools::finite::ops::exp_sq_value(state_up_z,state_up_z,hamiltonian_mpos, Ledge2,Redge2) - energy_up_z * energy_up_z;
    double variance_dn_z = mpstools::finite::ops::exp_sq_value(state_dn_z,state_dn_z,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_z * energy_dn_z;

    auto state_copy2 = state;
    mpstools::finite::ops::apply_mpos(state_copy2,hamiltonian_mpos,Ledge,Redge);
    double overlap_H = mpstools::finite::ops::overlap(state_copy2,state);
    mpstools::log->info("\tEnergy per site");
    mpstools::log->info("\t<psi     | H/L |psi     >       = {:0.16f}" ,  energy    /state.get_length());
    mpstools::log->info("\t<psi     | H/L  psi     >       = {:0.16f}" ,  overlap_H /state.get_length() );
    mpstools::log->info("\t<psi_up_x| H/L |psi_up_x>       = {:0.16f}" ,  energy_up_x / state.get_length());
    mpstools::log->info("\t<psi_dn_x| H/L |psi_dn_x>       = {:0.16f}" ,  energy_dn_x / state.get_length());
    mpstools::log->info("\t<psi_up_y| H/L |psi_up_y>       = {:0.16f}" ,  energy_up_y / state.get_length());
    mpstools::log->info("\t<psi_dn_y| H/L |psi_dn_y>       = {:0.16f}" ,  energy_dn_y / state.get_length());
    mpstools::log->info("\t<psi_up_z| H/L |psi_up_z>       = {:0.16f}" ,  energy_up_z / state.get_length());
    mpstools::log->info("\t<psi_dn_z| H/L |psi_dn_z>       = {:0.16f}" ,  energy_dn_z / state.get_length());

    mpstools::log->info("\tVariance per site");
    mpstools::log->info("\t<psi     | (H2-E2)/L |psi     > = {:0.16f} | log10 = {:0.16f}", variance / state.get_length()     , std::log10(variance / state.get_length())     );
    mpstools::log->info("\t<psi_up_x| (H2-E2)/L |psi_up_x> = {:0.16f} | log10 = {:0.16f}", variance_up_x / state.get_length(), std::log10(variance_up_x / state.get_length()));
    mpstools::log->info("\t<psi_dn_x| (H2-E2)/L |psi_dn_x> = {:0.16f} | log10 = {:0.16f}", variance_dn_x / state.get_length(), std::log10(variance_dn_x / state.get_length()));
    mpstools::log->info("\t<psi_up_y| (H2-E2)/L |psi_up_y> = {:0.16f} | log10 = {:0.16f}", variance_up_y / state.get_length(), std::log10(variance_up_y / state.get_length()));
    mpstools::log->info("\t<psi_dn_y| (H2-E2)/L |psi_dn_y> = {:0.16f} | log10 = {:0.16f}", variance_dn_y / state.get_length(), std::log10(variance_dn_y / state.get_length()));
    mpstools::log->info("\t<psi_up_z| (H2-E2)/L |psi_up_z> = {:0.16f} | log10 = {:0.16f}", variance_up_z / state.get_length(), std::log10(variance_up_z / state.get_length()));
    mpstools::log->info("\t<psi_dn_z| (H2-E2)/L |psi_dn_z> = {:0.16f} | log10 = {:0.16f}", variance_dn_z / state.get_length(), std::log10(variance_dn_z / state.get_length()));

    mpstools::log->info("\tMidchain entanglement entropies");
    mpstools::log->info("\tS(L/2) psi                      = {:0.16f}" ,  mpstools::finite::measure::entanglement_entropy_midchain(state     ));
    mpstools::log->info("\tS(L/2) psi_up_x                 = {:0.16f}" ,  mpstools::finite::measure::entanglement_entropy_midchain(state_up_x));
    mpstools::log->info("\tS(L/2) psi_dn_x                 = {:0.16f}" ,  mpstools::finite::measure::entanglement_entropy_midchain(state_dn_x));
    mpstools::log->info("\tS(L/2) psi_up_y                 = {:0.16f}" ,  mpstools::finite::measure::entanglement_entropy_midchain(state_up_y));
    mpstools::log->info("\tS(L/2) psi_dn_y                 = {:0.16f}" ,  mpstools::finite::measure::entanglement_entropy_midchain(state_dn_y));
    mpstools::log->info("\tS(L/2) psi_up_z                 = {:0.16f}" ,  mpstools::finite::measure::entanglement_entropy_midchain(state_up_z));
    mpstools::log->info("\tS(L/2) psi_dn_z                 = {:0.16f}" ,  mpstools::finite::measure::entanglement_entropy_midchain(state_dn_z));



//    rebuild_environments(state_copy);
//
//    rebuild_superblock(state_copy, superblock);


}



void mpstools::finite::debug::check_normalization_routine(const class_finite_chain_state &state){
    mpstools::log->info("Checking normalization routine");
    mpstools::log->info("\t Generating Pauli Identity mpo");

    auto [mpo,L,R] = class_mpo::pauli_mpo(3*qm::spinOneHalf::Id);
    auto state_3ID = state;
    mpstools::log->info("\t Measuring original norm");
    auto norm_3ID     = mpstools::finite::measure::norm(state_3ID);
    mpstools::log->info("\t Measuring original overlap");
    auto overlap_3ID  = mpstools::finite::ops::overlap(state,state_3ID);
    std::cout << std::setprecision(16) << std::endl;
    std::cout << "Norm 3ID    = " << norm_3ID << std::endl;
    std::cout << "Overlap 3ID = " << overlap_3ID << std::endl;

    mpstools::log->info("\t Applying Pauli Identity mpo");
    mpstools::finite::ops::apply_mpo(state_3ID,mpo,L,R);
    mpstools::log->info("\t Measuring new norm");
    norm_3ID     = mpstools::finite::measure::norm(state_3ID);
    mpstools::log->info("\t Measuring new overlap");

    overlap_3ID  = mpstools::finite::ops::overlap(state,state_3ID);
    std::cout << std::setprecision(16) << std::endl;
    std::cout << "Norm 3ID    = " << norm_3ID << std::endl;
    std::cout << "Overlap 3ID = " << overlap_3ID << std::endl;
    mpstools::log->info("\t Normalizing state");

    mpstools::finite::mps::normalize(state_3ID);
    mpstools::log->info("\t Measuring new norm");
    norm_3ID     = mpstools::finite::measure::norm(state_3ID);
    mpstools::log->info("\t Measuring new overlap");
    overlap_3ID  = mpstools::finite::ops::overlap(state,state_3ID);
    std::cout << "Norm 3ID    = " << norm_3ID << std::endl;
    std::cout << "Overlap 3ID = " << overlap_3ID << std::endl;



    mpstools::log->info("\t Generating Pauli  sx up/dn mpos");
    auto [mpo_up,L_up,R_up] = class_mpo::parity_projector_mpos(qm::spinOneHalf::sx, state.get_length() ,1);
    auto [mpo_dn,L_dn,R_dn] = class_mpo::parity_projector_mpos(qm::spinOneHalf::sx, state.get_length() ,-1);
    auto state_sx_up    = state;
    auto state_sx_dn    = state;
    auto state_sx_dn_up = state;
    auto state_sx_up_up = state;
    auto state_sx_dn_dn = state;
    mpstools::log->info("\t Applying Pauli sx up/dn mpos");
    mpstools::finite::ops::apply_mpos(state_sx_up,mpo_up,L_up,R_up);
    mpstools::finite::ops::apply_mpos(state_sx_dn,mpo_dn,L_dn,R_dn);

    mpstools::finite::ops::apply_mpos(state_sx_dn_up,mpo_up,L_dn,R_dn);
    mpstools::finite::ops::apply_mpos(state_sx_dn_up,mpo_dn,L_dn,R_dn);

    mpstools::finite::ops::apply_mpos(state_sx_up_up,mpo_up,L_dn,R_dn);
    mpstools::finite::ops::apply_mpos(state_sx_up_up,mpo_up,L_dn,R_dn);

    mpstools::finite::ops::apply_mpos(state_sx_dn_dn,mpo_dn,L_dn,R_dn);
    mpstools::finite::ops::apply_mpos(state_sx_dn_dn,mpo_dn,L_dn,R_dn);

    mpstools::log->info("\t Measuring new norms");
    auto norm_sx_up     = mpstools::finite::measure::norm(state_sx_up);
    auto norm_sx_dn     = mpstools::finite::measure::norm(state_sx_dn);
    auto norm_sx_dn_up  = mpstools::finite::measure::norm(state_sx_dn_up);
    auto norm_sx_up_up  = mpstools::finite::measure::norm(state_sx_up_up);
    auto norm_sx_dn_dn  = mpstools::finite::measure::norm(state_sx_dn_dn);
    auto overlap_sx_up_up_up = mpstools::finite::ops::overlap(state_sx_up,state_sx_up_up);
    auto overlap_sx_dn_dn_dn = mpstools::finite::ops::overlap(state_sx_dn,state_sx_dn_dn);
    std::cout << "<P+   psi | P+   psi>      = " << norm_sx_up << std::endl;
    std::cout << "<P-   psi | P-   psi>      = " << norm_sx_dn << std::endl;
    std::cout << "<P-P+ psi | P-P+ psi>      = " << norm_sx_dn_up << std::endl;
    std::cout << "<P+ psi   | P+P+ psi>      = " << overlap_sx_up_up_up << std::endl;
    std::cout << "<P- psi   | P-P- psi>      = " << overlap_sx_dn_dn_dn << std::endl;

    mpstools::log->info("\t Normalizing states");
    mpstools::finite::mps::normalize(state_sx_up);
    mpstools::finite::mps::normalize(state_sx_dn);
    mpstools::finite::mps::normalize(state_sx_dn_up);
    mpstools::finite::mps::normalize(state_sx_up_up);
    mpstools::finite::mps::normalize(state_sx_dn_dn);
    mpstools::log->info("\t Measuring new norms");
    norm_sx_up     = mpstools::finite::measure::norm(state_sx_up);
    norm_sx_dn     = mpstools::finite::measure::norm(state_sx_dn);
    norm_sx_dn_up  = mpstools::finite::measure::norm(state_sx_dn_up);
    norm_sx_up_up  = mpstools::finite::measure::norm(state_sx_up_up);
    norm_sx_dn_dn  = mpstools::finite::measure::norm(state_sx_dn_dn);
    std::cout << "<N(P+   psi) psi | N(P+   psi)>      = " << norm_sx_up << std::endl;
    std::cout << "<N(P-   psi) psi | N(P-   psi)>      = " << norm_sx_dn << std::endl;
    std::cout << "<N(P-P+ psi) psi | N(P-P+ psi)>      = " << norm_sx_dn_up << std::endl;
    std::cout << "<N(P+P+ psi) psi | N(P+P+ psi)>      = " << norm_sx_up_up << std::endl;
    std::cout << "<N(P-P- psi) psi | N(P-P- psi)>      = " << norm_sx_dn_dn << std::endl;

    mpstools::log->info("\t Measuring new overlap");
    auto overlap_sx_up  = mpstools::finite::ops::overlap(state,state_sx_up);
    auto overlap_sx_dn  = mpstools::finite::ops::overlap(state,state_sx_dn);
    auto overlap_sx_dn_up = mpstools::finite::ops::overlap(state,state_sx_dn_up);
    auto overlap_sx_up_up = mpstools::finite::ops::overlap(state,state_sx_up_up);
    auto overlap_sx_dn_dn = mpstools::finite::ops::overlap(state,state_sx_dn_dn);
    overlap_sx_up_up_up = mpstools::finite::ops::overlap(state_sx_up,state_sx_up_up);
    overlap_sx_dn_dn_dn = mpstools::finite::ops::overlap(state_sx_dn,state_sx_dn_dn);
    std::cout << "<psi | N(P+   psi)>      = " << overlap_sx_up << std::endl;
    std::cout << "<psi | N(P-   psi)>      = " << overlap_sx_dn << std::endl;
    std::cout << "<psi | N(P-P+ psi)>      = " << overlap_sx_dn_up << std::endl;
    std::cout << "<psi | N(P+P+ psi)>      = " << overlap_sx_up_up << std::endl;
    std::cout << "<psi | N(P-P- psi)>      = " << overlap_sx_dn_dn << std::endl;
    std::cout << "<N(P+ psi) | N(P+P+ psi)>      = " << overlap_sx_up_up_up << std::endl;
    std::cout << "<N(P- psi) | N(P-P- psi)>      = " << overlap_sx_dn_dn_dn << std::endl;
}
