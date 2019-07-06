//
// Created by david on 2019-02-17.
//
#include <iostream>
#include <iomanip>
#include <simulation/class_simulation_status.h>
#include <state/tools/nmspc_tools.h>
#include <state/class_finite_state.h>
#include <general/nmspc_quantum_mechanics.h>
#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>

void tools::finite::debug::check_integrity(const class_finite_state &state)
{
    tools::log->info("Checking integrity of state");
    state.unset_measurements();
    try{
        check_integrity_of_mps(state);
    }catch(std::exception & ex){
        tools::finite::print::print_state(state) ;
        throw std::runtime_error(fmt::format("Integrity check of MPS failed: {}", ex.what()));
    }
    try{
        check_integrity_of_mpo(state);
    }catch(std::exception & ex){
        tools::finite::print::print_state(state) ;
        throw std::runtime_error(fmt::format("Integrity check of MPO failed: {}", ex.what()));
    }
    try {
        tools::log->debug("Checking norms");
        auto norm_chain = tools::finite::measure::norm(state);
        if(std::abs(norm_chain - 1.0) > 1e-12) {
            throw std::runtime_error(fmt::format("Norm of state too far from unity: {}",norm_chain));
        }
    }
    catch(std::exception & ex){
        throw std::runtime_error(fmt::format("Integrity check of norm failed: {}", ex.what()));
    }
    state.unset_measurements();
}



void tools::finite::debug::check_integrity_of_mps(const class_finite_state &state){
    {
        tools::log->trace("Checking integrity of MPS");
        if(state.MPS_L.size() + state.MPS_R.size() != state.get_length() )
            throw std::runtime_error(fmt::format("Mismatch in MPS sizes: {} + {} != {}", state.MPS_L.size(), state.MPS_R.size(), state.get_length()));

        if(state.ENV_L.size() + state.ENV_R.size() != state.get_length())
            throw std::runtime_error(fmt::format("Mismatch in ENV sizes: {} + {} != {}", state.ENV_L.size(), state.ENV_R.size(),state.get_length()));

        if(state.ENV2_L.size() + state.ENV2_R.size() != state.get_length())
            throw std::runtime_error(fmt::format("Mismatch in ENV2 sizes: {} + {} != {}", state.ENV2_L.size(), state.ENV2_R.size(),state.get_length()));

        if( state.ENV_L.back().sites != state.get_position())
            throw std::runtime_error(fmt::format("Mismatch in ENV_L sites and position: {} != {}", state.ENV_L.back().sites, state.get_position()));

        if(state.ENV_R.front().sites != state.get_length() - state.get_position() - 2)
            throw std::runtime_error(fmt::format("Mismatch in ENV_R size+1 and length-position: {} != {}", state.ENV_R.front().sites, state.get_length() - state.get_position() - 2));

        using VectorType = const Eigen::Matrix<class_finite_state::Scalar,Eigen::Dynamic,1>;
        for (size_t pos = 1; pos < state.get_length(); pos++){

            auto & mps_left = state.get_MPS(pos-1);
            auto & mps_here = state.get_MPS(pos);
            auto & mpo_left = state.get_MPO(pos-1);
            auto & mpo_here = state.get_MPO(pos);

            // Check for validity
            if(not Eigen::Map<VectorType>(mps_here.get_G().data(), mps_here.get_G().size()).allFinite()) {
                std::cerr << "G: \n" << mps_here.get_G() << std::endl;
                throw std::runtime_error(fmt::format("Inf's or nan's in MPS G @ pos {}", pos));
            }

            if(not Eigen::Map<VectorType>(mps_here.get_L().data(), mps_here.get_L().size()).allFinite()){
                std::cerr << "L: \n" << mps_here.get_L() << std::endl;
                throw std::runtime_error(fmt::format("Inf's or nan's in MPS L @ pos {}", pos));
            }

            // Check positions
            if(mps_here.get_position() != mpo_here.get_position())
                throw std::runtime_error(fmt::format("Mismatch in MPS and MPO positions @ pos {}: {} != {}", pos, mps_here.get_position(), mpo_here.get_position()));

            if(mps_here.get_position() - mps_left.get_position() != 1)
                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS positions @ pos {}: {} - {} != 1", pos, mps_here.get_position() , mps_left.get_position()));

            if(mps_left.get_chiR() != mps_here.get_chiL())
                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS dimensions @ pos {}: {} != {}", pos, mps_left.get_chiR() , mps_here.get_chiL()));


            if(mpo_left.MPO().dimension(1) != mpo_here.MPO().dimension(0))
                throw std::runtime_error(fmt::format("Mismatch in adjacent MPO dimensions @ pos {}: {} != {}", pos, mpo_left.MPO().dimension(1)  , mpo_here.MPO().dimension(0)));

            if(mps_here.get_spin_dim() != mpo_here.MPO().dimension(2))
                throw std::runtime_error(fmt::format("Mismatch in MPS and MPO spin dimensions @ pos {}: {} != {}", pos, mps_here.get_spin_dim() , mpo_here.MPO().dimension(2)));



        }


        //Check left side of the state
        auto mps_it  = state.MPS_L.begin();
        auto mps_nx  = state.MPS_L.begin();
        auto env_it  = state.ENV_L.begin();
        auto env_nx  = state.ENV_L.begin();
        auto env2_it = state.ENV2_L.begin();
        auto mpo_it  = state.MPO_L.begin();
        std::advance(mps_nx,1);
        std::advance(env_nx,1);
        int i = 0;
        while(
            mps_it  != state.MPS_L.end() and
            mps_nx  != state.MPS_L.end() and
            env_it  != state.ENV_L.end() and
            env_nx  != state.ENV_L.end() and
            env2_it != state.ENV2_L.end() and
            mpo_it  != state.MPO_L.end()
            )
        {
            if(mps_it->get_position() != env_it->get_position())
                throw std::runtime_error(fmt::format("Mismatch in MPS and ENV positions (left side) @ site {}: {} != {}", i, mps_it->get_position(), env_it->get_position()));

            if(mps_it->get_position() != env_it->sites)
                throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV size (left side) @ site {}: {} != {}", i, mps_it->get_position(), env_it->sites));

            if(mps_it->get_chiL() != env_it->block.dimension(0))
                throw std::runtime_error(fmt::format("Mismatch in MPS and ENV dimensions (left side) @ site {}: {} != {}", i,mps_it->get_chiL() , env_it->block.dimension(0)));

            if(env_nx->get_position() - env_it->get_position() != 1)
                throw std::runtime_error(fmt::format("Mismatch in adjacent ENV positions (left side) @ site {}: {} - {} != 1", i, env_nx->get_position(), env_it->get_position()));

            if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
                throw std::runtime_error(fmt::format("Mismatch in ENV and MPO dimensions (left side) @ site {}: {} != {}", i,env_it->block.dimension(2), mpo_it->get()->MPO().dimension(0)));


            if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
                throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (left side) @ site {}: {} != {}", i,env2_it->block.dimension(2), mpo_it->get()->MPO().dimension(0)));

            if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(0))
                throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (left side) @ site {}: {} != {}", i,env2_it->block.dimension(3) ,mpo_it->get()->MPO().dimension(0)));


            if(env2_it->sites != env_it->sites)
                throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV sites (left side) @ site {}: {} != {}", i,env2_it->sites != env_it->sites));

            if(env2_it->get_position() != env_it->get_position())
                throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV positions (left side) @ site {}: {} != {}", i,env2_it->get_position(), env_it->get_position()));

            mps_it++;
            mps_nx++;
            env_it++;
            env_nx++;
            env2_it++;
            mpo_it++;
            i++;
        }
    }

    {
        //Check center
        if(state.MPS_C.dimension(0) != state.MPS_L.back().get_chiR())
            throw std::runtime_error(fmt::format("Mismatch in center bond matrix dimension: {} != {}",state.MPS_C.dimension(0) , state.MPS_L.back().get_chiR()));
        if(state.MPS_C.dimension(0) != state.MPS_R.front().get_chiL())
            throw std::runtime_error(fmt::format("Mismatch in center bond matrix dimension: {} != {}",state.MPS_C.dimension(0) , state.MPS_R.front().get_chiL()));
    }

    {
        auto mps_it  = state.MPS_R.rbegin();
        auto mps_nx  = state.MPS_R.rbegin();
        auto env_it  = state.ENV_R.rbegin();
        auto env_nx  = state.ENV_R.rbegin();
        auto env2_it = state.ENV2_R.rbegin();
        auto mpo_it  = state.MPO_R.rbegin();
        std::advance(mps_nx,1);
        std::advance(env_nx,1);
        auto i = state.get_length()-1;
        while(
            mps_it  != state.MPS_R.rend() and
            mps_nx  != state.MPS_R.rend() and
            env_it  != state.ENV_R.rend() and
            env_nx  != state.ENV_R.rend() and
            env2_it != state.ENV2_R.rend() and
            mpo_it  != state.MPO_R.rend())
        {
            if(mps_it->get_chiL() != mps_nx->get_chiR())
                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS dimensions (right side) @ site {}: {} != {}", i,  mps_nx->get_chiR(), mps_it->get_chiL()));

            if(mps_it->get_position() - mps_nx->get_position() != 1 )
                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS positions (right side) @ site {}: {} - {} != 1", i, mps_it->get_position(), mps_nx->get_position()));

            if(mps_it->get_position() != env_it->get_position())
                throw std::runtime_error(fmt::format("Mismatch in MPS and ENV positions (right side) @ site {}: {} != {}", i,  mps_it->get_position(), env_it->get_position()));

            if(mps_it->get_position() != state.get_length() - (env_it->sites + 1))
                throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV size + 1 (right side) @ site {}: {} != {}", i,  mps_it->get_position(), state.get_length() - (env_it->sites + 1)));

            if(mps_it->get_position() != env_it->get_position())
                throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV position (right side) @ site {}: {} != {}", i,  mps_it->get_position(), env_it->get_position()));

            if(mps_it->get_chiR() != env_it->block.dimension(0))
                throw std::runtime_error(fmt::format("Mismatch in MPS and ENV dimensions (right side) @ site {}: {} != {}", i,  mps_it->get_chiR(), env_it->block.dimension(0)));

            if(env_it->get_position() - env_nx->get_position() != 1)
                throw std::runtime_error(fmt::format("Mismatch in adjacent ENV positions (right side) @ site {}: {} - {} != 1", i, env_it->get_position(), env_nx->get_position()));

            if(env2_it->sites != env_it->sites)
                throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV sites (right side) @ site {}: {} != {}", i,  env2_it->sites, env_it->sites));

            if(env2_it->get_position() != env_it->get_position())
                throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV positions (right side) @ site {}: {} != {}", i,  env2_it->get_position(), env_it->get_position()));


            if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
                throw std::runtime_error(fmt::format("Mismatch in ENV and MPO dimensions (right side) @ site {}: {} != {}", i,  env_it->block.dimension(2), mpo_it->get()->MPO().dimension(1)));

            if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
                throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (right side) @ site {}: {} != {}", i,  env2_it->block.dimension(2), mpo_it->get()->MPO().dimension(1)));

            if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(1))
                throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (right side) @ site {}: {} != {}", i,  env2_it->block.dimension(3), mpo_it->get()->MPO().dimension(1)));


            mps_it++;
            mps_nx++;
            env_it++;
            env_nx++;
            env2_it++;
            mpo_it++;
            i--;
        }
    }
    tools::log->trace("MPS OK");
}


void tools::finite::debug::check_integrity_of_mpo(const class_finite_state &state) {

    for (auto &mpo : state.MPO_L){
        auto site = mpo->get_position();
        if (not mpo->all_mpo_parameters_have_been_set){throw std::runtime_error(fmt::format("All parameters have not been set on MPO_L site: {}" ,site));}
    }
    for (auto &mpo : state.MPO_R){
        auto site = mpo->get_position();
        if (not mpo->all_mpo_parameters_have_been_set){throw std::runtime_error(fmt::format("All parameters have not been set on MPO_R site: {}" ,site));}
    }
}


void tools::finite::debug::print_parity_properties(const class_finite_state &state) {
    tools::log->debug("Printing parity properties");

    tools::log->debug("\tComputing spin components");
    const auto sx = tools::finite::measure::spin_component(state,qm::spinOneHalf::sx);
    const auto sy = tools::finite::measure::spin_component(state,qm::spinOneHalf::sy);
    const auto sz = tools::finite::measure::spin_component(state,qm::spinOneHalf::sz);
    tools::log->debug("\t<psi | sx | psi>                = {:0.16f}", sx);
    tools::log->debug("\t<psi | sy | psi>                = {:0.16f}", sy);
    tools::log->debug("\t<psi | sz | psi>                = {:0.16f}", sz);

    tools::log->debug("\tComputing parity projected states");
    auto state_up_x = tools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sx, 1);
    auto state_dn_x = tools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sx, -1);
    auto state_up_y = tools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sy, 1);
    auto state_dn_y = tools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sy, -1);
    auto state_up_z = tools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sz, 1);
    auto state_dn_z = tools::finite::ops::get_parity_projected_state(state,qm::spinOneHalf::sz, -1);


    tools::log->debug("\tMore spin components");
    tools::log->debug("\t<psi_up_x | sx | psi_up_x>      = {:0.16f}", tools::finite::measure::spin_component(state_up_x,qm::spinOneHalf::sx));
    tools::log->debug("\t<psi_up_x | sy | psi_up_x>      = {:0.16f}", tools::finite::measure::spin_component(state_up_x,qm::spinOneHalf::sy));
    tools::log->debug("\t<psi_up_x | sz | psi_up_x>      = {:0.16f}", tools::finite::measure::spin_component(state_up_x,qm::spinOneHalf::sz));
    tools::log->debug("\t<psi_dn_x | sx | psi_dn_x>      = {:0.16f}", tools::finite::measure::spin_component(state_dn_x,qm::spinOneHalf::sx));
    tools::log->debug("\t<psi_dn_x | sy | psi_dn_x>      = {:0.16f}", tools::finite::measure::spin_component(state_dn_x,qm::spinOneHalf::sy));
    tools::log->debug("\t<psi_dn_x | sz | psi_dn_x>      = {:0.16f}", tools::finite::measure::spin_component(state_dn_x,qm::spinOneHalf::sz));
    tools::log->debug("\t<psi_up_y | sx | psi_up_y>      = {:0.16f}", tools::finite::measure::spin_component(state_up_y,qm::spinOneHalf::sx));
    tools::log->debug("\t<psi_up_y | sy | psi_up_y>      = {:0.16f}", tools::finite::measure::spin_component(state_up_y,qm::spinOneHalf::sy));
    tools::log->debug("\t<psi_up_y | sz | psi_up_y>      = {:0.16f}", tools::finite::measure::spin_component(state_up_y,qm::spinOneHalf::sz));
    tools::log->debug("\t<psi_dn_y | sx | psi_dn_y>      = {:0.16f}", tools::finite::measure::spin_component(state_dn_y,qm::spinOneHalf::sx));
    tools::log->debug("\t<psi_dn_y | sy | psi_dn_y>      = {:0.16f}", tools::finite::measure::spin_component(state_dn_y,qm::spinOneHalf::sy));
    tools::log->debug("\t<psi_dn_y | sz | psi_dn_y>      = {:0.16f}", tools::finite::measure::spin_component(state_dn_y,qm::spinOneHalf::sz));
    tools::log->debug("\t<psi_up_z | sx | psi_up_z>      = {:0.16f}", tools::finite::measure::spin_component(state_up_z,qm::spinOneHalf::sx));
    tools::log->debug("\t<psi_up_z | sy | psi_up_z>      = {:0.16f}", tools::finite::measure::spin_component(state_up_z,qm::spinOneHalf::sy));
    tools::log->debug("\t<psi_up_z | sz | psi_up_z>      = {:0.16f}", tools::finite::measure::spin_component(state_up_z,qm::spinOneHalf::sz));
    tools::log->debug("\t<psi_dn_z | sx | psi_dn_z>      = {:0.16f}", tools::finite::measure::spin_component(state_dn_z,qm::spinOneHalf::sx));
    tools::log->debug("\t<psi_dn_z | sy | psi_dn_z>      = {:0.16f}", tools::finite::measure::spin_component(state_dn_z,qm::spinOneHalf::sy));
    tools::log->debug("\t<psi_dn_z | sz | psi_dn_z>      = {:0.16f}", tools::finite::measure::spin_component(state_dn_z,qm::spinOneHalf::sz));


    tools::log->debug("\tNormalization check");
    tools::log->debug("\tComputing overlaps");
    tools::log->debug("\t<psi_up_x|psi_up_x>             = {:0.16f}", tools::finite::ops::overlap(state_up_x,state_up_x));
    tools::log->debug("\t<psi_dn_x|psi_dn_x>             = {:0.16f}", tools::finite::ops::overlap(state_dn_x,state_dn_x));
    tools::log->debug("\t<psi_up_y|psi_up_y>             = {:0.16f}", tools::finite::ops::overlap(state_up_y,state_up_y));
    tools::log->debug("\t<psi_dn_y|psi_dn_y>             = {:0.16f}", tools::finite::ops::overlap(state_dn_y,state_dn_y));
    tools::log->debug("\t<psi_up_z|psi_up_z>             = {:0.16f}", tools::finite::ops::overlap(state_up_z,state_up_z));
    tools::log->debug("\t<psi_dn_z|psi_dn_z>             = {:0.16f}", tools::finite::ops::overlap(state_dn_z,state_dn_z));

    tools::log->debug("\tOverlaps with original state");
    tools::log->debug("\t<psi|psi_up_x>                  = {:0.16f}", tools::finite::ops::overlap(state,state_up_x));
    tools::log->debug("\t<psi|psi_dn_x>                  = {:0.16f}", tools::finite::ops::overlap(state,state_dn_x));
    tools::log->debug("\t<psi|psi_up_y>                  = {:0.16f}", tools::finite::ops::overlap(state,state_up_y));
    tools::log->debug("\t<psi|psi_dn_y>                  = {:0.16f}", tools::finite::ops::overlap(state,state_dn_y));
    tools::log->debug("\t<psi|psi_up_z>                  = {:0.16f}", tools::finite::ops::overlap(state,state_up_z));
    tools::log->debug("\t<psi|psi_dn_z>                  = {:0.16f}", tools::finite::ops::overlap(state,state_dn_z));

    tools::log->debug("\tOverlaps between up/down sectors");
    tools::log->debug("\t<psi_dn_x|psi_up_x>             ={:0.16f}", tools::finite::ops::overlap(state_dn_x,state_up_x));
    tools::log->debug("\t<psi_dn_y|psi_up_y>             ={:0.16f}", tools::finite::ops::overlap(state_dn_y,state_up_y));
    tools::log->debug("\t<psi_dn_z|psi_up_z>             ={:0.16f}", tools::finite::ops::overlap(state_dn_z,state_up_z));
    tools::log->debug("\tOverlaps between different direction sectors");
    tools::log->debug("\t<psi_up_x|psi_up_y>             = {:0.16f}" ,tools::finite::ops::overlap(state_up_x,state_up_y));
    tools::log->debug("\t<psi_up_x|psi_up_z>             = {:0.16f}" ,tools::finite::ops::overlap(state_up_x,state_up_z));
    tools::log->debug("\tGenerating single hamiltonian MPO list");

    auto & Ledge = state.ENV_L.front().block;
    auto & Redge = state.ENV_R.back().block ;
    auto & Ledge2 = state.ENV2_L.front().block;
    auto & Redge2 = state.ENV2_R.back().block ;

    auto hamiltonian_mpos = tools::finite::ops::make_mpo_list(state.MPO_L,state.MPO_R);
    tools::log->debug("Computing expectation values");
    double energy      = tools::finite::ops::expectation_value(state     ,state     ,hamiltonian_mpos, Ledge,Redge);
    double energy_up_x = tools::finite::ops::expectation_value(state_up_x,state_up_x,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_x = tools::finite::ops::expectation_value(state_dn_x,state_dn_x,hamiltonian_mpos, Ledge,Redge);
    double energy_up_y = tools::finite::ops::expectation_value(state_up_y,state_up_y,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_y = tools::finite::ops::expectation_value(state_dn_y,state_dn_y,hamiltonian_mpos, Ledge,Redge);
    double energy_up_z = tools::finite::ops::expectation_value(state_up_z,state_up_z,hamiltonian_mpos, Ledge,Redge);
    double energy_dn_z = tools::finite::ops::expectation_value(state_dn_z,state_dn_z,hamiltonian_mpos, Ledge,Redge);
    tools::log->debug("Computing variances");
    double variance      = tools::finite::ops::exp_sq_value(state     ,state     ,hamiltonian_mpos, Ledge2,Redge2) - energy      * energy ;
    double variance_up_x = tools::finite::ops::exp_sq_value(state_up_x,state_up_x,hamiltonian_mpos, Ledge2,Redge2) - energy_up_x * energy_up_x;
    double variance_dn_x = tools::finite::ops::exp_sq_value(state_dn_x,state_dn_x,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_x * energy_dn_x;
    double variance_up_y = tools::finite::ops::exp_sq_value(state_up_y,state_up_y,hamiltonian_mpos, Ledge2,Redge2) - energy_up_y * energy_up_y;
    double variance_dn_y = tools::finite::ops::exp_sq_value(state_dn_y,state_dn_y,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_y * energy_dn_y;
    double variance_up_z = tools::finite::ops::exp_sq_value(state_up_z,state_up_z,hamiltonian_mpos, Ledge2,Redge2) - energy_up_z * energy_up_z;
    double variance_dn_z = tools::finite::ops::exp_sq_value(state_dn_z,state_dn_z,hamiltonian_mpos, Ledge2,Redge2) - energy_dn_z * energy_dn_z;

    auto state_copy2 = state;
    tools::finite::ops::apply_mpos(state_copy2,hamiltonian_mpos,Ledge,Redge);
    double overlap_H = tools::finite::ops::overlap(state_copy2,state);
    tools::log->debug("\tEnergy per site");
    tools::log->debug("\t<psi     | H/L |psi     >       = {:0.16f}" ,  energy    /state.get_length());
    tools::log->debug("\t<psi     | H/L  psi     >       = {:0.16f}" ,  overlap_H /state.get_length() );
    tools::log->debug("\t<psi_up_x| H/L |psi_up_x>       = {:0.16f}" ,  energy_up_x / state.get_length());
    tools::log->debug("\t<psi_dn_x| H/L |psi_dn_x>       = {:0.16f}" ,  energy_dn_x / state.get_length());
    tools::log->debug("\t<psi_up_y| H/L |psi_up_y>       = {:0.16f}" ,  energy_up_y / state.get_length());
    tools::log->debug("\t<psi_dn_y| H/L |psi_dn_y>       = {:0.16f}" ,  energy_dn_y / state.get_length());
    tools::log->debug("\t<psi_up_z| H/L |psi_up_z>       = {:0.16f}" ,  energy_up_z / state.get_length());
    tools::log->debug("\t<psi_dn_z| H/L |psi_dn_z>       = {:0.16f}" ,  energy_dn_z / state.get_length());

    tools::log->debug("\tVariance per site");
    tools::log->debug("\t<psi     | (H2-E2)/L |psi     > = {:0.16f} | log10 = {:0.16f}", variance / state.get_length()     , std::log10(variance / state.get_length())     );
    tools::log->debug("\t<psi_up_x| (H2-E2)/L |psi_up_x> = {:0.16f} | log10 = {:0.16f}", variance_up_x / state.get_length(), std::log10(variance_up_x / state.get_length()));
    tools::log->debug("\t<psi_dn_x| (H2-E2)/L |psi_dn_x> = {:0.16f} | log10 = {:0.16f}", variance_dn_x / state.get_length(), std::log10(variance_dn_x / state.get_length()));
    tools::log->debug("\t<psi_up_y| (H2-E2)/L |psi_up_y> = {:0.16f} | log10 = {:0.16f}", variance_up_y / state.get_length(), std::log10(variance_up_y / state.get_length()));
    tools::log->debug("\t<psi_dn_y| (H2-E2)/L |psi_dn_y> = {:0.16f} | log10 = {:0.16f}", variance_dn_y / state.get_length(), std::log10(variance_dn_y / state.get_length()));
    tools::log->debug("\t<psi_up_z| (H2-E2)/L |psi_up_z> = {:0.16f} | log10 = {:0.16f}", variance_up_z / state.get_length(), std::log10(variance_up_z / state.get_length()));
    tools::log->debug("\t<psi_dn_z| (H2-E2)/L |psi_dn_z> = {:0.16f} | log10 = {:0.16f}", variance_dn_z / state.get_length(), std::log10(variance_dn_z / state.get_length()));

    tools::log->debug("\tMidchain entanglement entropies");
    tools::log->debug("\tS(L/2) psi                      = {:0.16f}" ,  tools::finite::measure::entanglement_entropy_midchain(state     ));
    tools::log->debug("\tS(L/2) psi_up_x                 = {:0.16f}" ,  tools::finite::measure::entanglement_entropy_midchain(state_up_x));
    tools::log->debug("\tS(L/2) psi_dn_x                 = {:0.16f}" ,  tools::finite::measure::entanglement_entropy_midchain(state_dn_x));
    tools::log->debug("\tS(L/2) psi_up_y                 = {:0.16f}" ,  tools::finite::measure::entanglement_entropy_midchain(state_up_y));
    tools::log->debug("\tS(L/2) psi_dn_y                 = {:0.16f}" ,  tools::finite::measure::entanglement_entropy_midchain(state_dn_y));
    tools::log->debug("\tS(L/2) psi_up_z                 = {:0.16f}" ,  tools::finite::measure::entanglement_entropy_midchain(state_up_z));
    tools::log->debug("\tS(L/2) psi_dn_z                 = {:0.16f}" ,  tools::finite::measure::entanglement_entropy_midchain(state_dn_z));



//    rebuild_environments(state_copy);
//
//    rebuild_superblock(state_copy, state);


}



void tools::finite::debug::check_normalization_routine(const class_finite_state &state){
    tools::log->debug("Checking normalization routine");
    tools::log->debug("\t Generating Pauli Identity mpo");

    auto [mpo,L,R] = qm::mpo::pauli_mpo(3*qm::spinOneHalf::Id);
    auto state_3ID = state;
    tools::log->debug("\t Measuring original norm");
    auto norm_3ID     = tools::finite::measure::norm(state_3ID);
    tools::log->debug("\t Measuring original overlap");
    auto overlap_3ID  = tools::finite::ops::overlap(state,state_3ID);
    std::cout << std::setprecision(16) << std::endl;
    std::cout << "Norm 3ID    = " << norm_3ID << std::endl;
    std::cout << "Overlap 3ID = " << overlap_3ID << std::endl;

    tools::log->debug("\t Applying Pauli Identity mpo");
    tools::finite::ops::apply_mpo(state_3ID,mpo,L,R);
    tools::log->debug("\t Measuring new norm");
    norm_3ID     = tools::finite::measure::norm(state_3ID);
    tools::log->debug("\t Measuring new overlap");

    overlap_3ID  = tools::finite::ops::overlap(state,state_3ID);
    std::cout << std::setprecision(16) << std::endl;
    std::cout << "Norm 3ID    = " << norm_3ID << std::endl;
    std::cout << "Overlap 3ID = " << overlap_3ID << std::endl;
    tools::log->debug("\t Normalizing state");

    tools::finite::mps::normalize(state_3ID);
    tools::log->debug("\t Measuring new norm");
    norm_3ID     = tools::finite::measure::norm(state_3ID);
    tools::log->debug("\t Measuring new overlap");
    overlap_3ID  = tools::finite::ops::overlap(state,state_3ID);
    std::cout << "Norm 3ID    = " << norm_3ID << std::endl;
    std::cout << "Overlap 3ID = " << overlap_3ID << std::endl;



    tools::log->debug("\t Generating Pauli  sx up/dn mpos");
    auto [mpo_up,L_up,R_up] = qm::mpo::parity_projector_mpos(qm::spinOneHalf::sx, state.get_length() ,1);
    auto [mpo_dn,L_dn,R_dn] = qm::mpo::parity_projector_mpos(qm::spinOneHalf::sx, state.get_length() ,-1);
    auto state_sx_up    = state;
    auto state_sx_dn    = state;
    auto state_sx_dn_up = state;
    auto state_sx_up_up = state;
    auto state_sx_dn_dn = state;
    tools::log->debug("\t Applying Pauli sx up/dn mpos");
    tools::finite::ops::apply_mpos(state_sx_up,mpo_up,L_up,R_up);
    tools::finite::ops::apply_mpos(state_sx_dn,mpo_dn,L_dn,R_dn);

    tools::finite::ops::apply_mpos(state_sx_dn_up,mpo_up,L_dn,R_dn);
    tools::finite::ops::apply_mpos(state_sx_dn_up,mpo_dn,L_dn,R_dn);

    tools::finite::ops::apply_mpos(state_sx_up_up,mpo_up,L_dn,R_dn);
    tools::finite::ops::apply_mpos(state_sx_up_up,mpo_up,L_dn,R_dn);

    tools::finite::ops::apply_mpos(state_sx_dn_dn,mpo_dn,L_dn,R_dn);
    tools::finite::ops::apply_mpos(state_sx_dn_dn,mpo_dn,L_dn,R_dn);

    tools::log->debug("\t Measuring new norms");
    auto norm_sx_up     = tools::finite::measure::norm(state_sx_up);
    auto norm_sx_dn     = tools::finite::measure::norm(state_sx_dn);
    auto norm_sx_dn_up  = tools::finite::measure::norm(state_sx_dn_up);
    auto norm_sx_up_up  = tools::finite::measure::norm(state_sx_up_up);
    auto norm_sx_dn_dn  = tools::finite::measure::norm(state_sx_dn_dn);
    auto overlap_sx_up_up_up = tools::finite::ops::overlap(state_sx_up,state_sx_up_up);
    auto overlap_sx_dn_dn_dn = tools::finite::ops::overlap(state_sx_dn,state_sx_dn_dn);
    std::cout << "<P+   psi | P+   psi>      = " << norm_sx_up << std::endl;
    std::cout << "<P-   psi | P-   psi>      = " << norm_sx_dn << std::endl;
    std::cout << "<P-P+ psi | P-P+ psi>      = " << norm_sx_dn_up << std::endl;
    std::cout << "<P+ psi   | P+P+ psi>      = " << overlap_sx_up_up_up << std::endl;
    std::cout << "<P- psi   | P-P- psi>      = " << overlap_sx_dn_dn_dn << std::endl;

    tools::log->debug("\t Normalizing states");
    tools::finite::mps::normalize(state_sx_up);
    tools::finite::mps::normalize(state_sx_dn);
    tools::finite::mps::normalize(state_sx_dn_up);
    tools::finite::mps::normalize(state_sx_up_up);
    tools::finite::mps::normalize(state_sx_dn_dn);
    tools::log->debug("\t Measuring new norms");
    norm_sx_up     = tools::finite::measure::norm(state_sx_up);
    norm_sx_dn     = tools::finite::measure::norm(state_sx_dn);
    norm_sx_dn_up  = tools::finite::measure::norm(state_sx_dn_up);
    norm_sx_up_up  = tools::finite::measure::norm(state_sx_up_up);
    norm_sx_dn_dn  = tools::finite::measure::norm(state_sx_dn_dn);
    std::cout << "<N(P+   psi) psi | N(P+   psi)>      = " << norm_sx_up << std::endl;
    std::cout << "<N(P-   psi) psi | N(P-   psi)>      = " << norm_sx_dn << std::endl;
    std::cout << "<N(P-P+ psi) psi | N(P-P+ psi)>      = " << norm_sx_dn_up << std::endl;
    std::cout << "<N(P+P+ psi) psi | N(P+P+ psi)>      = " << norm_sx_up_up << std::endl;
    std::cout << "<N(P-P- psi) psi | N(P-P- psi)>      = " << norm_sx_dn_dn << std::endl;

    tools::log->debug("\t Measuring new overlap");
    auto overlap_sx_up  = tools::finite::ops::overlap(state,state_sx_up);
    auto overlap_sx_dn  = tools::finite::ops::overlap(state,state_sx_dn);
    auto overlap_sx_dn_up = tools::finite::ops::overlap(state,state_sx_dn_up);
    auto overlap_sx_up_up = tools::finite::ops::overlap(state,state_sx_up_up);
    auto overlap_sx_dn_dn = tools::finite::ops::overlap(state,state_sx_dn_dn);
    overlap_sx_up_up_up = tools::finite::ops::overlap(state_sx_up,state_sx_up_up);
    overlap_sx_dn_dn_dn = tools::finite::ops::overlap(state_sx_dn,state_sx_dn_dn);
    std::cout << "<psi | N(P+   psi)>      = " << overlap_sx_up << std::endl;
    std::cout << "<psi | N(P-   psi)>      = " << overlap_sx_dn << std::endl;
    std::cout << "<psi | N(P-P+ psi)>      = " << overlap_sx_dn_up << std::endl;
    std::cout << "<psi | N(P+P+ psi)>      = " << overlap_sx_up_up << std::endl;
    std::cout << "<psi | N(P-P- psi)>      = " << overlap_sx_dn_dn << std::endl;
    std::cout << "<N(P+ psi) | N(P+P+ psi)>      = " << overlap_sx_up_up_up << std::endl;
    std::cout << "<N(P- psi) | N(P-P- psi)>      = " << overlap_sx_dn_dn_dn << std::endl;
}
