//
// Created by david on 2019-02-17.
//
#include <general/nmspc_quantum_mechanics.h>
#include <iostream>
#include <math/nmspc_math.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/edges/class_env_ene.h>
#include <tensors/edges/class_env_var.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/class_tensors_finite.h>
#include <config/nmspc_settings.h>

#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/debug.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/ops.h>
#include <tools/finite/print.h>

void tools::finite::debug::check_integrity(const class_tensors_finite & tensors)
{
    if(not settings::debug) return;

    tools::log->trace("Checking integrity of state");
//    state.clear_measurements();

    try{
        check_integrity(*tensors.state);
        check_integrity(*tensors.model);
        check_integrity(*tensors.edges);
        check_integrity(*tensors.state,*tensors.model, *tensors.edges);
    }
    catch(std::exception & ex){
        tools::finite::print::dimensions(tensors);
//        tools::finite::debug::print_parity_properties(state) ;
        throw std::runtime_error(fmt::format("Check failed: {}", ex.what()));
    }
//    tensors.clear_measurements();
}



void tools::finite::debug::check_integrity(const class_state_finite & state){
//    if constexpr (not settings::debug) return;
//    tools::log->trace("Checking integrity of MPS");
//    tools::common::profile::t_chk->tic();
//    try{
//
//        if(state.MPS_L.size() + state.MPS_R.size() != state.get_length() )
//            throw std::runtime_error(fmt::format("Mismatch in MPS sizes: {} + {} != {}", state.MPS_L.size(), state.MPS_R.size(), state.get_length()));
//
//        if(state.ENV_L.size() + state.ENV_R.size() != state.get_length())
//            throw std::runtime_error(fmt::format("Mismatch in ENV sizes: {} + {} != {}", state.ENV_L.size(), state.ENV_R.size(),state.get_length()));
//
//        if(state.ENV2_L.size() + state.ENV2_R.size() != state.get_length())
//            throw std::runtime_error(fmt::format("Mismatch in ENV2 sizes: {} + {} != {}", state.ENV2_L.size(), state.ENV2_R.size(),state.get_length()));
//
//        if( state.ENV_L.back().sites != state.get_position())
//            throw std::runtime_error(fmt::format("Mismatch in ENV_L sites and position: {} != {}", state.ENV_L.back().sites, state.get_position()));
//
//        if(state.ENV_R.front().sites != state.get_length() - state.get_position() - 2)
//            throw std::runtime_error(fmt::format("Mismatch in ENV_R size+1 and length-position: {} != {}", state.ENV_R.front().sites, state.get_length() - state.get_position() - 2));
//
//        for(size_t pos = 0; pos < state.get_length(); pos++){
//            if(state.get_mps(pos).get_L().size() == 0)
//                throw std::runtime_error(fmt::format("Bond dimension is zero at position {}", pos));
//            if(pos != state.get_mps(pos).get_position())
//                throw std::runtime_error(fmt::format("Position mismatch {} != {}", pos, state.get_mps(pos).get_position()));
//            if(state.get_mps(pos).isCenter()){
//                if(state.get_mps(pos).get_LC().size() == 0)
//                    throw std::runtime_error(fmt::format("Center bond dimension is zero at position {}", pos));
//                if(pos != state.MPS_L.back().get_position())
//                    throw std::runtime_error(fmt::format("Center position mismatch {} != {}", pos, state.MPS_L.back().get_position()));
//                if(state.MPS_L.back().get_chiR() != state.MPS_L.back().get_LC().dimension(0) )
//                    throw std::runtime_error(fmt::format("Center and left dimension mismatch {} != {}", state.MPS_L.back().get_chiR(), state.MPS_L.back().get_LC().dimension(0)));
//                if(state.MPS_R.front().get_chiL() != state.MPS_L.back().get_LC().dimension(0) )
//                    throw std::runtime_error(fmt::format("Center and right dimension mismatch {} != {}", state.MPS_L.back().get_chiR(), state.MPS_L.back().get_LC().dimension(0)));
//            }
//        }
//
//        using VectorType = const Eigen::Matrix<class_state_finite::Scalar,Eigen::Dynamic,1>;
//        for (size_t pos = 1; pos < state.get_length(); pos++){
//
//            auto & mps_left = state.get_mps(pos - 1);
//            auto & mps_here = state.get_mps(pos);
//            auto & mpo_left = state.get_mpo(pos-1);
//            auto & mpo_here = state.get_mpo(pos);
//
//            // Check for validity
//            if(not Eigen::Map<VectorType>(mps_here.get_M().data(), mps_here.get_M().size()).allFinite()) {
//                std::cerr << "M: \n" << mps_here.get_M() << std::endl;
//                throw std::runtime_error(fmt::format("Inf's or nan's in MPS G @ pos {}", pos));
//            }
//
//            if(not Eigen::Map<VectorType>(mps_here.get_L().data(), mps_here.get_L().size()).allFinite()){
//                std::cerr << "L: \n" << mps_here.get_L() << std::endl;
//                throw std::runtime_error(fmt::format("Inf's or nan's in MPS L @ pos {}", pos));
//            }
//
//            // Check positions
//            if(mps_here.get_position() != mpo_here.get_position())
//                throw std::runtime_error(fmt::format("Mismatch in MPS and MPO positions @ pos {}: {} != {}", pos, mps_here.get_position(), mpo_here.get_position()));
//
//            if(mps_here.get_position() - mps_left.get_position() != 1)
//                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS positions @ pos {}: {} - {} != 1", pos, mps_here.get_position() , mps_left.get_position()));
//
//            if(mps_left.get_chiR() != mps_here.get_chiL())
//                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS dimensions @ pos {}: {} != {}", pos, mps_left.get_chiR() , mps_here.get_chiL()));
//
//
//            if(mpo_left.MPO().dimension(1) != mpo_here.MPO().dimension(0))
//                throw std::runtime_error(fmt::format("Mismatch in adjacent MPO dimensions @ pos {}: {} != {}", pos, mpo_left.MPO().dimension(1)  , mpo_here.MPO().dimension(0)));
//
//            if(mps_here.spin_dim() != mpo_here.MPO().dimension(2))
//                throw std::runtime_error(fmt::format("Mismatch in MPS and MPO spin dimensions @ pos {}: {} != {}", pos, mps_here.spin_dim() , mpo_here.MPO().dimension(2)));
//        }
//
//        {
//            //Check left side of the state
//            auto mps_it  = state.MPS_L.begin();
//            auto mps_nx  = state.MPS_L.begin();
//            auto env_it  = state.ENV_L.begin();
//            auto env_nx  = state.ENV_L.begin();
//            auto env2_it = state.ENV2_L.begin();
//            auto mpo_it  = state.MPO_L.begin();
//            std::advance(mps_nx,1);
//            std::advance(env_nx,1);
//            int i = 0;
//            while(
//                mps_it  != state.MPS_L.end() and
//                mps_nx  != state.MPS_L.end() and
//                env_it  != state.ENV_L.end() and
//                env_nx  != state.ENV_L.end() and
//                env2_it != state.ENV2_L.end() and
//                mpo_it  != state.MPO_L.end()
//                )
//            {
//                if(mps_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in MPS and ENV positions (left side) @ site {}: {} != {}", i, mps_it->get_position(), env_it->get_position()));
//
//                if(mps_it->get_position() != env_it->sites)
//                    throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV size (left side) @ site {}: {} != {}", i, mps_it->get_position(), env_it->sites));
//
//                if(mps_it->get_chiL() != env_it->block.dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in MPS and ENV dimensions (left side) @ site {}: {} != {}", i,mps_it->get_chiL() , env_it->block.dimension(0)));
//
//                if(env_nx->get_position() - env_it->get_position() != 1)
//                    throw std::runtime_error(fmt::format("Mismatch in adjacent ENV positions (left side) @ site {}: {} - {} != 1", i, env_nx->get_position(), env_it->get_position()));
//
//                if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV and MPO dimensions (left side) @ site {}: {} != {}", i,env_it->block.dimension(2), mpo_it->get()->MPO().dimension(0)));
//
//
//                if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (left side) @ site {}: {} != {}", i,env2_it->block.dimension(2), mpo_it->get()->MPO().dimension(0)));
//
//                if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (left side) @ site {}: {} != {}", i,env2_it->block.dimension(3) ,mpo_it->get()->MPO().dimension(0)));
//
//
//                if(env2_it->sites != env_it->sites)
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV sites (left side) @ site {}: {} != {}", i,env2_it->sites != env_it->sites));
//
//                if(env2_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV positions (left side) @ site {}: {} != {}", i,env2_it->get_position(), env_it->get_position()));
//
//                mps_it++;
//                mps_nx++;
//                env_it++;
//                env_nx++;
//                env2_it++;
//                mpo_it++;
//                i++;
//            }
//        }
//
//
//        {
//            auto mps_it  = state.MPS_R.rbegin();
//            auto mps_nx  = state.MPS_R.rbegin();
//            auto env_it  = state.ENV_R.rbegin();
//            auto env_nx  = state.ENV_R.rbegin();
//            auto env2_it = state.ENV2_R.rbegin();
//            auto mpo_it  = state.MPO_R.rbegin();
//            std::advance(mps_nx,1);
//            std::advance(env_nx,1);
//            auto i = state.get_length()-1;
//            while(
//                mps_it  != state.MPS_R.rend() and
//                mps_nx  != state.MPS_R.rend() and
//                env_it  != state.ENV_R.rend() and
//                env_nx  != state.ENV_R.rend() and
//                env2_it != state.ENV2_R.rend() and
//                mpo_it  != state.MPO_R.rend())
//            {
//                if(mps_it->get_chiL() != mps_nx->get_chiR())
//                    throw std::runtime_error(fmt::format("Mismatch in adjacent MPS dimensions (right side) @ site {}: {} != {}", i,  mps_nx->get_chiR(), mps_it->get_chiL()));
//
//                if(mps_it->get_position() - mps_nx->get_position() != 1 )
//                    throw std::runtime_error(fmt::format("Mismatch in adjacent MPS positions (right side) @ site {}: {} - {} != 1", i, mps_it->get_position(), mps_nx->get_position()));
//
//                if(mps_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in MPS and ENV positions (right side) @ site {}: {} != {}", i,  mps_it->get_position(), env_it->get_position()));
//
//                if(mps_it->get_position() != state.get_length() - (env_it->sites + 1))
//                    throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV size + 1 (right side) @ site {}: {} != {}", i,  mps_it->get_position(), state.get_length() - (env_it->sites + 1)));
//
//                if(mps_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV position (right side) @ site {}: {} != {}", i,  mps_it->get_position(), env_it->get_position()));
//
//                if(mps_it->get_chiR() != env_it->block.dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in MPS and ENV dimensions (right side) @ site {}: {} != {}", i,  mps_it->get_chiR(), env_it->block.dimension(0)));
//
//                if(env_it->get_position() - env_nx->get_position() != 1)
//                    throw std::runtime_error(fmt::format("Mismatch in adjacent ENV positions (right side) @ site {}: {} - {} != 1", i, env_it->get_position(), env_nx->get_position()));
//
//                if(env2_it->sites != env_it->sites)
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV sites (right side) @ site {}: {} != {}", i,  env2_it->sites, env_it->sites));
//
//                if(env2_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV positions (right side) @ site {}: {} != {}", i,  env2_it->get_position(), env_it->get_position()));
//
//
//                if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV and MPO dimensions (right side) @ site {}: {} != {}", i,  env_it->block.dimension(2), mpo_it->get()->MPO().dimension(1)));
//
//                if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (right side) @ site {}: {} != {}", i,  env2_it->block.dimension(2), mpo_it->get()->MPO().dimension(1)));
//
//                if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(1))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (right side) @ site {}: {} != {}", i,  env2_it->block.dimension(3), mpo_it->get()->MPO().dimension(1)));
//
//
//                mps_it++;
//                mps_nx++;
//                env_it++;
//                env_nx++;
//                env2_it++;
//                mpo_it++;
//                i--;
//            }
//        }
//
//
//        tools::log->trace("Checking norms");
//        auto norm_chain = tools::finite::measure::norm(state);
//        if(std::abs(norm_chain - 1.0) > settings::precision::max_norm_error) {
//            throw std::runtime_error(fmt::format("Norm of state too far from unity: {:.16f}",norm_chain));
//        }
//
//    }
//    catch(std::exception &ex){
//        throw std::runtime_error(fmt::format("Integrity check of MPS failed: {}", ex.what()));
//    }
//    tools::log->trace("MPS OK");
//    tools::common::profile::t_chk->toc();
}


void tools::finite::debug::check_integrity(const class_model_finite & model) {
//    if constexpr (not settings::debug) return;
//    tools::common::profile::t_chk->tic();
    model.assert_validity();


//
//    try{
//        for (auto &mpo : state.MPO_L){
//            auto site = mpo->get_position();
//            if (not mpo->all_mpo_parameters_have_been_set){throw std::runtime_error(fmt::format("All parameters have not been set on MPO_L site: {}" ,site));}
//        }
//        for (auto &mpo : state.MPO_R){
//            auto site = mpo->get_position();
//            if (not mpo->all_mpo_parameters_have_been_set){throw std::runtime_error(fmt::format("All parameters have not been set on MPO_R site: {}" ,site));}
//        }
//    }
//    catch(std::exception &ex){
//        throw std::runtime_error(fmt::format("Integrity check of MPO failed: {}", ex.what()));
//    }
//    tools::common::profile::t_chk->toc();

}


void tools::finite::debug::check_integrity(const class_edges_finite & edges) {
//    if constexpr (not settings::debug) return;
//    tools::log->trace("Checking integrity of edges");
//    tools::common::profile::t_chk->tic();
//    try{
//
//        if(not math::all_equal(state.get_length(),model.get_length(), edges.get_length()))
//            throw std::runtime_error(fmt::format("Lengths not all equal: state {} | model {} | edges {}", state.get_length(),model.get_length(), edges.get_length()));
//        for(const auto &ene : edges.ene)
//            if (not ene.)
//
//        if( state.ENV_L.back().sites != state.get_position())
//            throw std::runtime_error(fmt::format("Mismatch in ENV_L sites and position: {} != {}", state.ENV_L.back().sites, state.get_position()));
//
//        if(state.ENV_R.front().sites != state.get_length() - state.get_position() - 2)
//            throw std::runtime_error(fmt::format("Mismatch in ENV_R size+1 and length-position: {} != {}", state.ENV_R.front().sites, state.get_length() - state.get_position() - 2));
//
//        for(size_t pos = 0; pos < state.get_length(); pos++){
//            if(state.get_mps(pos).get_L().size() == 0)
//                throw std::runtime_error(fmt::format("Bond dimension is zero at position {}", pos));
//            if(pos != state.get_mps(pos).get_position())
//                throw std::runtime_error(fmt::format("Position mismatch {} != {}", pos, state.get_mps(pos).get_position()));
//            if(state.get_mps(pos).isCenter()){
//                if(state.get_mps(pos).get_LC().size() == 0)
//                    throw std::runtime_error(fmt::format("Center bond dimension is zero at position {}", pos));
//                if(pos != state.MPS_L.back().get_position())
//                    throw std::runtime_error(fmt::format("Center position mismatch {} != {}", pos, state.MPS_L.back().get_position()));
//                if(state.MPS_L.back().get_chiR() != state.MPS_L.back().get_LC().dimension(0) )
//                    throw std::runtime_error(fmt::format("Center and left dimension mismatch {} != {}", state.MPS_L.back().get_chiR(), state.MPS_L.back().get_LC().dimension(0)));
//                if(state.MPS_R.front().get_chiL() != state.MPS_L.back().get_LC().dimension(0) )
//                    throw std::runtime_error(fmt::format("Center and right dimension mismatch {} != {}", state.MPS_L.back().get_chiR(), state.MPS_L.back().get_LC().dimension(0)));
//            }
//        }
//
//        using VectorType = const Eigen::Matrix<class_state_finite::Scalar,Eigen::Dynamic,1>;
//        for (size_t pos = 1; pos < state.get_length(); pos++){
//
//            auto & mps_left = state.get_mps(pos - 1);
//            auto & mps_here = state.get_mps(pos);
//            auto & mpo_left = state.get_mpo(pos-1);
//            auto & mpo_here = state.get_mpo(pos);
//
//            // Check for validity
//            if(not Eigen::Map<VectorType>(mps_here.get_M().data(), mps_here.get_M().size()).allFinite()) {
//                std::cerr << "M: \n" << mps_here.get_M() << std::endl;
//                throw std::runtime_error(fmt::format("Inf's or nan's in MPS G @ pos {}", pos));
//            }
//
//            if(not Eigen::Map<VectorType>(mps_here.get_L().data(), mps_here.get_L().size()).allFinite()){
//                std::cerr << "L: \n" << mps_here.get_L() << std::endl;
//                throw std::runtime_error(fmt::format("Inf's or nan's in MPS L @ pos {}", pos));
//            }
//
//            // Check positions
//            if(mps_here.get_position() != mpo_here.get_position())
//                throw std::runtime_error(fmt::format("Mismatch in MPS and MPO positions @ pos {}: {} != {}", pos, mps_here.get_position(), mpo_here.get_position()));
//
//            if(mps_here.get_position() - mps_left.get_position() != 1)
//                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS positions @ pos {}: {} - {} != 1", pos, mps_here.get_position() , mps_left.get_position()));
//
//            if(mps_left.get_chiR() != mps_here.get_chiL())
//                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS dimensions @ pos {}: {} != {}", pos, mps_left.get_chiR() , mps_here.get_chiL()));
//
//
//            if(mpo_left.MPO().dimension(1) != mpo_here.MPO().dimension(0))
//                throw std::runtime_error(fmt::format("Mismatch in adjacent MPO dimensions @ pos {}: {} != {}", pos, mpo_left.MPO().dimension(1)  , mpo_here.MPO().dimension(0)));
//
//            if(mps_here.spin_dim() != mpo_here.MPO().dimension(2))
//                throw std::runtime_error(fmt::format("Mismatch in MPS and MPO spin dimensions @ pos {}: {} != {}", pos, mps_here.spin_dim() , mpo_here.MPO().dimension(2)));
//        }
//
//        {
//            //Check left side of the state
//            auto mps_it  = state.MPS_L.begin();
//            auto mps_nx  = state.MPS_L.begin();
//            auto env_it  = state.ENV_L.begin();
//            auto env_nx  = state.ENV_L.begin();
//            auto env2_it = state.ENV2_L.begin();
//            auto mpo_it  = state.MPO_L.begin();
//            std::advance(mps_nx,1);
//            std::advance(env_nx,1);
//            int i = 0;
//            while(
//                mps_it  != state.MPS_L.end() and
//                mps_nx  != state.MPS_L.end() and
//                env_it  != state.ENV_L.end() and
//                env_nx  != state.ENV_L.end() and
//                env2_it != state.ENV2_L.end() and
//                mpo_it  != state.MPO_L.end()
//                )
//            {
//                if(mps_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in MPS and ENV positions (left side) @ site {}: {} != {}", i, mps_it->get_position(), env_it->get_position()));
//
//                if(mps_it->get_position() != env_it->sites)
//                    throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV size (left side) @ site {}: {} != {}", i, mps_it->get_position(), env_it->sites));
//
//                if(mps_it->get_chiL() != env_it->block.dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in MPS and ENV dimensions (left side) @ site {}: {} != {}", i,mps_it->get_chiL() , env_it->block.dimension(0)));
//
//                if(env_nx->get_position() - env_it->get_position() != 1)
//                    throw std::runtime_error(fmt::format("Mismatch in adjacent ENV positions (left side) @ site {}: {} - {} != 1", i, env_nx->get_position(), env_it->get_position()));
//
//                if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV and MPO dimensions (left side) @ site {}: {} != {}", i,env_it->block.dimension(2), mpo_it->get()->MPO().dimension(0)));
//
//
//                if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (left side) @ site {}: {} != {}", i,env2_it->block.dimension(2), mpo_it->get()->MPO().dimension(0)));
//
//                if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (left side) @ site {}: {} != {}", i,env2_it->block.dimension(3) ,mpo_it->get()->MPO().dimension(0)));
//
//
//                if(env2_it->sites != env_it->sites)
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV sites (left side) @ site {}: {} != {}", i,env2_it->sites != env_it->sites));
//
//                if(env2_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV positions (left side) @ site {}: {} != {}", i,env2_it->get_position(), env_it->get_position()));
//
//                mps_it++;
//                mps_nx++;
//                env_it++;
//                env_nx++;
//                env2_it++;
//                mpo_it++;
//                i++;
//            }
//        }
//
//
//        {
//            auto mps_it  = state.MPS_R.rbegin();
//            auto mps_nx  = state.MPS_R.rbegin();
//            auto env_it  = state.ENV_R.rbegin();
//            auto env_nx  = state.ENV_R.rbegin();
//            auto env2_it = state.ENV2_R.rbegin();
//            auto mpo_it  = state.MPO_R.rbegin();
//            std::advance(mps_nx,1);
//            std::advance(env_nx,1);
//            auto i = state.get_length()-1;
//            while(
//                mps_it  != state.MPS_R.rend() and
//                mps_nx  != state.MPS_R.rend() and
//                env_it  != state.ENV_R.rend() and
//                env_nx  != state.ENV_R.rend() and
//                env2_it != state.ENV2_R.rend() and
//                mpo_it  != state.MPO_R.rend())
//            {
//                if(mps_it->get_chiL() != mps_nx->get_chiR())
//                    throw std::runtime_error(fmt::format("Mismatch in adjacent MPS dimensions (right side) @ site {}: {} != {}", i,  mps_nx->get_chiR(), mps_it->get_chiL()));
//
//                if(mps_it->get_position() - mps_nx->get_position() != 1 )
//                    throw std::runtime_error(fmt::format("Mismatch in adjacent MPS positions (right side) @ site {}: {} - {} != 1", i, mps_it->get_position(), mps_nx->get_position()));
//
//                if(mps_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in MPS and ENV positions (right side) @ site {}: {} != {}", i,  mps_it->get_position(), env_it->get_position()));
//
//                if(mps_it->get_position() != state.get_length() - (env_it->sites + 1))
//                    throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV size + 1 (right side) @ site {}: {} != {}", i,  mps_it->get_position(), state.get_length() - (env_it->sites + 1)));
//
//                if(mps_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV position (right side) @ site {}: {} != {}", i,  mps_it->get_position(), env_it->get_position()));
//
//                if(mps_it->get_chiR() != env_it->block.dimension(0))
//                    throw std::runtime_error(fmt::format("Mismatch in MPS and ENV dimensions (right side) @ site {}: {} != {}", i,  mps_it->get_chiR(), env_it->block.dimension(0)));
//
//                if(env_it->get_position() - env_nx->get_position() != 1)
//                    throw std::runtime_error(fmt::format("Mismatch in adjacent ENV positions (right side) @ site {}: {} - {} != 1", i, env_it->get_position(), env_nx->get_position()));
//
//                if(env2_it->sites != env_it->sites)
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV sites (right side) @ site {}: {} != {}", i,  env2_it->sites, env_it->sites));
//
//                if(env2_it->get_position() != env_it->get_position())
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV positions (right side) @ site {}: {} != {}", i,  env2_it->get_position(), env_it->get_position()));
//
//
//                if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV and MPO dimensions (right side) @ site {}: {} != {}", i,  env_it->block.dimension(2), mpo_it->get()->MPO().dimension(1)));
//
//                if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (right side) @ site {}: {} != {}", i,  env2_it->block.dimension(2), mpo_it->get()->MPO().dimension(1)));
//
//                if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(1))
//                    throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (right side) @ site {}: {} != {}", i,  env2_it->block.dimension(3), mpo_it->get()->MPO().dimension(1)));
//
//
//                mps_it++;
//                mps_nx++;
//                env_it++;
//                env_nx++;
//                env2_it++;
//                mpo_it++;
//                i--;
//            }
//        }
//
//
//        tools::log->trace("Checking norms");
//        auto norm_chain = tools::finite::measure::norm(state);
//        if(std::abs(norm_chain - 1.0) > settings::precision::max_norm_error) {
//            throw std::runtime_error(fmt::format("Norm of state too far from unity: {:.16f}",norm_chain));
//        }
//
//    }
//    catch(std::exception &ex){
//        throw std::runtime_error(fmt::format("Integrity check of MPS failed: {}", ex.what()));
//    }
//    tools::log->trace("MPS OK");
//    tools::common::profile::t_chk->toc();
}

void tools::finite::debug::check_integrity(const class_state_finite &state, const class_model_finite &model, const class_edges_finite &edges) {
//    if constexpr (not settings::debug) return;
//    tools::common::profile::t_chk->tic();
//    tools::log->trace("Integrity check...");
//    if(not math::all_equal(state.get_length(),model.get_length(), edges.get_length()))
//        throw std::runtime_error(fmt::format("Lengths not all equal: state {} | model {} | edges {}", state.get_length(),model.get_length(), edges.get_length()));
//
//    if(not math::all_equal(state.active_sites,model.active_sites, edges.active_sites))
//        throw std::runtime_error(fmt::format("Active sites not all equal: state {} | model {} | edges {}", state.active_sites,model.active_sites, edges.active_sites));
//
//    if(not math::all_equal(state.active_sites,model.active_sites, edges.active_sites))
//        throw std::runtime_error(fmt::format("Active sites not all equal: state {} | model {} | edges {}", state.active_sites,model.active_sites, edges.active_sites));
//
//
//
//    if(not math::all_equal(edges.ene.back().get_position(),model.get_length(), edges.get_length()))
//        throw std::runtime_error(fmt::format("Lengths not all equal: state {} | model {} | edges {}", state.get_length(),model.get_length(), edges.get_length()));
//
//    if( edges.ene.back().sites != state.get_position())
//        throw std::runtime_error(fmt::format("Mismatch in ENV_L sites and position: {} != {}", state.ENV_L.back().sites, state.get_position()));
//
//    if(state.ENV_R.front().sites != state.get_length() - state.get_position() - 2)
//        throw std::runtime_error(fmt::format("Mismatch in ENV_R size+1 and length-position: {} != {}", state.ENV_R.front().sites, state.get_length() - state.get_position() - 2));
//
//    for(size_t pos = 0; pos < state.get_length(); pos++){
//        if(state.get_mps(pos).get_L().size() == 0)
//            throw std::runtime_error(fmt::format("Bond dimension is zero at position {}", pos));
//        if(pos != state.get_mps(pos).get_position())
//            throw std::runtime_error(fmt::format("Position mismatch {} != {}", pos, state.get_mps(pos).get_position()));
//        if(state.get_mps(pos).isCenter()){
//            if(state.get_mps(pos).get_LC().size() == 0)
//                throw std::runtime_error(fmt::format("Center bond dimension is zero at position {}", pos));
//            if(pos != state.MPS_L.back().get_position())
//                throw std::runtime_error(fmt::format("Center position mismatch {} != {}", pos, state.MPS_L.back().get_position()));
//            if(state.MPS_L.back().get_chiR() != state.MPS_L.back().get_LC().dimension(0) )
//                throw std::runtime_error(fmt::format("Center and left dimension mismatch {} != {}", state.MPS_L.back().get_chiR(), state.MPS_L.back().get_LC().dimension(0)));
//            if(state.MPS_R.front().get_chiL() != state.MPS_L.back().get_LC().dimension(0) )
//                throw std::runtime_error(fmt::format("Center and right dimension mismatch {} != {}", state.MPS_L.back().get_chiR(), state.MPS_L.back().get_LC().dimension(0)));
//        }
//    }
//
//    using VectorType = const Eigen::Matrix<class_state_finite::Scalar,Eigen::Dynamic,1>;
//    for (size_t pos = 1; pos < state.get_length(); pos++){
//
//        auto & mps_left = state.get_mps(pos - 1);
//        auto & mps_here = state.get_mps(pos);
//        auto & mpo_left = state.get_mpo(pos-1);
//        auto & mpo_here = state.get_mpo(pos);
//
//        // Check for validity
//        if(not Eigen::Map<VectorType>(mps_here.get_M().data(), mps_here.get_M().size()).allFinite()) {
//            std::cerr << "M: \n" << mps_here.get_M() << std::endl;
//            throw std::runtime_error(fmt::format("Inf's or nan's in MPS G @ pos {}", pos));
//        }
//
//        if(not Eigen::Map<VectorType>(mps_here.get_L().data(), mps_here.get_L().size()).allFinite()){
//            std::cerr << "L: \n" << mps_here.get_L() << std::endl;
//            throw std::runtime_error(fmt::format("Inf's or nan's in MPS L @ pos {}", pos));
//        }
//
//        // Check positions
//        if(mps_here.get_position() != mpo_here.get_position())
//            throw std::runtime_error(fmt::format("Mismatch in MPS and MPO positions @ pos {}: {} != {}", pos, mps_here.get_position(), mpo_here.get_position()));
//
//        if(mps_here.get_position() - mps_left.get_position() != 1)
//            throw std::runtime_error(fmt::format("Mismatch in adjacent MPS positions @ pos {}: {} - {} != 1", pos, mps_here.get_position() , mps_left.get_position()));
//
//        if(mps_left.get_chiR() != mps_here.get_chiL())
//            throw std::runtime_error(fmt::format("Mismatch in adjacent MPS dimensions @ pos {}: {} != {}", pos, mps_left.get_chiR() , mps_here.get_chiL()));
//
//
//        if(mpo_left.MPO().dimension(1) != mpo_here.MPO().dimension(0))
//            throw std::runtime_error(fmt::format("Mismatch in adjacent MPO dimensions @ pos {}: {} != {}", pos, mpo_left.MPO().dimension(1)  , mpo_here.MPO().dimension(0)));
//
//        if(mps_here.spin_dim() != mpo_here.MPO().dimension(2))
//            throw std::runtime_error(fmt::format("Mismatch in MPS and MPO spin dimensions @ pos {}: {} != {}", pos, mps_here.spin_dim() , mpo_here.MPO().dimension(2)));
//    }
//
//    {
//        //Check left side of the state
//        auto mps_it  = state.MPS_L.begin();
//        auto mps_nx  = state.MPS_L.begin();
//        auto env_it  = state.ENV_L.begin();
//        auto env_nx  = state.ENV_L.begin();
//        auto env2_it = state.ENV2_L.begin();
//        auto mpo_it  = state.MPO_L.begin();
//        std::advance(mps_nx,1);
//        std::advance(env_nx,1);
//        int i = 0;
//        while(
//            mps_it  != state.MPS_L.end() and
//            mps_nx  != state.MPS_L.end() and
//            env_it  != state.ENV_L.end() and
//            env_nx  != state.ENV_L.end() and
//            env2_it != state.ENV2_L.end() and
//            mpo_it  != state.MPO_L.end()
//            )
//        {
//            if(mps_it->get_position() != env_it->get_position())
//                throw std::runtime_error(fmt::format("Mismatch in MPS and ENV positions (left side) @ site {}: {} != {}", i, mps_it->get_position(), env_it->get_position()));
//
//            if(mps_it->get_position() != env_it->sites)
//                throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV size (left side) @ site {}: {} != {}", i, mps_it->get_position(), env_it->sites));
//
//            if(mps_it->get_chiL() != env_it->block.dimension(0))
//                throw std::runtime_error(fmt::format("Mismatch in MPS and ENV dimensions (left side) @ site {}: {} != {}", i,mps_it->get_chiL() , env_it->block.dimension(0)));
//
//            if(env_nx->get_position() - env_it->get_position() != 1)
//                throw std::runtime_error(fmt::format("Mismatch in adjacent ENV positions (left side) @ site {}: {} - {} != 1", i, env_nx->get_position(), env_it->get_position()));
//
//            if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
//                throw std::runtime_error(fmt::format("Mismatch in ENV and MPO dimensions (left side) @ site {}: {} != {}", i,env_it->block.dimension(2), mpo_it->get()->MPO().dimension(0)));
//
//
//            if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(0))
//                throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (left side) @ site {}: {} != {}", i,env2_it->block.dimension(2), mpo_it->get()->MPO().dimension(0)));
//
//            if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(0))
//                throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (left side) @ site {}: {} != {}", i,env2_it->block.dimension(3) ,mpo_it->get()->MPO().dimension(0)));
//
//
//            if(env2_it->sites != env_it->sites)
//                throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV sites (left side) @ site {}: {} != {}", i,env2_it->sites != env_it->sites));
//
//            if(env2_it->get_position() != env_it->get_position())
//                throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV positions (left side) @ site {}: {} != {}", i,env2_it->get_position(), env_it->get_position()));
//
//            mps_it++;
//            mps_nx++;
//            env_it++;
//            env_nx++;
//            env2_it++;
//            mpo_it++;
//            i++;
//        }
//    }
//
//
//    {
//        auto mps_it  = state.MPS_R.rbegin();
//        auto mps_nx  = state.MPS_R.rbegin();
//        auto env_it  = state.ENV_R.rbegin();
//        auto env_nx  = state.ENV_R.rbegin();
//        auto env2_it = state.ENV2_R.rbegin();
//        auto mpo_it  = state.MPO_R.rbegin();
//        std::advance(mps_nx,1);
//        std::advance(env_nx,1);
//        auto i = state.get_length()-1;
//        while(
//            mps_it  != state.MPS_R.rend() and
//            mps_nx  != state.MPS_R.rend() and
//            env_it  != state.ENV_R.rend() and
//            env_nx  != state.ENV_R.rend() and
//            env2_it != state.ENV2_R.rend() and
//            mpo_it  != state.MPO_R.rend())
//        {
//            if(mps_it->get_chiL() != mps_nx->get_chiR())
//                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS dimensions (right side) @ site {}: {} != {}", i,  mps_nx->get_chiR(), mps_it->get_chiL()));
//
//            if(mps_it->get_position() - mps_nx->get_position() != 1 )
//                throw std::runtime_error(fmt::format("Mismatch in adjacent MPS positions (right side) @ site {}: {} - {} != 1", i, mps_it->get_position(), mps_nx->get_position()));
//
//            if(mps_it->get_position() != env_it->get_position())
//                throw std::runtime_error(fmt::format("Mismatch in MPS and ENV positions (right side) @ site {}: {} != {}", i,  mps_it->get_position(), env_it->get_position()));
//
//            if(mps_it->get_position() != state.get_length() - (env_it->sites + 1))
//                throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV size + 1 (right side) @ site {}: {} != {}", i,  mps_it->get_position(), state.get_length() - (env_it->sites + 1)));
//
//            if(mps_it->get_position() != env_it->get_position())
//                throw std::runtime_error(fmt::format("Mismatch in MPS position and ENV position (right side) @ site {}: {} != {}", i,  mps_it->get_position(), env_it->get_position()));
//
//            if(mps_it->get_chiR() != env_it->block.dimension(0))
//                throw std::runtime_error(fmt::format("Mismatch in MPS and ENV dimensions (right side) @ site {}: {} != {}", i,  mps_it->get_chiR(), env_it->block.dimension(0)));
//
//            if(env_it->get_position() - env_nx->get_position() != 1)
//                throw std::runtime_error(fmt::format("Mismatch in adjacent ENV positions (right side) @ site {}: {} - {} != 1", i, env_it->get_position(), env_nx->get_position()));
//
//            if(env2_it->sites != env_it->sites)
//                throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV sites (right side) @ site {}: {} != {}", i,  env2_it->sites, env_it->sites));
//
//            if(env2_it->get_position() != env_it->get_position())
//                throw std::runtime_error(fmt::format("Mismatch in ENV2 position and ENV positions (right side) @ site {}: {} != {}", i,  env2_it->get_position(), env_it->get_position()));
//
//
//            if(env_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
//                throw std::runtime_error(fmt::format("Mismatch in ENV and MPO dimensions (right side) @ site {}: {} != {}", i,  env_it->block.dimension(2), mpo_it->get()->MPO().dimension(1)));
//
//            if(env2_it->block.dimension(2) != mpo_it->get()->MPO().dimension(1))
//                throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (right side) @ site {}: {} != {}", i,  env2_it->block.dimension(2), mpo_it->get()->MPO().dimension(1)));
//
//            if(env2_it->block.dimension(3) != mpo_it->get()->MPO().dimension(1))
//                throw std::runtime_error(fmt::format("Mismatch in ENV2 and MPO dimensions (right side) @ site {}: {} != {}", i,  env2_it->block.dimension(3), mpo_it->get()->MPO().dimension(1)));
//
//
//            mps_it++;
//            mps_nx++;
//            env_it++;
//            env_nx++;
//            env2_it++;
//            mpo_it++;
//            i--;
//        }
//    }
//
//
//    tools::log->trace("Checking norms");
//    auto norm_chain = tools::finite::measure::norm(state);
//    if(std::abs(norm_chain - 1.0) > settings::precision::max_norm_error) {
//        throw std::runtime_error(fmt::format("Norm of state too far from unity: {:.16f}",norm_chain));
//    }
//
//
//    tools::log->trace("Integrity check... OK");

}




void tools::finite::debug::print_parity_properties(const class_state_finite &state) {
    tools::log->debug("Printing parity properties");

    tools::log->debug("\tComputing spin components");
    const auto sx = tools::finite::measure::spin_component(state,qm::spinOneHalf::sx);
    const auto sy = tools::finite::measure::spin_component(state,qm::spinOneHalf::sy);
    const auto sz = tools::finite::measure::spin_component(state,qm::spinOneHalf::sz);
    tools::log->debug("\t<psi | sx | psi>                = {:0.16f}", sx);
    tools::log->debug("\t<psi | sy | psi>                = {:0.16f}", sy);
    tools::log->debug("\t<psi | sz | psi>                = {:0.16f}", sz);

    tools::log->debug("\tComputing parity projected states");
    auto state_sx = tools::finite::ops::get_projection_to_closest_parity_sector(state, "x");
    auto state_sy = tools::finite::ops::get_projection_to_closest_parity_sector(state, "x");
    auto state_sz = tools::finite::ops::get_projection_to_closest_parity_sector(state, "x");



    tools::log->debug("\tMore spin components");
    tools::log->debug("\t<psi_sx | sx | psi_sx>      = {:0.16f}", tools::finite::measure::spin_component(state_sx,qm::spinOneHalf::sx));
    tools::log->debug("\t<psi_sy | sy | psi_sy>      = {:0.16f}", tools::finite::measure::spin_component(state_sy,qm::spinOneHalf::sy));
    tools::log->debug("\t<psi_sz | sz | psi_sz>      = {:0.16f}", tools::finite::measure::spin_component(state_sz,qm::spinOneHalf::sz));



    tools::log->debug("\tNormalization check");
    tools::log->debug("\tComputing overlaps");
    tools::log->debug("\t<psi_sx|psi_sx>             = {:0.16f}", tools::finite::ops::overlap(state_sx,state_sx));
    tools::log->debug("\t<psi_sy|psi_sy>             = {:0.16f}", tools::finite::ops::overlap(state_sy,state_sy));
    tools::log->debug("\t<psi_sz|psi_sz>             = {:0.16f}", tools::finite::ops::overlap(state_sz,state_sz));


    tools::log->debug("\tOverlaps with original state");
    tools::log->debug("\t<psi|psi_sx>                  = {:0.16f}", tools::finite::ops::overlap(state,state_sx));
    tools::log->debug("\t<psi|psi_sy>                  = {:0.16f}", tools::finite::ops::overlap(state,state_sy));
    tools::log->debug("\t<psi|psi_sz>                  = {:0.16f}", tools::finite::ops::overlap(state,state_sz));


    tools::log->debug("\tOverlaps between different direction sectors");
    tools::log->debug("\t<psi_sx|psi_sy>             = {:0.16f}" ,tools::finite::ops::overlap(state_sx,state_sy));
    tools::log->debug("\t<psi_sx|psi_sz>             = {:0.16f}" ,tools::finite::ops::overlap(state_sx,state_sz));
}



void tools::finite::debug::check_normalization_routine(const class_state_finite &state){
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
