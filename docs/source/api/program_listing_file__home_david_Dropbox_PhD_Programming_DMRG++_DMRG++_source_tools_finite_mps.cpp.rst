
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_mps.cpp:

Program Listing for File mps.cpp
================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_mps.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/mps.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-01-29.
   //
   
   
   #include <tools/nmspc_tools.h>
   #include <state/class_state_finite.h>
   #include <state/class_mps_2site.h>
   #include <state/class_environment.h>
   #include <general/nmspc_random_numbers.h>
   #include <general/nmspc_quantum_mechanics.h>
   #include <simulation/nmspc_settings.h>
   
   
   void tools::finite::mps::initialize(class_state_finite &state, const size_t length){
       log->info("Initializing mps");
       using Scalar = class_state_finite::Scalar;
       //Generate MPS
       Eigen::Tensor<Scalar,3> A;
       Eigen::Tensor<Scalar,3> B;
       Eigen::Tensor<Scalar,1> L = {1};
       size_t pos = 0;
       state.MPS_L.emplace_back(class_mps_site(A, L, pos++));
       state.MPS_L.back().set_LC(L);
       while(true){
           state.MPS_R.emplace_back(class_mps_site(B, L, pos++));
           if(state.MPS_L.size() + state.MPS_R.size() >= length){break;}
       }
       state.site_update_tags = std::vector<bool>(length,false);
   }
   
   
   
   void tools::finite::mps::randomize(class_state_finite &state, const std::string &parity_sector, int seed_state, bool use_pauli_eigenstates, bool enumeration)
   {
       tools::log->trace("Randomizing mps");
       state.clear_measurements();
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
   
   
   void tools::finite::mps::rebuild_environments(class_state_finite &state){
       tools::log->trace("Rebuilding environments");
       if (state.MPS_L.size() != state.MPO_L.size())
           throw std::runtime_error(fmt::format("Size mismatch in MPSL and MPOL: {} != {}",state.MPS_L.size(), state.MPO_L.size()));
       if (state.MPS_R.size() != state.MPO_R.size())
           throw std::runtime_error(fmt::format("Size mismatch in MPSR and MPOR: {} != {}",state.MPS_R.size(), state.MPO_R.size()));
       // Generate new environments
   
       {
           state.ENV_L.clear();
           state.ENV2_L.clear();
   
   
           auto mpsL_it   = state.MPS_L.begin();
           auto mpoL_it   = state.MPO_L.begin();
           auto ENV_L     = class_environment    ("L",*mpsL_it, *mpoL_it->get()); //Initialized envs
           auto ENV2_L    = class_environment_var("L",*mpsL_it, *mpoL_it->get()); //Initialized envs
           while(mpsL_it != state.MPS_L.end() and mpoL_it != state.MPO_L.end()) {
               state.ENV_L .emplace_back(ENV_L );
               state.ENV2_L.emplace_back(ENV2_L);
               if (mpsL_it->get_position() != state.ENV_L.back().get_position())
                   throw std::runtime_error(fmt::format("Size mismatch in MPSL and ENVL: {} != {}",mpsL_it->get_position(), state.ENV_L.back().get_position()));
               if (mpsL_it->get_chiL() != state.ENV_L.back().block.dimension(0))
                   throw std::runtime_error(fmt::format("Size mismatch in MPSL and ENVL dimensions {} != {}",mpsL_it->get_chiL(), state.ENV_L.back().block.dimension(2)));
   
               ENV_L  = ENV_L .enlarge(*mpsL_it,*mpoL_it->get());
               ENV2_L = ENV2_L.enlarge(*mpsL_it,*mpoL_it->get());
               mpsL_it++;
               mpoL_it++;
           }
       }
   
       {
   
           state.ENV_R.clear();
           state.ENV2_R.clear();
   
           auto mpsR_it   = state.MPS_R.rbegin();
           auto mpoR_it   = state.MPO_R.rbegin();
           auto ENV_R  = class_environment    ("R",*mpsR_it,*mpoR_it->get());
           auto ENV2_R = class_environment_var("R",*mpsR_it,*mpoR_it->get());
           while(mpsR_it != state.MPS_R.rend() and mpoR_it != state.MPO_R.rend()){
               state.ENV_R .emplace_front(ENV_R );
               state.ENV2_R.emplace_front(ENV2_R);
               if (mpsR_it->get_position() != state.ENV_R.front().get_position())
                   throw std::runtime_error(fmt::format("Size mismatch in MPSR and ENVR: {} != {}",mpsR_it->get_position(), state.ENV_R.front().get_position()));
               if (mpsR_it->get_chiR() != state.ENV_R.front().block.dimension(0))
                   throw std::runtime_error(fmt::format("Size mismatch in MPSR and ENVR dimensions {} != {}",mpsR_it->get_chiR() , state.ENV_R.front().block.dimension(2)));
               ENV_R  = ENV_R .enlarge(*mpsR_it,*mpoR_it->get());
               ENV2_R = ENV2_R.enlarge(*mpsR_it,*mpoR_it->get());
               mpsR_it++;
               mpoR_it++;
   
           }
       }
   }
   
   
   
   void tools::finite::mps::project_to_closest_parity_sector   (class_state_finite & state, std::string parity_sector){
       state = tools::finite::ops::get_projection_to_closest_parity_sector(state, parity_sector);
   }
