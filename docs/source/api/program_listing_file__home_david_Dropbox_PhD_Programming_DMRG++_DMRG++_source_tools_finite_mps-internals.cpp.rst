
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_mps-internals.cpp:

Program Listing for File mps-internals.cpp
==========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_mps-internals.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/mps-internals.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-08-12.
   //
   
   #include <tools/nmspc_tools.h>
   #include <state/class_finite_state.h>
   #include <state/class_mps_2site.h>
   #include <state/class_environment.h>
   #include <general/nmspc_random_numbers.h>
   #include <general/nmspc_quantum_mechanics.h>
   #include <simulation/nmspc_settings.h>
   #include <bitset>
   
   
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
       if (sign == 0 ) return rn::uniform_integer_1();
       throw std::runtime_error("Invalid sign in get_elem.");
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
   
   void tools::finite::mps::internals::set_product_state_in_parity_sector_from_bitset(class_finite_state & state, const std::string &parity_sector, const int seed_state){
       if (seed_state < 0){
           throw std::runtime_error(fmt::format("Can't set sector from bitset with negative seed_state: {}", seed_state));
       }
       std::vector<std::string> ok_parity_sectors = {"x","+x","-x","y","+y","-y", "z","+z","-z"};
       bool parity_sector_is_defined = std::find(ok_parity_sectors.begin(), ok_parity_sectors.end(), parity_sector) != ok_parity_sectors.end();
       if (not parity_sector_is_defined)
           throw std::logic_error(fmt::format("Can't use seed_state as enumeration when parity_sector is not well defined. Got: {}", parity_sector));
   
   
       constexpr int maxbits = 128;
       if (maxbits > state.get_length()) throw std::range_error("Max supported state length for bitset is 128");
       std::bitset<maxbits> bs (seed_state);
   
       std::string axis = get_axis(parity_sector);
       int sector       = get_sign(parity_sector);
       Eigen::Tensor<Scalar,1> L (1);
       L.setConstant(1);
       int carry_sign = 1;
       for (auto &mpsL : state.MPS_L ){
           int sign = 2*bs[mpsL.get_position()] - 1;
           carry_sign *= sign;
           auto G = Textra::Matrix_to_Tensor(get_eigvec(axis,sign).normalized(),2,1,1);
           mpsL.set_mps(G,L);
       }
       state.MPS_C = L;
       for (auto &mpsR : state.MPS_R ){
           int sign = 2*bs[mpsR.get_position()] - 1;
           carry_sign *= sign;
           auto G = Textra::Matrix_to_Tensor(get_eigvec(axis,sign).normalized(),2,1,1);
           mpsR.set_mps(G,L);
       }
   
       if(sector * carry_sign == -1){
           //Flip the last spin to get the correct total sign.
           auto &mpsR = state.MPS_R.back();
           int sign = 2*bs[mpsR.get_position()] - 1;
           sign *= -1;
           auto G = Textra::Matrix_to_Tensor(get_eigvec(axis,sign).normalized(),2,1,1);
           mpsR.set_mps(G,L);
       }
   }
   
   
   
   void tools::finite::mps::internals::set_product_state_in_parity_sector_randomly(class_finite_state & state, const std::string &parity_sector){
   
       Eigen::Tensor<Scalar,1> L (1);
       std::string axis = get_axis(parity_sector);
       int sector       = get_sign(parity_sector);
       int carry_sign = 1;
       int last_sign  = 1;
   
       L.setConstant(1);
       for (auto &mpsL : state.MPS_L ){
           int sign = 2*rn::uniform_integer_1()-1;
           carry_sign *= sign;
           auto G = Textra::Matrix_to_Tensor(get_eigvec(axis,sign).normalized(), 2, 1, 1);
           mpsL.set_mps(G,L);
       }
       state.MPS_C = L;
       for (auto &mpsR : state.MPS_R ){
           int sign = 2*rn::uniform_integer_1()-1;
           carry_sign *= sign;
           last_sign = sign;
           auto G = Textra::Matrix_to_Tensor(get_eigvec(axis, sign).normalized(), 2, 1, 1);
           mpsR.set_mps(G,L);
       }
   
       if(sector * carry_sign == -1){
           //Flip the last spin to get the correct total sign.
           auto &mpsR = state.MPS_R.back();
           int sign = -last_sign;
           sign *= -1;
           auto G = Textra::Matrix_to_Tensor(get_eigvec(axis, sign).normalized(), 2, 1, 1);
           mpsR.set_mps(G,L);
       }
   }
   
   
   
   void tools::finite::mps::internals::set_product_state_randomly(class_finite_state & state,const std::string &parity_sector,bool use_pauli_eigenstates){
       std::vector<std::string> ok_parity_sectors = {"x","+x","-x","y","+y","-y", "z","+z","-z"};
       bool parity_sector_is_defined = std::find(ok_parity_sectors.begin(), ok_parity_sectors.end(), parity_sector) != ok_parity_sectors.end();
       if (parity_sector_is_defined and use_pauli_eigenstates){
           // Case a)
           set_product_state_in_parity_sector_randomly(state,parity_sector);
       }
       else if (parity_sector_is_defined and not use_pauli_eigenstates){
           set_product_state_randomly(state,"random",false);
           state = tools::finite::ops::get_projection_to_closest_parity_sector(state,parity_sector,false);
       }
       else if (parity_sector == "randomAxis") {
           std::vector<std::string> possibilities = {"x", "y", "z"};
           std::string chosen_axis = possibilities[rn::uniform_integer(0, 2)];
           set_product_state_in_parity_sector_randomly(state, chosen_axis);
       }else if (parity_sector == "random") {
           Eigen::Tensor<Scalar,1> L (1);
           L.setConstant(1);
           for (auto &mpsL : state.MPS_L ){
               auto G = Textra::Matrix_to_Tensor(Eigen::VectorXcd::Random(2).normalized(),2,1,1);
               mpsL.set_mps(G,L);
           }
           state.MPS_C = L;
           for (auto &mpsR : state.MPS_R ){
               auto G = Textra::Matrix_to_Tensor(Eigen::VectorXcd::Random(2).normalized(),2,1,1);
               mpsR.set_mps(G,L);
           }
       }else if (parity_sector == "none"){
           return;
       }else{
           throw std::runtime_error(fmt::format(R"(Wrong pauli string. Expected one of (+-) "x","y","z", "randomAxis", "random" or "none". Got: )" + parity_sector));
       }
   
   }
