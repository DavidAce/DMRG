
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_print.cpp:

Program Listing for File print.cpp
==================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_print.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/print.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-01-29.
   //
   
   
   #include <tools/nmspc_tools.h>
   #include <state/class_state_finite.h>
   #include <state/class_environment.h>
   #include <general/nmspc_tensor_extra.h>
   #include <string>
   #include <iomanip>
   
   
   using Scalar         = std::complex<double>;
   
   void tools::finite::print::print_full_state(const class_state_finite &state) {
       
       for (auto & mps : state.MPS_L){
           std::cout << "MPS " << mps.get_position() << "  :\n";
           std::cout << "  L:\n"<< mps.get_L() << '\n';
           std::cout << "  M:\n"<< mps.get_M() << '\n';
       }
       std::cout << "  LC:\n"<< state.MPS_L.back().get_LC() << '\n';
       for (auto & mps : state.MPS_R){
           std::cout << "MPS " << mps.get_position() << "  :\n";
           std::cout << "  M:\n"<< mps.get_M() << '\n';
           std::cout << "  L:\n"<< mps.get_L() << '\n';
       }
   }
   
   
   
   void tools::finite::print::print_state(const class_state_finite &state){
       using namespace Textra;
       auto & MPS_L  = state.MPS_L;
       auto & MPS_R  = state.MPS_R;
       auto & MPS_C  = state.midchain_bond();
       auto & ENV_L  = state.ENV_L;
       auto & ENV_R  = state.ENV_R;
   
       int i = 0;
       std::cout << std::setprecision(10);
       std::cout << "State length              : " << state.get_length() << std::endl;
       std::cout << "State position            : "    << state.get_position() << std::endl;
       std::cout << "Environment L size        : "    << ENV_L.size() << std::endl;
       std::cout << "Environment R size        : "    << ENV_R.size() << std::endl;
   
       auto envitL = ENV_L.begin();
       std::cout << std::left;
       for(auto &it : MPS_L){
           std::cout << "L[" << std::setw(3) << i  <<  "]: " << it.get_L().dimensions()<< std::setw(5) << "   "
                     << "M[" << std::setw(3) << i  <<  "]: " << it.get_M().dimensions()<< std::setw(5) << " pos: " << it.get_position() << "   ";
           if (envitL != ENV_L.end()){std::cout << " ENV_" << envitL->side << ": " << envitL->block.dimensions() << " pos: " << envitL->get_position() << "   " << " env spins: " << envitL++->sites << " ";}
           if(&it == &MPS_L.back()){
               std::cout << " <--- Position A";
           }
           std::cout << std::endl;
           i++;
       }
       std::cout << "L[" << std::setw(3) << '*'  <<  "]: " << MPS_C.dimensions() << std::setw(80) << std::right << "<--- Center" << std::left << std::endl;
       auto envitR = ENV_R.begin();
       for(auto &it : MPS_R){
           std::cout << "M[" << std::setw(3) << i  <<  "]: " << it.get_M().dimensions() << std::setw(5) << "  "
                     << "L[" << std::setw(3) << i  <<  "]: " << it.get_L().dimensions() << std::setw(5) << " pos: " << it.get_position() << "  ";
           if (envitR != ENV_R.end()){std::cout << " ENV_" << envitR->side << ": " << envitR->block.dimensions() << " pos: " << envitR->get_position()  << "   "<< " env spins: " << envitR++->sites << " ";}
           if(&it == &MPS_R.front()){
               std::cout << " <--- Position B" ;
           }
           std::cout << std::endl;
           i++;
       }
   
   }
   
   
   void tools::finite::print::print_state_compact(const class_state_finite &state){
       using namespace Textra;
       auto & MPS_L  = state.MPS_L;
       auto & MPS_R  = state.MPS_R;
       auto & MPS_C  = state.midchain_bond();
       auto & ENV_L  = state.ENV_L;
       auto & ENV_R  = state.ENV_R;
   
       std::cout << std::setprecision(10);
   
       std::cout << "State length              : "    << state.get_length() << std::endl;
       std::cout << "State position            : "    << state.get_position() << std::endl;
       std::cout << "Environment L size        : "    << ENV_L.size() << std::endl;
       std::cout << "Environment R size        : "    << ENV_R.size() << std::endl;
       if(!ENV_L.empty()){std::cout << "ENV_L[" <<std::setw(3) << ENV_L.size()-1 << "]: " << ENV_L.back().block.dimensions() << " Particles: " << ENV_L.back().sites << "  <--- Also current environment L" << std::endl;}
       if(!MPS_L.empty()){std::cout << "MPS_L[" <<std::setw(3) << MPS_L.size()-1 << "]: " << MPS_L.back().get_M().dimensions() <<  "   <--- Also current M" << std::endl;}
       std::cout << "L[" << std::setw(3) << '*'  <<  "]: " << MPS_C.dimensions() << "                    <--- Center" << std::endl;
       if(!MPS_R.empty()){std::cout << "MPS_R[" <<std::setw(3) << MPS_R.size()-1 << "]: " << MPS_R.front().get_M().dimensions() << "   <--- Also current M" << std::endl;}
       if(!ENV_R.empty()){std::cout << "ENV_R[" <<std::setw(3) << ENV_R.size()-1 << "]: " << ENV_R.front().block.dimensions() << " Particles: " << ENV_R.front().sites << " <--- Also current environment R"  << std::endl;}
   }
   
   
   
   
   void tools::finite::print::print_hamiltonians(const class_state_finite &state) {
       auto & MPO_L  = state.MPO_L;
       auto & MPO_R  = state.MPO_R;
       if (MPO_L.empty()) throw std::runtime_error("MPO_L is empty. Can't print hamiltonian");
       if (MPO_R.empty()) throw std::runtime_error("MPO_R is empty. Can't print hamiltonian");
   
       MPO_L.begin()->get()->print_parameter_names();
       for(auto &it : MPO_L){
           it->print_parameter_values();
       }
       for(auto &it : MPO_R){
           it->print_parameter_values();
       }
   }
