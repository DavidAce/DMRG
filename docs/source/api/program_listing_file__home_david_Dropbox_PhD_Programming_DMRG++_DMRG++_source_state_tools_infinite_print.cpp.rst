
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_tools_infinite_print.cpp:

Program Listing for File print.cpp
==================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_tools_infinite_print.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/state/tools/infinite/print.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-02-20.
   //
   
   #include <state/tools/nmspc_tools.h>
   #include <state/class_infinite_state.h>
   #include <state/class_mps_2site.h>
   #include <model/class_model_base.h>
   
   void tools::infinite::print::print_hamiltonians(const class_infinite_state & state){
       state.HA->print_parameter_names();
       state.HA->print_parameter_values();
       state.HB->print_parameter_values();
   
   }
   
   void tools::infinite::print::print_state(const class_infinite_state & state){
       using namespace Textra;
       std::cout << std::setprecision(10);
       std::cout << "State length              : "    << state.get_length()   << std::endl;
       std::cout << "State position            : "    << state.get_position() << std::endl;
   }
