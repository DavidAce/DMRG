
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_infinite_env.cpp:

Program Listing for File env.cpp
================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_infinite_env.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/infinite/env.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-10-24.
   //
   
   #include <tools/nmspc_tools.h>
   #include <state/class_state_infinite.h>
   #include <state/class_environment.h>
   #include <state/class_mps_2site.h>
   #include <state/class_mps_site.h>
   #include <model/class_model_base.h>
   void tools::infinite::env::initialize  (class_state_infinite & state){
       state.Lblock   = std::make_unique<class_environment>("L",*state.MPS->MPS_A, *state.HA);
       state.Rblock   = std::make_unique<class_environment>("R",*state.MPS->MPS_B, *state.HB);
       state.Lblock2  = std::make_unique<class_environment_var>("L",*state.MPS->MPS_A, *state.HA);
       state.Rblock2  = std::make_unique<class_environment_var>("R",*state.MPS->MPS_B, *state.HB);
   }
