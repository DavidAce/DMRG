
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_tools_infinite_mps.cpp:

Program Listing for File mps.cpp
================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_tools_infinite_mps.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/state/tools/infinite/mps.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-06-25.
   //
   
   #include <state/tools/nmspc_tools.h>
   #include <state/class_infinite_state.h>
   
   
   class_infinite_state tools::infinite::mps::set_random_state(const class_infinite_state & state, std::string parity){
       throw std::runtime_error("You need to implement set random state for infinite state");
       return state;
   }
