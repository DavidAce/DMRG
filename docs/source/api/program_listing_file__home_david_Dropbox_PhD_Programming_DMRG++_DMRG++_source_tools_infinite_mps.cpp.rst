
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_infinite_mps.cpp:

Program Listing for File mps.cpp
================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_infinite_mps.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/infinite/mps.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-06-25.
   //
   
   #include <tools/nmspc_tools.h>
   #include <state/class_state_infinite.h>
   #include <state/class_mps_2site.h>
   #include <state/class_mps_site.h>
   
   
   void tools::infinite::mps::initialize  (class_state_infinite & state, std::string model_type_str){
       tools::log->trace("Constructing 2site MPS");
       long spin_dimension = 2;
       if      (model_type_str == "tf_ising")             spin_dimension = 2;
       else if (model_type_str == "tf_nn_ising")          spin_dimension = 2;
       else if (model_type_str == "selfdual_tf_rf_ising") spin_dimension = 2;
   
       Eigen::Tensor<Scalar,3> M(spin_dimension,1,1);
       Eigen::Tensor<Scalar,1> L(1);
   
       // Default is a product state, spins pointing up in z.
       M.setZero();
       M(0,0,0) = 1;
       L.setConstant(1.0);
       state.MPS->MPS_A = std::make_unique<class_mps_site>(M,L,0);
       state.MPS->MPS_B = std::make_unique<class_mps_site>(M,L,1);
       state.MPS->MPS_A->set_LC(L);
   }
   
   //void tools::infinite::mps::rebuild_environments(class_state_infinite & state){
   //
   //}
   
   
   
   
   class_state_infinite tools::infinite::mps::set_random_state(const class_state_infinite & state, [[maybe_unused]] std::string parity, [[maybe_unused]]  int seed_state){
       throw std::runtime_error("You need to implement set random state for infinite state");
       return state;
   }
