
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_model_factory.cpp:

Program Listing for File class_model_factory.cpp
================================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_model_factory.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/model/class_model_factory.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-07-04.
   //
   
   
   
   #include <math/nmspc_math.h>
   #include "class_model_factory.h"
   #include "class_tf_ising.h"
   #include "class_model_base.h"
   #include "class_hamiltonian_h5tables.h"
   #include "class_selfdual_tf_rf_ising.h"
   
   std::unique_ptr<class_model_base> class_model_factory::create_mpo(size_t position, std::string model_type_str){
   
       if (model_type_str == std::string("tf_ising")){
           return std::make_unique<class_tf_ising>(position,model_type_str);
       }
       else
       if (model_type_str == std::string("tf_nn_ising")){
           return std::make_unique<class_tf_ising>(position,model_type_str);
       }
       else
       if (model_type_str == std::string("selfdual_tf_rf_ising")){
           return std::make_unique<class_selfdual_tf_rf_ising>(position,model_type_str);
       }
       else{
           throw std::runtime_error("Wrong model: [ "  + model_type_str + " ]");
       }
   }
   
   
   
   std::unique_ptr<class_model_base> class_model_factory::clone(std::unique_ptr<class_model_base> other){
       return other->clone();
   }
   
