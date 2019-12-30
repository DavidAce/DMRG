
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_model_factory.h:

Program Listing for File class_model_factory.h
==============================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_model_factory.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/model/class_model_factory.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-07-04.
   //
   
   #pragma once
   
   #include <iostream>
   #include <memory>
   #include "class_model_base.h"
   #include "class_hamiltonian_h5tables.h"
   
   class class_model_factory{
   public:
       static std::unique_ptr<class_model_base>         create_mpo(size_t position,std::string model_type_str);
   
       template <typename... T>
       static std::unique_ptr<class_model_base>         create_mpo(size_t position,std::string model_type_str, T... args){
           auto mpo = create_mpo(position,model_type_str);
           mpo->set_hamiltonian(args...);
           return mpo;
       }
   //    static std::unique_ptr<class_hamiltonian_h5log_base> create_table(std::string model_type_str);
       static std::unique_ptr<class_model_base>         clone(std::unique_ptr<class_model_base> other);
   };
   
   
   
