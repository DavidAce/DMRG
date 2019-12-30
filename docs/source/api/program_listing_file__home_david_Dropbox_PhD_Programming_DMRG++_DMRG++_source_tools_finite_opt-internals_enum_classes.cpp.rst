
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_opt-internals_enum_classes.cpp:

Program Listing for File enum_classes.cpp
=========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_opt-internals_enum_classes.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/opt-internals/enum_classes.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-10-03.
   //
   
   #include "enum_classes.h"
   
   void tools::finite::opt::OptType::print(std::ostream& str)const {
       switch (this->option){
       case opt::TYPE::REAL : str << "REAL"; break;
       case opt::TYPE::CPLX : str << "CPLX"; break;
       default: throw std::runtime_error("Enum TYPE not set");
       }
   }
   
   void tools::finite::opt::OptMode::print(std::ostream& str)const {
       switch (this->option){
       case opt::MODE::OVERLAP  : str << "OVERLAP"; break;
       case opt::MODE::VARIANCE : str << "VARIANCE"; break;
       default: throw std::runtime_error("Enum MODE not set");
       }
   }
   
   
   void tools::finite::opt::OptSpace::print(std::ostream& str)const {
       switch (this->option){
       case opt::SPACE::SUBSPACE : str << "SUBSPACE"; break;
       case opt::SPACE::DIRECT   : str << "DIRECT"  ; break;
       default: throw std::runtime_error("Enum SPACE not set");
       }
   }
