
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_model_base.cpp:

Program Listing for File class_model_base.cpp
=============================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_model_class_model_base.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/model/class_model_base.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-07-04.
   //
   
   #include "class_model_base.h"
   #include <general/nmspc_tensor_extra.h>
   #include <general/nmspc_quantum_mechanics.h>
   #include <general/nmspc_random_numbers.h>
   #include <iomanip>
   
   using namespace qm;
   using Scalar = std::complex<double>;
   
   const Eigen::Tensor<Scalar,4> & class_model_base::MPO() const{
       if (all_mpo_parameters_have_been_set){
           return mpo_internal;
       }else{
           throw std::runtime_error("All MPO parameters haven't been set yet.");
       }
   }
   
   bool class_model_base::isReal()const{
       return Textra::isReal(MPO(),"MPO");
   }
   
   class_model_base::class_model_base(size_t position_, std::string logName){
       position = position_;
       log = Logger::setLogger(logName);
   }
   
   void class_model_base::set_position(size_t position_){
       position = position_;
   }
   
   
   
   size_t class_model_base::get_position() const{
       if (position){return position.value();}
       else{throw std::runtime_error("Position of MPO has not been set");}
   }
   
   Eigen::MatrixXcd class_model_base::MPO_matrix_view(){
       auto rows = MPO().dimension(0)*MPO().dimension(2);
       auto cols = MPO().dimension(1)*MPO().dimension(3);
       Eigen::Tensor<Scalar,4> MPO_temp = MPO().shuffle(Textra::array4{0,2,1,3});
       return Textra::Tensor_to_Matrix<Scalar>(MPO_temp, rows ,cols);
   }
