
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_environment.h:

Program Listing for File class_environment.h
============================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_environment.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/state/class_environment.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 7/21/17.
   //
   
   #ifndef DMRG_CLASS_ENVIRONMENT_H
   #define DMRG_CLASS_ENVIRONMENT_H
   
   #include <memory>
   
   #include "general/nmspc_tensor_extra.h"
   #include <spdlog/fmt/fmt.h>
   
   //using namespace Textra;
   //using namespace std;
   
   class class_vidal_site;
   
   class class_environment{
   public:
       using Scalar = std::complex<double>;
   private:
       std::optional<size_t> position;
       bool edge_has_been_set = false;
       void enlarge      (const Eigen::Tensor<Scalar,3>  & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
   public:
       std::string side;
       size_t sites = 0;                               
       Eigen::Tensor<Scalar,3> block;                 
       explicit class_environment(std::string side_):side(side_){};
       class_environment(
               std::string side_,
               int mpsDim,
               int mpoDim)
               :side(side_)
       {
           set_edge_dims(mpsDim,mpoDim);
       }
       bool isReal () const;
       void enlarge      (const class_vidal_site & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
       void set_edge_dims(const class_vidal_site & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
       void set_edge_dims(const Eigen::Tensor<Scalar,3> & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
       void set_edge_dims(const int mpsDim, const int mpoDim);
       void set_position(const size_t position_){position = position_;}
       size_t get_position() const {
           if(position) {return position.value();}
           else{throw std::runtime_error(fmt::format("Position hasn't been set on environment side {}", side));}
       }
   
   };
   
   
   
   class class_environment_var{
   public:
       using Scalar = std::complex<double>;
   private:
       std::optional<size_t> position;
       bool edge_has_been_set = false;
       void enlarge(const Eigen::Tensor<Scalar,3>  & MPS, const Eigen::Tensor<Scalar,4> &MPO);
   public:
       size_t sites = 0;                                      
       std::string side;
       Eigen::Tensor<Scalar,4> block;                         
       explicit class_environment_var(std::string side_):side(side_){};
       class_environment_var(
               std::string side_,
               int mpsDim,
               int mpoDim)
               :side(side_)
       {
           set_edge_dims(mpsDim,mpoDim);
       }
       bool isReal () const;
       void enlarge(const class_vidal_site & MPS, const Eigen::Tensor<Scalar,4> &MPO);
       void set_edge_dims(const class_vidal_site & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
       void set_edge_dims(const Eigen::Tensor<Scalar,3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO);
       void set_edge_dims(const int mpsDim, const int mpoDim);
       void set_position(const size_t position_){position = position_;}
       size_t get_position() const {
           if(position) {return position.value();}
           else{throw std::runtime_error("Position hasn't been set on environment var " + side);}
       }};
   
   #endif //DMRG_CLASS_ENVIRONMENT_H
