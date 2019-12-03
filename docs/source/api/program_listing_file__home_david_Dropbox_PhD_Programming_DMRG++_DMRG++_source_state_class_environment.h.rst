
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_environment.h:

Program Listing for File class_environment.h
============================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_environment.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/state/class_environment.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 7/21/17.
   //
   
   #pragma once
   
   #include <memory>
   
   #include "general/nmspc_tensor_extra.h"
   #include <spdlog/fmt/fmt.h>
   
   //using namespace Textra;
   //using namespace std;
   
   class class_mps_site;
   class class_model_base;
   
   class class_environment_base{
   public:
       using Scalar = std::complex<double>;
   protected:
       std::optional<size_t> position;
       bool edge_has_been_set = false;
       virtual void enlarge      (const Eigen::Tensor<Scalar,3>  & MPS, const Eigen::Tensor<Scalar, 4> &MPO) = 0;
   public:
       std::string side;
       size_t sites = 0;                               
       explicit class_environment_base(std::string side_, size_t position_);
       explicit class_environment_base(std::string side_, const class_mps_site & MPS, const class_model_base &MPO);
   
       virtual bool isReal () const = 0;
       virtual void set_edge_dims(const class_mps_site & MPS, const class_model_base &MPO)                   = 0;
       void set_position(const size_t position_){position = position_;}
       size_t get_position() const {
           if(position) {return position.value();}
           else{throw std::runtime_error(fmt::format("Position hasn't been set on environment side {}", side));}
       }
   };
   
   class class_environment final : public class_environment_base{
   private:
       void enlarge      (const Eigen::Tensor<Scalar,3>  & MPS, const Eigen::Tensor<Scalar, 4> &MPO)   final;
   public:
       Eigen::Tensor<Scalar,3> block;                 
       using class_environment_base::class_environment_base;
       explicit class_environment(std::string side_, const class_mps_site & MPS, const class_model_base &MPO);
       [[nodiscard]] class_environment enlarge(const class_mps_site & MPS, const class_model_base &MPO);
   
       bool isReal () const                                                                            final;
       void set_edge_dims(const class_mps_site & MPS, const class_model_base &MPO)                     final;
   
   };
   
   
   class class_environment_var final : public class_environment_base{
   private:
       void enlarge      (const Eigen::Tensor<Scalar,3>  & MPS, const Eigen::Tensor<Scalar, 4> &MPO)   final;
   public:
       Eigen::Tensor<Scalar,4> block;                 
       using class_environment_base::class_environment_base;
       explicit class_environment_var(std::string side_, const class_mps_site & MPS, const class_model_base &MPO);
       [[nodiscard]] class_environment_var enlarge(const class_mps_site & MPS, const class_model_base &MPO);
       bool isReal () const                                                                            final;
       void set_edge_dims(const class_mps_site & MPS, const class_model_base &MPO)                     final;
   };
   
   //
   //
   //class class_environment_var{
   //public:
   //    using Scalar = std::complex<double>;
   //private:
   //    std::optional<size_t> position;
   //    bool edge_has_been_set = false;
   //    void enlarge(const Eigen::Tensor<Scalar,3>  & MPS, const Eigen::Tensor<Scalar,4> &MPO);
   //    void set_edge_dims(const class_mps_site & MPS, const class_model_base &MPO);
   //    void set_edge_dims(const Eigen::Tensor<Scalar,3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO);
   //    void set_edge_dims(const int mpsDim, const int mpoDim);
   //public:
   //    size_t sites = 0;                                      /*!< Number of particles that have been contracted into this left environment. */
   //    std::string side;
   //    Eigen::Tensor<Scalar,4> block;                         /*!< The environment block. */
   //    explicit class_environment_var(std::string side_, size_t position_):position(position_),side(side_){};
   //
   //    bool isReal () const;
   //    class_environment_var enlarge(const class_mps_site & MPS, const class_model_base &MPO);
   //
   //
   //    void set_position(const size_t position_){position = position_;}
   //    size_t get_position() const {
   //        if(position) {return position.value();}
   //        else{throw std::runtime_error("Position hasn't been set on environment var " + side);}
   //    }};
   //
