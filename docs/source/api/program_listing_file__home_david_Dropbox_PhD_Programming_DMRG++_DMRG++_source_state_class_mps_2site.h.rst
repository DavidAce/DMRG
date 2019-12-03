
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_mps_2site.h:

Program Listing for File class_mps_2site.h
==========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_state_class_mps_2site.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/state/class_mps_2site.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #pragma once
   
   
   #include <memory>
   #include "general/nmspc_tensor_extra.h"
   
   
   
   
   
   class class_mps_site;
   
   class class_mps_2site
   {
   public:
       using Scalar = std::complex<double>;
   private:
   
   //    long spin_dimension = 2;                      /*!< Local (or physical) spin or qubit dimension, usually denoted \f$d\f$ elsewhere. */
       Eigen::Tensor<Scalar,3> tmp3;                 
       Eigen::Tensor<Scalar,1> tmp1;                 
       template< class T >
       std::unique_ptr<T> copy_unique(const std::unique_ptr<T>& source)
       {
           return source ? std::make_unique<T>(*source) : nullptr;
       }
   public:
   //    explicit class_mps_2site(std::string model_type_str);
       class_mps_2site() = default;
       explicit class_mps_2site(const class_mps_2site &other);
   
   
       bool swapped = false;                                  
       double truncation_error = 0;
   
       std::unique_ptr<class_mps_site> MPS_A;
       std::unique_ptr<class_mps_site> MPS_B;
       bool isReal()const;
       long chiA () const;
       long chiB () const;
       long chiC () const;
       long spindim () const;
       const Eigen::Tensor<Scalar,3> & A_bare()  const;
       const Eigen::Tensor<Scalar,3> & A()  const;
       const Eigen::Tensor<Scalar,3> & B()  const;
       Eigen::Tensor<Scalar,2> LC()     const;
       Eigen::Tensor<Scalar,3> GA()    const;
       Eigen::Tensor<Scalar,3> GB()    const;
       Eigen::Tensor<Scalar,2> LA()    const;
       Eigen::Tensor<Scalar,2> LB()    const;
       void set_mps(const Eigen::Tensor<Scalar,1> &LA,
                    const Eigen::Tensor<Scalar,3> &A,
                    const Eigen::Tensor<Scalar,1> &LC_,
                    const Eigen::Tensor<Scalar,3> &B,
                    const Eigen::Tensor<Scalar,1> &LB);
   
       Eigen::DSizes<long,4> dimensions() const;
   
       void initialize(int spin_dim);                                  
       void swap_AB();                                                 
       Eigen::Tensor<Scalar,4> get_theta (Scalar norm = 1.0)  const;   
   };
   
