
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_opt-internals_ceres_subspace_functor.h:

Program Listing for File ceres_subspace_functor.h
=================================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_opt-internals_ceres_subspace_functor.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/opt-internals/ceres_subspace_functor.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-07-15.
   //
   
   #pragma once
   
   
   #include <tools/finite/opt.h>
   
   namespace tools::finite::opt{
       namespace internal{
           template<typename Scalar>
           class ceres_subspace_functor : public ceres_base_functor{
           private:
               using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
               using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
   //            const Eigen::MatrixXcd &eigvecs;
               const MatrixType       &H2;
               const Eigen::VectorXd  &eigvals;
           public:
   //            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
               explicit ceres_subspace_functor(const class_state_finite & state,
                                               const class_simulation_status & sim_status,
                                               const MatrixType & H2_subspace_,
                                               const Eigen::VectorXd  & eigvals_);
               bool Evaluate(const double* v_double_double,
                             double* fx,
                             double* grad_double_double) const override;
           };
       }
   }
