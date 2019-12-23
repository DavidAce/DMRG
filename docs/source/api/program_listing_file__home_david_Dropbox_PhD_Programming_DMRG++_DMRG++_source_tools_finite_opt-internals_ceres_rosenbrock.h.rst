
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_opt-internals_ceres_rosenbrock.h:

Program Listing for File ceres_rosenbrock.h
===========================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_opt-internals_ceres_rosenbrock.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/opt-internals/ceres_rosenbrock.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-12-10.
   //
   
   #pragma once
   
   #include <tools/finite/opt.h>
   
   namespace tools::finite::opt::internal {
           class ceres_rosenbrock_functor : public ceres::FirstOrderFunction {
           public:
               bool Evaluate(const double *parameters,
                                     double *cost,
                                     double *gradient) const final
               {
                   const double x = parameters[0];
                   const double y = parameters[1];
   
                   cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
                   if (gradient != nullptr) {
                       gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
                       gradient[1] = 200.0 * (y - x * x);
                   }
                   return true;
               }
   
               int NumParameters() const final { return 2; }
           };
       }
   
   
