
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_opt-internals_ceres_rosenbrock.cpp:

Program Listing for File ceres_rosenbrock.cpp
=============================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_tools_finite_opt-internals_ceres_rosenbrock.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/tools/finite/opt-internals/ceres_rosenbrock.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2019-12-10.
   //
   #include "ceres_rosenbrock.h"
   #include <state/class_state_finite.h>
   
   
   Eigen::Tensor<class_state_finite::Scalar,3>
   tools::finite::opt::internal::ceres_rosenbrock_optimization(const class_state_finite & state) {
   
       Eigen::VectorXd parameters(2);
       parameters << -1.2, 1.0;
   
       ceres::GradientProblemSolver::Options options;
       options.minimizer_progress_to_stdout = true;
   
       ceres::GradientProblemSolver::Summary summary;
       auto *functor = new ceres_rosenbrock_functor();
       ceres::GradientProblem problem(functor);
       ceres::Solve(options, problem, parameters.data(), &summary);
   
       std::cout << summary.FullReport() << "\n";
       std::cout << "Initial x: " << -1.2 << " y: " << 1.0 << "\n";
       std::cout << "Final   x: " << parameters[0]
                 << " y: " << parameters[1] << "\n";
   
       return state.get_multitheta();
   }
