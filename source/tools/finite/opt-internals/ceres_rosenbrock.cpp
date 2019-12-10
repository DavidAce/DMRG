//
// Created by david on 2019-12-10.
//

#include "ceres_rosenbrock.h"
#include <state/class_state_finite.h>
Eigen::Tensor<class_state_finite::Scalar,3>
tools::finite::opt::internal::ceres_rosenbrock_optimization(const class_state_finite & state) {

    double parameters[2] = {-1.2, 1.0};

    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout = true;

    ceres::GradientProblemSolver::Summary summary;
    ceres::GradientProblem problem(new ceres_rosenbrock_functor());
    ceres::Solve(options, problem, parameters, &summary);

    std::cout << summary.FullReport() << "\n";
    std::cout << "Initial x: " << -1.2 << " y: " << 1.0 << "\n";
    std::cout << "Final   x: " << parameters[0]
              << " y: " << parameters[1] << "\n";

    return state.get_multitheta();
}
