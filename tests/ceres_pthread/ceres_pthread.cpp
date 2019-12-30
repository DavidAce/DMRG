//#define EIGEN_MALLOC_ALREADY_ALIGNED 0



#include "ceres_pthread.h"
void opt::SolveRosenbrock(){
    google::InitGoogleLogging("ceres_pthread");
    ceres::GradientProblemSolver::Summary summary;

    ceres::GradientProblemSolver::Options options = ceres_default_options;

    options.line_search_type = ceres::LineSearchType::WOLFE;
    options.line_search_interpolation_type = ceres::LineSearchInterpolationType::CUBIC;
    options.line_search_direction_type = ceres::LineSearchDirectionType::LBFGS;
    options.nonlinear_conjugate_gradient_type = ceres::NonlinearConjugateGradientType::POLAK_RIBIERE;
    options.max_num_iterations = 500;
    options.max_lbfgs_rank     = 250;
    options.use_approximate_eigenvalue_bfgs_scaling = true;  // True makes a huge difference, takes longer steps at each iteration!!
    options.max_line_search_step_expansion = 1e4;// 100.0;
    options.min_line_search_step_size = 1e-256;//  std::numeric_limits<double>::epsilon();
    options.max_line_search_step_contraction = 1e-3;
    options.min_line_search_step_contraction = 0.6;
    options.max_num_line_search_step_size_iterations  = 50;//20;
    options.max_num_line_search_direction_restarts    = 20;//2;
    options.line_search_sufficient_function_decrease  = 1e-2;
    options.line_search_sufficient_curvature_decrease = 0.4; //0.5;
    options.max_solver_time_in_seconds = 60*10;//60*2;
    options.function_tolerance = 1e-4;
    options.gradient_tolerance = 1e-4;
    options.parameter_tolerance = 1e-256;//std::numeric_limits<double>::epsilon();//1e-12;
    options.minimizer_progress_to_stdout = true;
    options.logging_type = ceres::LoggingType::PER_MINIMIZER_ITERATION;
    Eigen::VectorXd parameters = Eigen::VectorXd::Random(100);
    Eigen::MatrixXd H(100,100);
    H.setRandom();
    H = H.selfadjointView<Eigen::Upper>();
    Eigen::MatrixXd H2 = H*H;


//    auto * functor = new opt::Rosenbrock<double>(H);
    auto * functor = new opt::Rosenbrock<double>(H,H2);
    ceres::GradientProblem problem(functor);
    ceres::Solve(options, problem, parameters.data(), &summary);
    std::cout << summary.FullReport() << "\n";
    int iter  = (int)summary.iterations.size();

    std::cout << "Initial x: " << -1.2 << " y: " << 1.0 << "\n";
    std::cout << "Final   x: " << parameters[0]
              << " y: " << parameters[1] << "\n";
    std::cout << "iterations: " << iter << std::endl;
    std::cout << "Finished LBFGS after " << summary.total_time_in_seconds << " seconds " << " and " << summary.iterations.size() << " iters " << std::endl;
    std::cout << "Exit status " << ceres::TerminationTypeToString(summary.termination_type) << std::endl;
    std::cout << "Message: " << summary.message.c_str() << std::endl;
}


int main(){

    opt::SolveRosenbrock();
    return 0;

}