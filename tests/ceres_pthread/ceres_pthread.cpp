#include <ceres/ceres.h>
#include "glog/logging.h"

class Rosenbrock : public ceres::FirstOrderFunction {
public:
    virtual bool Evaluate(const double* parameters,
                          double* cost,
                          double* gradient) const {
        const double x = parameters[0];
        const double y = parameters[1];

        cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
        if (gradient != NULL) {
            gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
            gradient[1] = 200.0 * (y - x * x);
        }
        return true;
    }

    virtual int NumParameters() const { return 2; }
};

int main(){
    google::InitGoogleLogging("ceres_pthread");
    double parameters[2] = {-1.2, 1.0};
    ceres::GradientProblemSolver::Options options;
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


    ceres::GradientProblemSolver::Summary summary;
    ceres::GradientProblem problem(new Rosenbrock());
    ceres::Solve(options, problem, parameters, &summary);
    std::cout << summary.FullReport() << "\n";
    std::cout << "Initial x: " << -1.2 << " y: " << 1.0 << "\n";
    std::cout << "Final   x: " << parameters[0]
              << " y: " << parameters[1] << "\n";
    return 0;

}