#include <ceres/ceres.h>

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

int main (){
    double parameters[2] = {-1.2, 1.0};

    ceres::GradientProblem problem(new Rosenbrock());

    ceres::GradientProblemSolver::Options options;
    options.minimizer_progress_to_stdout = true;
    ceres::GradientProblemSolver::Summary summary;
    ceres::Problem::Options::cost_function_ownership = ceres::Ownership::DO_NOT_TAKE_OWNERSHIP;
    ceres::Solve(options, problem, parameters, &summary);

    std::cout << summary.FullReport() << "\n";
}