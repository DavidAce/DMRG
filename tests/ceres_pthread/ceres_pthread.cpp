#include <ceres/ceres.h>
#include "glog/logging.h"
#include <general/nmspc_omp.h>

class RosenbrockBase : public ceres::FirstOrderFunction {
protected:
    Eigen::MatrixXd H;
    Eigen::MatrixXd H2;
    OMP omp;
public:
    explicit RosenbrockBase(Eigen::MatrixXd & H_):H(H_), H2(H_*H_){

    }
    int NumParameters() const final { return H.rows(); }

};


template<typename T>
class Rosenbrock : public RosenbrockBase {
private:
    using MatrixType = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> ;
    using VectorType = Eigen::Matrix<T,Eigen::Dynamic,1> ;
public:
    explicit Rosenbrock(MatrixType &H_): RosenbrockBase(H_){}
    bool Evaluate(const double* v_ptr,
                  double* fx,
                  double* grad_ptr) const final
    {
        Eigen::Map<const VectorType>  v    (v_ptr, NumParameters());
        double vv   = v.squaredNorm();
        double norm = std::sqrt(vv);
        VectorType Hv = H*v;
        VectorType H2v = H2*v;
        double ene             = v.adjoint()*Hv;
        double ene2            = v.adjoint()*H2v;

        double var             = std::abs(ene2 - ene);
        double norm_offset     = std::abs(1-norm);
        double log10var        = std::log10(var);
        if(fx != nullptr){
            fx[0] = log10var + norm_offset;
        }

        if (grad_ptr != nullptr){
            Eigen::Map<VectorType>  grad (grad_ptr, NumParameters());
            auto vv_1  = std::pow(vv,-1);
            auto var_1 = 1.0/var/std::log(10);
            grad = var_1 * vv_1 * 2.0*(H2v - 2.0*ene*Hv - (ene2 - 2.0*ene*ene)*v);
            grad += 2.0*norm_offset * v;
        }
        return true;
    }

};

int main(){
    google::InitGoogleLogging("ceres_pthread");
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
    Eigen::VectorXd parameters = Eigen::VectorXd::Random(100);
    Eigen::MatrixXd H(100,100);
    H.setRandom();
    H = H.selfadjointView<Eigen::Upper>();
    ceres::GradientProblemSolver::Summary summary;
    ceres::GradientProblem problem(new Rosenbrock<double>(H));
    ceres::Solve(options, problem, parameters.data(), &summary);
    std::cout << summary.FullReport() << "\n";

    return 0;

}