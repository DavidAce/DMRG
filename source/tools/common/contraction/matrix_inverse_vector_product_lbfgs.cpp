#include <math/tenx/fwd_decl.h>
// Eigen goes first
#include "debug/exceptions.h"
#include "tools/common/contraction.h"
#include "tools/finite/opt/lbfgs_simps_functor.h"
#include "tools/finite/opt/opt-internal.h"
#include <ceres/gradient_problem.h>
#include <ceres/gradient_problem_solver.h>
#include <complex>

/* clang-format off */
template<typename Scalar>
void tools::common::contraction::matrix_inverse_vector_product_lbfgs(
                                    Scalar * res_ptr,
                                    const Scalar * const mps_ptr, std::array<long,3> mps_dims,
                                    const Scalar * const mpo_ptr, std::array<long,4> mpo_dims,
                                    const Scalar * const envL_ptr, std::array<long,3> envL_dims,
                                    const Scalar * const envR_ptr, std::array<long,3> envR_dims){
    /* clang-format on */

    // Here we apply the inverse operation y = A⁻¹ * x
    // In this case we use lbfgs for unconstrained minimization of f = |Aφ - ψ|², where
    //      φ = res
    //      ψ = mps
    //      A = (H-E) = effective hamiltonian (from mpo and env)
    //
    // The gradient is ∇f = (H-E)²φ - (H-E)ψ, which is why we also need mpo2 and env2
    // After minimization we have φ ~ (H-E)⁻¹ψ = A⁻¹ * x

    auto res  = Eigen::TensorMap<Eigen::Tensor<Scalar, 3>>(res_ptr, mps_dims);
    auto mps  = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>(mps_ptr, mps_dims);
    auto mpo  = Eigen::TensorMap<const Eigen::Tensor<Scalar, 4>>(mpo_ptr, mpo_dims);
    auto envL = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>(envL_ptr, envL_dims);
    auto envR = Eigen::TensorMap<const Eigen::Tensor<Scalar, 3>>(envR_ptr, envR_dims);

    if(mps.dimension(1) != envL.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envL {}", mps.dimensions(), envL.dimensions());
    if(mps.dimension(2) != envR.dimension(0)) throw except::runtime_error("Dimension mismatch mps {} and envR {}", mps.dimensions(), envR.dimensions());
    if(mps.dimension(0) != mpo.dimension(2)) throw except::runtime_error("Dimension mismatch mps {} and mpo {}", mps.dimensions(), mpo.dimensions());
    if(envL.dimension(2) != mpo.dimension(0)) throw except::runtime_error("Dimension mismatch envL {} and mpo {}", envL.dimensions(), mpo.dimensions());
    if(envR.dimension(2) != mpo.dimension(1)) throw except::runtime_error("Dimension mismatch envR {} and mpo {}", envR.dimensions(), mpo.dimensions());

    auto  summary               = ceres::GradientProblemSolver::Summary();
    auto *functor               = new tools::finite::opt::internal::lbfgs_simps_functor<Scalar>(mps, envL, envR, mpo);
    auto  options               = tools::finite::opt::internal::lbfgs_default_options;
    options.max_lbfgs_rank      = 8;
    options.max_num_iterations  = 10000;
    options.function_tolerance  = 1e-6;
    options.gradient_tolerance  = 1e-2;
    options.parameter_tolerance = 1e-6;

    ceres::GradientProblem problem(functor);
    if constexpr(std::is_same_v<Scalar, double>) {
        ceres::Solve(options, problem, res.data(), &summary);
    } else {
        auto y_cplx_as_2x_real = Eigen::Map<Eigen::VectorXd>(reinterpret_cast<double *>(res.data()), 2 * res.size());
        ceres::Solve(options, problem, y_cplx_as_2x_real.data(), &summary);
    }
    if(not summary.iterations.empty()) {
        tools::log->info("LBFGS Preconditioner: f: {:8.5e} | size {} | time {:8.5e} | it {} | exit: {}. Message: {}", summary.final_cost,
                         summary.num_parameters, summary.total_time_in_seconds, summary.iterations.size(),
                         ceres::TerminationTypeToString(summary.termination_type), summary.message.c_str());
    }

    /* clang-format off */

}

using namespace tools::common::contraction;
template void tools::common::contraction::matrix_inverse_vector_product_lbfgs(
                                            real * res_ptr,
                                            const real * const mps_ptr, std::array<long,3> mps_dims,
                                            const real * const mpo_ptr, std::array<long,4> mpo_dims,
                                            const real * const envL_ptr, std::array<long,3> envL_dims,
                                            const real * const envR_ptr, std::array<long,3> envR_dims);
template void tools::common::contraction::matrix_inverse_vector_product_lbfgs(
                                            cplx * res_ptr,
                                            const cplx * const mps_ptr, std::array<long,3> mps_dims,
                                            const cplx * const mpo_ptr, std::array<long,4> mpo_dims,
                                            const cplx * const envL_ptr, std::array<long,3> envL_dims,
                                            const cplx * const envR_ptr, std::array<long,3> envR_dims);
