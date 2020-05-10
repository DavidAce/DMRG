//
// Created by david on 2019-06-24.
//

#include <tools/finite/opt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <state/class_state_finite.h>
#include <simulation/nmspc_settings.h>
#include <math/arpack_extra/matrix_product_hamiltonian.h>
#include <math/class_eigsolver.h>

Eigen::Tensor<class_state_finite::Scalar,4> tools::finite::opt::internal::ground_state_optimization(const class_state_finite & state, std::string ritzstring){
    tools::log->trace("Starting ground state optimization");
//    using Scalar = std::complex<double>;
    using namespace internal;
    using namespace settings::precision;
    using namespace eigutils::eigSetting;

    Ritz ritz = stringToRitz(ritzstring);
    auto dimsL = state.MPS_L.back().get_M().dimensions();
    auto dimsR = state.MPS_R.front().get_M().dimensions();
    std::array<long,4> shape_theta4  = {dimsL[0], dimsL[1], dimsR[0], dimsR[2]};
    std::array<long,4> shape_mpo4    =  state.MPO_L.back()->MPO().dimensions();

    tools::common::profile::t_eig->tic();
    int nev = 1;
    DenseHamiltonianProduct<Scalar>  matrix (
            state.ENV_L.back().block.data(),
            state.ENV_R.front().block.data(),
            state.MPO_L.back()->MPO().data(),
            state.MPO_R.front()->MPO().data(),
            shape_theta4,
            shape_mpo4,
            settings::threading::num_threads);

    class_eigsolver solver;
    solver.eigs_dense(matrix, nev, static_cast<int>(eig_max_ncv), NAN, Form::SYMMETRIC, ritz, Side::R, true, true);

    [[maybe_unused]] auto eigvals           = Eigen::TensorMap<const Eigen::Tensor<double,1>>  (solver.solution.get_eigvals<Form::SYMMETRIC>().data() ,solver.solution.meta.cols);
    [[maybe_unused]] auto eigvecs           = Eigen::TensorMap<const Eigen::Tensor<Scalar,1>>  (solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(),solver.solution.meta.rows);

    tools::common::profile::t_eig->toc();

    return eigvecs.reshape(shape_theta4);
}
