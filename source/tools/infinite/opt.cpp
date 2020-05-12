//
// Created by david on 2019-07-06.
//

#include <math/arpack_extra/matrix_product_hamiltonian.h>
#include <math/class_eigsolver.h>
#include <math/class_svd_wrapper.h>
#include <model/class_mpo_base.h>
#include <state/class_environment.h>
#include <state/class_mps_2site.h>
#include <state/class_mps_site.h>
#include <state/class_state_infinite.h>
#include <tools/common/prof.h>
#include <tools/infinite/opt.h>

using Scalar = tools::infinite::opt::Scalar;

Eigen::Tensor<Scalar, 4> tools::infinite::opt::find_ground_state(const class_state_infinite &state, StateRitz ritz){
    return tools::infinite::opt::find_ground_state(state,enum2str(ritz));
}

Eigen::Tensor<Scalar, 4> tools::infinite::opt::find_ground_state(const class_state_infinite &state, std::string_view ritzstring) {
    auto                theta        = state.get_theta();
    std::array<long, 4> shape_theta4 = theta.dimensions();
    std::array<long, 4> shape_mpo4   = state.HA->MPO().dimensions();

    tools::common::profile::t_eig->tic();
    int nev = 1;
    using namespace settings::precision;
    using namespace eigutils::eigSetting;
    Ritz ritz = stringToRitz(ritzstring);

    DenseHamiltonianProduct<Scalar> matrix(state.Lblock->block.data(), state.Rblock->block.data(), state.HA->MPO().data(), state.HB->MPO().data(), shape_theta4,
                                           shape_mpo4, settings::threading::num_threads);
    class_eigsolver                 solver;
    solver.eigs_dense(matrix, nev, static_cast<int>(eig_max_ncv), NAN, Form::SYMMETRIC, ritz, Side::R, true, true);
    auto eigvec =
        Eigen::TensorMap<const Eigen::Tensor<Scalar, 1>>(solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(), solver.solution.meta.rows);

    tools::common::profile::t_eig->toc();

    return eigvec.reshape(theta.dimensions());
}

//============================================================================//
// Do unitary evolution on an MPS
//============================================================================//

Eigen::Tensor<Scalar, 4> tools::infinite::opt::time_evolve_theta(const class_state_infinite &state, const Eigen::Tensor<Scalar, 4> &U)
/*!
@verbatim
  1--[ Θ ]--3
     |   |
     0   2
                   1--[ Θ ]--3
     0   1   --->     |   |
     |   |            0   2
     [ U ]
     |   |
     2   3
@endverbatim
*/
{
    return U.contract(state.get_theta(), Textra::idx({0, 1}, {0, 2})).shuffle(Textra::array4{0, 2, 1, 3});
}

void tools::infinite::opt::truncate_theta(Eigen::Tensor<Scalar, 4> &theta, class_state_infinite &state) {
    class_SVD SVD;
    SVD.setThreshold(settings::precision::svd_threshold);
    auto [U, S, V]              = SVD.schmidt(theta, state.get_chi_lim());
    state.MPS->truncation_error = SVD.get_truncation_error();
    state.MPS->MPS_A->set_LC(S);
    state.MPS->MPS_A->set_M(U);
    state.MPS->MPS_B->set_M(V);
}
