//
// Created by david on 2019-07-06.
//

#include <math/arpack_extra/matrix_product_hamiltonian.h>
#include <math/class_eigsolver.h>
#include <math/class_svd_wrapper.h>
#include <tensors/class_tensors_infinite.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_environment.h>
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/prof.h>
#include <tools/common/views.h>
#include <tools/infinite/opt.h>

using Scalar = tools::infinite::opt::Scalar;

Eigen::Tensor<Scalar, 3> tools::infinite::opt::find_ground_state(const class_tensors_infinite &state, StateRitz ritz){
    return tools::infinite::opt::find_ground_state(state,enum2str(ritz));
}

Eigen::Tensor<Scalar, 3> tools::infinite::opt::find_ground_state(const class_tensors_infinite &tensors, std::string_view ritzstring) {
//    const auto                mps = state.get_mps_site();
//    std::array<long, 3> shape_mps   = tensors.state->dimensions();
//    std::array<long, 4> shape_mpo   = tensors.model->dimensions();
//    tools::common::profile::t_eig->tic();
//    int nev = 1;
//    using namespace settings::precision;
//    using namespace eigutils::eigSetting;
//    Ritz ritz = stringToRitz(ritzstring);
//
//    std::cerr << "split up superblock!" << std::endl;
////    DenseHamiltonianProduct<Scalar> matrix(state.Lblock->block.data(), state.Rblock->block.data(), state.HA->MPO().data(), state.HB->MPO().data(), shape_theta4,
////                                           shape_mpo4, settings::threading::num_threads);
//    class_eigsolver                 solver;
////    solver.eigs_dense(matrix, nev, static_cast<int>(eig_max_ncv), NAN, Form::SYMMETRIC, ritz, Side::R, true, true);
//    auto eigvec =
//        Eigen::TensorMap<const Eigen::Tensor<Scalar, 1>>(solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(), solver.solution.meta.rows);
//
//    tools::common::profile::t_eig->toc();
//    return eigvec.reshape(theta.dimensions());
//






    tools::log->trace("Starting ground state optimization");

//
    eigutils::eigSetting::Ritz ritz = eigutils::eigSetting::stringToRitz(ritzstring);
    auto shape_mps = tensors.state->dimensions();
    auto shape_mpo = tensors.model->dimensions();
    const auto & mpo = tensors.model->get_2site_tensor();
    const auto & env = tensors.edges->get_ene_blk();

    tools::common::profile::t_eig->tic();
    int nev = 1;
    DenseHamiltonianProduct<Scalar>  matrix (
        env.L.data(),
        env.R.data(),
        mpo.data(),
        shape_mps,
        shape_mpo,
        settings::threading::num_threads);

    class_eigsolver solver;
    solver.eigs_dense(matrix, nev, static_cast<int>(settings::precision::eig_max_ncv), NAN, eigutils::eigSetting::Form::SYMMETRIC, ritz, eigutils::eigSetting::Side::R, true, true);

    [[maybe_unused]] auto eigvals = Eigen::TensorMap<const Eigen::Tensor<double,1>>  (solver.solution.get_eigvals<eigutils::eigSetting::Form::SYMMETRIC>().data() ,solver.solution.meta.cols);
    [[maybe_unused]] auto eigvecs = Eigen::TensorMap<const Eigen::Tensor<Scalar,1>>  (solver.solution.get_eigvecs<eigutils::eigSetting::Type::CPLX, eigutils::eigSetting::Form::SYMMETRIC>().data(),solver.solution.meta.rows);

    tools::common::profile::t_eig->toc();
    return eigvecs.reshape(shape_mps);
}

//============================================================================//
// Do unitary evolution on an MPS
//============================================================================//

Eigen::Tensor<Scalar, 3> tools::infinite::opt::time_evolve_state(const class_state_infinite &state, const Eigen::Tensor<Scalar, 2> &U)
/*!
@verbatim
  1--[ mps ]--2
        |
        0
                         1--[ mps ]--2
        0         --->         |
        |                      0
      [ U ]
        |
        1
@endverbatim
*/
{
    return U.contract(state.get_2site_tensor(), Textra::idx({0}, {0}));
}

