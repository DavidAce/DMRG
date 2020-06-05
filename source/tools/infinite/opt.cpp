//
// Created by david on 2019-07-06.
//

#include <eig/arpack_extra/matrix_product_hamiltonian.h>
#include <eig/eig.h>
#include <tensors/class_tensors_infinite.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/prof.h>
#include <tools/infinite/opt.h>

using Scalar = tools::infinite::opt::Scalar;

Eigen::Tensor<Scalar, 3> tools::infinite::opt::find_ground_state(const class_tensors_infinite &state, StateRitz ritz) {
    return tools::infinite::opt::find_ground_state(state, enum2str(ritz));
}

Eigen::Tensor<Scalar, 3> tools::infinite::opt::find_ground_state(const class_tensors_infinite &tensors, std::string_view ritzstring) {
    tools::log->trace("Starting ground state optimization");

    eig::Ritz ritz = eig::stringToRitz(ritzstring);

    auto        shape_mps = tensors.state->dimensions();
    auto        shape_mpo = tensors.model->dimensions();
    const auto &mpo       = tensors.model->get_2site_tensor();
    const auto &env       = tensors.edges->get_ene_blk();
    auto        nev       = static_cast<eig::size_type>(1);
    auto        ncv       = static_cast<eig::size_type>(settings::precision::eig_max_ncv);

    tools::common::profile::t_eig->tic();
    MatrixProductHamiltonian<Scalar> matrix(env.L.data(), env.R.data(), mpo.data(), shape_mps, shape_mpo, settings::threading::num_threads);
    eig::solver                      solver;
    solver.eigs(matrix, nev, ncv, ritz, eig::Form::SYMM, eig::Side::R, std::nullopt, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::ON);
    tools::common::profile::t_eig->toc();
    return eig::view::get_eigvec<Scalar>(solver.result, shape_mps);
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
