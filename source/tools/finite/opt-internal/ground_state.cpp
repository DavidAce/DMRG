//
// Created by david on 2019-06-24.
//

#include <config/nmspc_settings.h>
#include <math/eig.h>
#include <math/eig/arpack_solver/matrix_product_hamiltonian.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt.h>

Eigen::Tensor<class_state_finite::Scalar,3> tools::finite::opt::internal::ground_state_optimization(const class_tensors_finite & tensors, StateRitz ritz){
    return ground_state_optimization(tensors,enum2str(ritz));
}

Eigen::Tensor<class_state_finite::Scalar,3> tools::finite::opt::internal::ground_state_optimization(const class_tensors_finite & tensors, std::string_view ritzstring){
    tools::log->trace("Ground state optimization with ritz {} ...", ritzstring);
    tools::log->info("clearing caches");
    tensors.clear_cache();
//    using Scalar = std::complex<double>;
    using namespace internal;
    using namespace settings::precision;

    eig::Ritz ritz = eig::stringToRitz(ritzstring);
    auto shape_mps = tensors.state->active_dimensions();
    auto shape_mpo = tensors.model->active_dimensions();
    const auto & mpo = tensors.model->get_multisite_tensor();
    const auto & env = tensors.edges->get_multisite_ene_blk();
//    int nev = std::min(4l,(long)(tensors.state->active_problem_size()/2));
    auto        nev       = static_cast<eig::size_type>(1);
    auto        ncv       = static_cast<eig::size_type>(settings::precision::eig_max_ncv);
    tools::common::profile::t_eig->tic();
    MatrixProductHamiltonian<Scalar>  matrix (
            env.L.data(),
            env.R.data(),
            mpo.data(),
            shape_mps,
            shape_mpo);

    eig::solver solver(0);
    // Since we use reduced energy mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
    // would otherwise cause trouble for Arpack. This equates to subtracting sigma * identity from the bottom corner of the mpo.
    // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
    solver.eigs(matrix, nev, ncv, ritz, eig::Form::SYMM, eig::Side::R, 1.0,  eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF);
    tools::common::profile::t_eig->toc();
    return eig::view::get_eigvec<Scalar>(solver.result,shape_mps, 0);
}
