//
// Created by david on 2019-06-24.
//

#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <math/eig.h>
#include <math/eig/arpack_solver/matrix_product_hamiltonian.h>
#include <tensors/class_tensors_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt-internal/opt-internal.h>

Eigen::Tensor<class_tensors_finite::Scalar, 3> tools::finite::opt::internal::ground_state_optimization(const class_tensors_finite &tensors, StateRitz ritz) {
    return ground_state_optimization(tensors, enum2str(ritz));
}

Eigen::Tensor<class_tensors_finite::Scalar, 3> tools::finite::opt::internal::ground_state_optimization(const class_tensors_finite &tensors,
                                                                                                       std::string_view            ritzstring) {
    tools::log->debug("Ground state optimization with ritz {} ...", ritzstring);
    using namespace internal;
    using namespace settings::precision;
    auto t_eig = tools::common::profile::get_default_prof()["t_eig"]->tic_token();

    eig::Ritz   ritz      = eig::stringToRitz(ritzstring);
    const auto &mpo       = tensors.get_multisite_mpo();
    const auto &env       = tensors.get_multisite_ene_blk();
    auto        shape_mps = tensors.active_problem_dims();
    auto        shape_mpo = mpo.dimensions();
    //    int nev = std::min(4l,(long)(tensors.state->active_problem_size()/2));
    auto nev = static_cast<eig::size_type>(1);
    auto ncv = static_cast<eig::size_type>(settings::precision::eig_max_ncv);
    tools::log->trace("Defining Hamiltonian matrix-vector product");
    MatrixProductHamiltonian<Scalar> matrix(env.L.data(), env.R.data(), mpo.data(), shape_mps, shape_mpo);
    tools::log->trace("Defining eigenvalue solver");
    eig::solver solver;

    if(tensors.measurements.energy_variance)
        solver.config.eigThreshold = std::clamp( 1e-2 * tensors.measurements.energy_variance.value(), 1e-16 , settings::precision::eig_threshold);

    // Since we use reduced energy mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
    // would otherwise cause trouble for Arpack. This equates to subtracting sigma * identity from the bottom corner of the mpo.
    // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
    tools::log->trace("Finding ground state");
    try{
        solver.eigs(matrix, nev, ncv, ritz, eig::Form::SYMM, eig::Side::R, 1.0, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF);
        return Textra::asNormalized(eig::view::get_eigvec<Scalar>(solver.result, shape_mps, 0));
    }catch(const std::exception & ex){
        tools::log->warn("Arpack failed. shape_mps {} | shape_mpo {} | sites {}", shape_mps, shape_mpo, tensors.active_sites);
        return Textra::asNormalized(tensors.get_multisite_mps());
    }

}
