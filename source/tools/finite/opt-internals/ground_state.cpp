//
// Created by david on 2019-06-24.
//

#include <tools/finite/opt.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tensors/state/class_state_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/class_tensors_finite.h>
#include <config/nmspc_settings.h>
#include <math/arpack_extra/matrix_product_hamiltonian.h>
#include <math/class_eigsolver.h>

Eigen::Tensor<class_state_finite::Scalar,3> tools::finite::opt::internal::ground_state_optimization(const class_tensors_finite & tensors, StateRitz ritz){
    return ground_state_optimization(tensors,enum2str(ritz));
}

Eigen::Tensor<class_state_finite::Scalar,3> tools::finite::opt::internal::ground_state_optimization(const class_tensors_finite & tensors, std::string_view ritzstring){
    tools::log->trace("Starting ground state optimization");
//    using Scalar = std::complex<double>;
    using namespace internal;
    using namespace settings::precision;
    using namespace eigutils::eigSetting;

    Ritz ritz = stringToRitz(ritzstring);
    auto shape_mps = tensors.state->active_dimensions();
    auto shape_mpo = tensors.model->active_dimensions();
    const auto & mpo = tensors.model->get_multisite_tensor();
    const auto & env = tensors.edges->get_multisite_ene_blk();

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
    solver.eigs_dense(matrix, nev, static_cast<int>(eig_max_ncv), NAN, Form::SYMMETRIC, ritz, Side::R, true, true);

    [[maybe_unused]] auto eigvals = Eigen::TensorMap<const Eigen::Tensor<double,1>>  (solver.solution.get_eigvals<Form::SYMMETRIC>().data() ,solver.solution.meta.cols);
    [[maybe_unused]] auto eigvecs = Eigen::TensorMap<const Eigen::Tensor<Scalar,1>>  (solver.solution.get_eigvecs<Type::CPLX, Form::SYMMETRIC>().data(),solver.solution.meta.rows);

    tools::common::profile::t_eig->toc();
    return eigvecs.reshape(shape_mps);
}
