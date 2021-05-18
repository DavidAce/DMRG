
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <general/nmspc_tensor_extra.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_mpo.h>
#include <math/linalg/matrix.h>
#include <tensors/class_tensors_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/measure.h>
#include <tools/finite/opt-internal/opt-internal.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_mps.h>
tools::finite::opt::opt_mps tools::finite::opt::internal::krylov_energy_optimization(const class_tensors_finite &tensors, const opt_mps &initial_mps,
                                                                                     const class_algorithm_status &status, OptType optType, OptMode optMode,
                                                                                     OptSpace optSpace) {
    using namespace internal;
    using namespace settings::precision;
    auto t_eig = tools::common::profile::get_default_prof()["t_eig"]->tic_token();




    eig::Ritz   ritz      = eig::Ritz::SM;
    const auto &mpo       = tensors.get_multisite_mpo();
    const auto &env       = tensors.get_multisite_ene_blk();
    auto        shape_mps = tensors.active_problem_dims();
    auto        shape_mpo = mpo.dimensions();
    auto size = shape_mps[0]*shape_mps[1]*shape_mps[2];
    auto        nev       = std::clamp<eig::size_type>(size / 8, 4, 16);
    auto        ncv       = static_cast<eig::size_type>(settings::precision::eig_default_ncv);
    tools::log->info("Excited state energy optimization with ritz SM | dims {} = {}", shape_mps, size);

    tools::log->trace("Defining reduced Hamiltonian matrix-vector product");
    MatVecMPO<Scalar> matrix(env.L.data(), env.R.data(), mpo.data(), shape_mps, shape_mpo);
    tools::log->trace("Defining eigenvalue solver");
    eig::solver solver;
    solver.config.tol = 1e-12;

    // Since we use reduced energy mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
    // would otherwise cause trouble for Arpack. This equates to subtracting sigma * identity from the bottom corner of the mpo.
    // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
    tools::log->trace("Finding excited state");
    auto residual = initial_mps.get_tensor();
    solver.eigs(matrix, nev, ncv, ritz, eig::Form::SYMM, eig::Side::R, initial_mps.get_eigval(), eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF, residual.data());
    t_eig.toc();
    auto eigvecs = eig::view::get_eigvecs<Scalar>(solver.result);
    auto eigvals = eig::view::get_eigvals<double>(solver.result);
    std::vector<opt_mps> candidate_list;
    candidate_list.reserve(static_cast<size_t>(eigvals.size()));
    for(long idx = 0; idx < eigvals.size(); idx++) {
        auto eigvec_i = Textra::TensorCast(eigvecs.col(idx), shape_mps);
        auto overlap = std::abs(initial_mps.get_vector().dot(eigvecs.col(idx)));
        candidate_list.emplace_back(fmt::format("eigenvector {}", idx), eigvec_i, tensors.active_sites, eigvals(idx), initial_mps.get_energy_reduced(), std::nullopt,
                                    overlap, tensors.get_length());

        candidate_list.back().set_time(tools::common::profile::get_default_prof()["t_eig"]->get_last_interval());
        candidate_list.back().set_counter(static_cast<size_t>(solver.result.meta.counter));
        candidate_list.back().set_iter(static_cast<size_t>(solver.result.meta.iter));
        candidate_list.back().is_basis_vector = true;
        candidate_list.back().validate_candidate();

//        candidate_list.back().set_energy(tools::finite::measure::energy(candidate_list.back().get_tensor(),tensors));
        auto energy_check = tools::finite::measure::energy(candidate_list.back().get_tensor(),tensors);
        candidate_list.back().set_variance(tools::finite::measure::energy_variance(candidate_list.back().get_tensor(),tensors));
        tools::log->info("Candidate {:<2} overlap {:.16f} | E {:>20.16f} | E/L {:>20.16f} | eigval {:>20.16f} | Er {:>20.16f} | e_check  {:>20.16f} | variance: {:>20.16f}", idx,
                         candidate_list.back().get_overlap(),
                         candidate_list.back().get_energy(),
                         candidate_list.back().get_energy_per_site(),
                         candidate_list.back().get_eigval(),
                         candidate_list.back().get_energy_reduced(),
                         energy_check,
                         std::log10(candidate_list.back().get_variance()));
    }
    return internal::ceres_direct_optimization(tensors,initial_mps,status, optType,optMode,optSpace);

    // Construct H² as a matrix (expensive operation!)
    Eigen::MatrixXcd H2_subspace = internal::get_multisite_hamiltonian_squared_subspace_matrix<Scalar>(*tensors.model, *tensors.edges, candidate_list);
//    if(optType == OptType::REAL) H2_subspace = H2_subspace.real();
    fmt::print("H2_sub: \n{}\n",linalg::matrix::to_string(H2_subspace));
    auto current_level = tools::log->level();
    tools::log->set_level(spdlog::level::trace);
    auto optimized_tensor = tools::finite::opt::internal::ceres_optimize_subspace(tensors, initial_mps, candidate_list, H2_subspace, status, optType, optMode, optSpace);
    tools::finite::opt::internal::reports::print_bfgs_report();
    tools::finite::opt::internal::reports::print_time_report();
    tools::log->set_level(current_level);
    return optimized_tensor;
}
