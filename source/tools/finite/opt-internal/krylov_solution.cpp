//
// Created by david on 2021-04-17.
//
#include <tools/finite/opt-internal/opt-internal.h>
#include <math/eig/solver.h>
#include <math/eig/view.h>
#include <tools/finite/opt_mps.h>
#include <general/nmspc_tensor_extra.h>
#include <tools/finite/measure.h>
#include <tensors/class_tensors_finite.h>

void tools::finite::opt::internal::krylov_extract_solutions(const opt_mps &initial_mps, const class_tensors_finite &tensors, eig::solver &solver,
                                           std::vector<tools::finite::opt::opt_mps> &eigvecs_mps, const std::string &tag) {
    if(solver.result.meta.eigvals_found) {
        auto dims_mps = initial_mps.get_tensor().dimensions();
        auto eigvecs  = eig::view::get_eigvecs<Scalar>(solver.result, eig::Side::R, true);
        auto eigvals  = eig::view::get_eigvals<double>(solver.result, true);
        eigvecs_mps.reserve(eigvecs_mps.size() + static_cast<size_t>(eigvals.size()));
        for(long idx = 0; idx < eigvals.size(); idx++) {
            auto eigvec_i = Textra::asNormalized(Textra::TensorCast(eigvecs.col(idx), dims_mps));
            auto overlap  = std::abs(initial_mps.get_vector().dot(eigvecs.col(idx)));
            auto energy   = tools::finite::measure::energy(eigvec_i, tensors);
            auto eigval   = energy - initial_mps.get_energy_reduced();
            auto variance = tools::finite::measure::energy_variance(eigvec_i, tensors);
            eigvecs_mps.emplace_back(fmt::format("{:<8} eigenvector {}", tag,idx), eigvec_i, tensors.active_sites, eigval, initial_mps.get_energy_reduced(),
                                     variance, overlap, tensors.get_length());
            auto & mps = eigvecs_mps.back();
            mps.set_time(solver.result.meta.time_total);
            mps.set_counter(static_cast<size_t>(solver.result.meta.counter));
            mps.set_iter(static_cast<size_t>(solver.result.meta.iter));
            mps.is_basis_vector = true;
            mps.validate_candidate();
            mps.set_krylov_nev(solver.result.meta.nev_converged);
            mps.set_krylov_ncv(solver.result.meta.ncv);
            mps.set_krylov_tol(solver.result.meta.tol);
            mps.set_krylov_eigval(eigvals(idx));
            mps.set_krylov_ritz(solver.result.meta.ritz);
            mps.set_krylov_shift(solver.result.meta.sigma);
        }
    }
}