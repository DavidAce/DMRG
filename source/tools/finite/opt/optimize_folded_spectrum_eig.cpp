#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/eig.h"
#include "math/num.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/log.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"

namespace tools::finite::opt::internal {

    template<typename Scalar>
    void optimize_folded_spectrum_eig_executor(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, OptMeta &meta) {
        // Solve the folded spectrum problem H²v=E²x (H² is positive definite)
        if(meta.optRitz == OptRitz::NONE) return;
        eig::solver solver;
        auto        matrix = tensors.get_effective_hamiltonian_squared<Scalar>();
        auto        nev    = std::min<int>(static_cast<int>(matrix.dimension(0)), meta.eigs_nev.value_or(1));
        auto        il     = 1;
        auto        iu     = nev;
        switch(meta.optRitz) {
            case OptRitz::LR: [[fallthrough]];
            case OptRitz::LM: {
                il = static_cast<int>(matrix.dimension(0) - (nev - 1));
                iu = static_cast<int>(matrix.dimension(0));
                break;
            }
            default: break;
        }
        solver.eig(matrix.data(), matrix.dimension(0), 'I', il, iu, 0.0, 1.0);
        extract_results(tensors, initial_mps, meta, solver, results, true);

        // if(solver.result.meta.eigvals_found and solver.result.meta.eigvecsR_found) {
        //     // tools::log->info("optimize_variance_eig_executor: vl {:.3e} | vu {:.3e}", vl, vu);
        //     tools::log->info("Found {} eigvals ({} converged)", solver.result.meta.nev, solver.result.meta.nev_converged);
        //     auto eigvals = eig::view::get_eigvals<real>(solver.result, false);
        //     auto indices = num::range<long>(0l, eigvals.size());
        //     auto eigComp = EigIdxComparator(OptRitz::LM, 0, eigvals.data(), eigvals.size());
        //     std::sort(indices.begin(), indices.end(), eigComp); // Should sort them according to distance from eigval
        //     indices.resize(std::min(eigvals.size(), 10l));      // We only need the first few indices, say 4
        //     for(auto idx : indices) { tools::log->info(" -- idx {}: {:.16f}", idx, eigvals(idx)); }
        // }

        // if(meta.chosen_sites.size() <= 4 and matrix.dimension(0) <= 8192) {
        //     // Find all eigenvalues within a thin band
        // auto eigval = initial_mps.get_energy(); // The current energy
        // auto eigvar = initial_mps.get_variance();
        // auto eshift = initial_mps.get_eshift();                    // The energy shift is our target energy for excited states
        // auto vl     = eshift - std::abs(eigval) - 2 * std::sqrt(eigvar); // Find energies at most two sigma away from the band
        // auto vu     = eshift + std::abs(eigval) + 2 * std::sqrt(eigvar); // Find energies at most two sigma away from the band
        // solver.eig(matrix.data(), matrix.dimension(0), 'V', 1, 1, vl, vu);

        // if(solver.result.meta.eigvals_found and solver.result.meta.eigvecsR_found) {
        //     // tools::log->info("optimize_variance_eig_executor: vl {:.3e} | vu {:.3e}", vl, vu);
        //     tools::log->info("Found {} eigvals ({} converged)", solver.result.meta.nev, solver.result.meta.nev_converged);
        //     auto eigvals = eig::view::get_eigvals<real>(solver.result, false);
        //     auto indices = num::range<long>(0l, eigvals.size());
        //     auto eigComp = EigIdxComparator(meta.optRitz, 0, eigvals.data(), eigvals.size());
        //     std::sort(indices.begin(), indices.end(), eigComp); // Should sort them according to distance from eigval
        //     indices.resize(std::min(eigvals.size(), 10l));      // We only need the first few indices, say 4
        //     for(auto idx : indices) { tools::log->info(" -- idx {}: {:.16f}", idx, eigvals(idx)); }
        // }

        // eig::solver solver1, solver2;
        // auto        H1 = tensors.get_effective_hamiltonian<Scalar>();
        // auto        H2 = tensors.get_effective_hamiltonian_squared<Scalar>();
        // solver1.eig<eig::Form::SYMM>(H1.data(), H1.dimension(0));
        // solver2.eig<eig::Form::SYMM>(H2.data(), H2.dimension(0));
        // auto evals1 = eig::view::get_eigvals<real>(solver1.result);
        // auto evals2 = eig::view::get_eigvals<real>(solver2.result);
        // for(long idx = 0; idx < std::min(evals1.size(), evals2.size()); ++idx) {
        // fmt::print("idx {:2}: H {:20.16f}  H² {:20.16f}\n", idx, evals1[idx], evals2[idx]);
        // }
        // fmt::print("\n");
        // }
    }

    opt_mps optimize_folded_spectrum_eig(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status,
                                         OptMeta &meta) {
        if(meta.optSolver == OptSolver::EIGS) return optimize_folded_spectrum(tensors, initial_mps, status, meta);

        initial_mps.validate_initial_mps();

        const auto problem_size = tensors.active_problem_size();
        if(problem_size > settings::precision::eig_max_size)
            throw except::logic_error("optimize_folded_spectrum_eig: the problem size is too large for eig: {}", problem_size);

        tools::log->debug("optimize_folded_spectrum_eig: ritz {} | type {} | algo {}", enum2sv(meta.optRitz), enum2sv(meta.optType), enum2sv(meta.optAlgo));

        reports::eigs_add_entry(initial_mps, spdlog::level::debug);
        auto                 t_var = tid::tic_scope("eig-xdmrg", tid::level::higher);
        std::vector<opt_mps> results;
        switch(meta.optType) {
            case OptType::REAL: optimize_folded_spectrum_eig_executor<real>(tensors, initial_mps, results, meta); break;
            case OptType::CPLX: optimize_folded_spectrum_eig_executor<cplx>(tensors, initial_mps, results, meta); break;
        }
        auto t_post = tid::tic_scope("post");
        if(results.empty()) {
            meta.optExit = OptExit::FAIL_ERROR;
            return initial_mps; // The solver failed
        }

        if(results.size() >= 2) {
            std::sort(results.begin(), results.end(), Comparator(meta)); // Smallest eigenvalue (i.e. variance) wins
        }

        for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::debug);
        return results.front();
    }
}
