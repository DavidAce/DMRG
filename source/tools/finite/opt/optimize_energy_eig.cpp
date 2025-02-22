#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/cast.h"
#include "math/eig.h"
#include "math/eig/matvec/matvec_mpo.h"
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
    void optimize_energy_eig_executor(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        eig::solver solver;
        auto        matrix = tensors.get_effective_hamiltonian<Scalar>();
        int         SR_il  = 1; // min nev index (starts from 1)
        int         SR_iu  = 1; // max nev index
        int         LR_il  = safe_cast<int>(matrix.dimension(0));
        int         LR_iu  = safe_cast<int>(matrix.dimension(0));
        switch(meta.optRitz) {
            case OptRitz::NONE: throw std::logic_error("optimize_energy_eig_executor: Invalid: OptRitz::NONE");
            case OptRitz::SR: solver.eig(matrix.data(), matrix.dimension(0), 'I', SR_il, SR_iu, 0.0, 1.0); break;
            case OptRitz::LR: solver.eig(matrix.data(), matrix.dimension(0), 'I', LR_il, LR_iu, 0.0, 1.0); break;
            case OptRitz::LM: solver.eig(matrix.data(), matrix.dimension(0)); break; // Find all eigenvalues
            case OptRitz::SM:
            case OptRitz::TE:
            case OptRitz::IS: {
                // Find all eigenvalues within a thin energy band
                auto eigval = initial_mps.get_energy(); // The current energy
                auto eigvar = initial_mps.get_variance();
                auto eshift = initial_mps.get_eshift();                          // The energy shift is our target energy for excited states
                auto vl     = eshift - std::abs(eigval) - 2 * std::sqrt(eigvar); // Find energies at most two sigma away from the band
                auto vu     = eshift + std::abs(eigval) + 2 * std::sqrt(eigvar); // Find energies at most two sigma away from the band
                solver.eig(matrix.data(), matrix.dimension(0), 'V', 1, 1, vl, vu);
                // tools::log->info("optimize_energy_eig_executor: vl {:.3e} | vu {:.3e}", vl, vu);
                // Filter the results
                if(solver.result.meta.eigvals_found and solver.result.meta.eigvecsR_found) {
                    // tools::log->info("Found {} eigvals ({} converged)", solver.result.meta.nev, solver.result.meta.nev_converged);
                    auto eigvals = eig::view::get_eigvals<fp64>(solver.result, false);
                    auto indices = num::range<long>(0l, eigvals.size());
                    auto eigComp = EigIdxComparator(meta.optRitz, eigval, eigvals.data(), eigvals.size());
                    std::sort(indices.begin(), indices.end(), eigComp); // Should sort them according to distance from eigval
                    indices.resize(safe_cast<size_t>(std::min(eigvals.size(), 10l)));      // We only need the first few indices, say 4
                    // for(auto idx : indices) { tools::log->info(" -- idx {}: {:.16f}", idx, eigvals(idx)); }
                    extract_results(tensors, initial_mps, meta, solver, results, false, indices);
                }
                return;
            }
        }

        extract_results(tensors, initial_mps, meta, solver, results, false);
    }

    opt_mps optimize_energy_eig(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status, OptMeta &meta) {
        if(meta.optAlgo != OptAlgo::DMRG)
            throw except::logic_error("optimize_energy_eig: Expected OptAlgo [{}] | Got [{}]", enum2sv(OptAlgo::DMRG), enum2sv(meta.optAlgo));

        const auto problem_size = tensors.active_problem_size();
        if(problem_size > settings::precision::eig_max_size)
            throw except::logic_error("optimize_energy_eig: the problem size is too large for eig: {} > {}(max)", problem_size,
                                      settings::precision::eig_max_size);

        tools::log->debug("optimize_energy_eig: ritz {} | type {} | algo {}", enum2sv(meta.optRitz), enum2sv(meta.optType), enum2sv(meta.optAlgo));

        initial_mps.validate_initial_mps();
        // if(not tensors.model->is_shifted()) throw std::runtime_error("optimize_variance_eigs requires energy-shifted MPO²");
        reports::eigs_add_entry(initial_mps, spdlog::level::debug);

        auto                 t_gs = tid::tic_scope("eig-ene", tid::level::higher);
        std::vector<opt_mps> results;
        switch(meta.optType) {
            case OptType::REAL: optimize_energy_eig_executor<fp64>(tensors, initial_mps, results, meta); break;
            case OptType::CPLX: optimize_energy_eig_executor<cx64>(tensors, initial_mps, results, meta); break;
        }
        if(results.empty()) {
            meta.optExit = OptExit::FAIL_ERROR;
            return initial_mps; // The solver failed
        }
        // Smallest energy wins (because they are shifted!)
        tools::log->debug("Sorting eigenpairs | initial energy {}", initial_mps.get_energy());
        if(results.size() >= 2) std::sort(results.begin(), results.end(), Comparator(meta, initial_mps.get_energy()));
        for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::debug);
        return results.front();
    }

}
