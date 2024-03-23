#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/cast.h"
#include "math/eig.h"
#include "math/eig/matvec/matvec_mpo.h"
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
            case OptRitz::SR: solver.eig(matrix.data(), matrix.dimension(0), 'I', SR_il, SR_iu, 0.0, 1.0); break;
            case OptRitz::LR: solver.eig(matrix.data(), matrix.dimension(0), 'I', LR_il, LR_iu, 0.0, 1.0); break;
            case OptRitz::SM:
            case OptRitz::CR: {
                auto eigval = initial_mps.get_eigval();
                auto eigvar = initial_mps.get_variance();
                auto vl     = eigval - 2 * eigvar;
                auto vu     = eigval + 2 * eigvar;
                solver.eig(matrix.data(), matrix.dimension(0), 'V', 1, 1, vl, vu);
                break;
            }
        }
        extract_results(tensors, initial_mps, meta, solver, results, false);
    }

    opt_mps optimize_energy_eig(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status, OptMeta &meta) {
        if(meta.optFunc != OptFunc::ENERGY)
            throw except::logic_error("optimize_energy_eig: Expected OptFunc [{}] | Got [{}]", enum2sv(OptFunc::ENERGY), enum2sv(meta.optFunc));
        if(meta.optRitz == OptRitz::SM and not tensors.model->has_energy_shifted_mpo())
            throw std::logic_error("optimize_energy_eig: Ritz [SM] requires energy-shifted MPO ");

        const auto problem_size = tensors.active_problem_size();
        if(problem_size > settings::solver::eig_max_size)
            throw except::logic_error("optimize_energy_eig: the problem size is too large for eig: {}", problem_size);

        tools::log->trace("Full diagonalization of H");
        tools::log->debug("optimize_energy_eig: ritz {} | type {} | func {} | algo {}", enum2sv(meta.optRitz), enum2sv(meta.optType), enum2sv(meta.optFunc),
                          enum2sv(meta.optAlgo));

        auto                 t_gs = tid::tic_scope("eig-ene", tid::level::higher);
        std::vector<opt_mps> results;
        switch(meta.optType) {
            case OptType::REAL: optimize_energy_eig_executor<real>(tensors, initial_mps, results, meta); break;
            case OptType::CPLX: optimize_energy_eig_executor<cplx>(tensors, initial_mps, results, meta); break;
        }
        if(results.empty()) {
            meta.optExit = OptExit::FAIL_ERROR;
            return initial_mps; // The solver failed
        }
        // Smallest energy wins (because they are shifted!)
        if(results.size() >= 2) std::sort(results.begin(), results.end(), Comparator(meta, initial_mps.get_energy()));
        for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::debug);
        return results.front();
    }

}
