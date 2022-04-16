#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
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

namespace tools::finite::opt {
    template<typename Scalar>
    void internal::eig_executor(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        const auto problem_size = tensors.active_problem_size();
        if(problem_size > settings::precision::max_size_full_diag)
            throw except::logic_error("eig_executor: the Hamiltonian is too large for eig: {}", problem_size);
        eig::solver solver;
        if(meta.optMode == OptMode::VARIANCE) {
            tools::log->trace("Full diagonalization of local (H-E)²");
            const auto &matrix = tensors.get_effective_hamiltonian_squared<Scalar>();
            solver.eig(matrix.data(), matrix.dimension(0), 'I', 1, 1, 0.0, 1.0, 1);
        } else if(meta.optMode == OptMode::ENERGY) {
            tools::log->trace("Full diagonalization of local H");
            const auto &matrix = tensors.get_effective_hamiltonian<Scalar>();
            int         SR_il  = 1; // min nev index (starts from 1)
            int         SR_iu  = 1; // max nev index
            int         LR_il  = static_cast<int>(matrix.dimension(0));
            int         LR_iu  = static_cast<int>(matrix.dimension(0));
            switch(meta.optRitz) {
                case OptRitz::SR: solver.eig(matrix.data(), matrix.dimension(0), 'I', SR_il, SR_iu, 0.0, 1.0, 1); break;
                case OptRitz::LR: solver.eig(matrix.data(), matrix.dimension(0), 'I', LR_il, LR_iu, 0.0, 1.0, 1); break;
                case OptRitz::SM: {
                    auto eigval = initial_mps.get_eigval();
                    auto eigvar = initial_mps.get_variance();
                    auto vl     = eigval - 2 * eigvar;
                    auto vu     = eigval + 2 * eigvar;
                    solver.eig(matrix.data(), matrix.dimension(0), 'V', 1, 1, vl, vu, 10);
                    break;
                }
            }
        } else {
            tools::log->trace("Full diagonalization of local H");
            const auto &matrix = tensors.get_effective_hamiltonian<Scalar>();
            solver.eig(matrix.data(), matrix.dimension(0), eig::Vecs::ON);
        }
        tools::finite::opt::internal::eigs_extract_results(tensors, initial_mps, meta, solver, results, false);
    }
    template void internal::eig_executor<real>(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta);
    template void internal::eig_executor<cplx>(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta);

    opt_mps internal::eig_optimize_energy(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status,
                                          OptMeta &meta) {
        tools::log->debug("Energy optimization with ritz {} | type {}", enum2sv(meta.optRitz), enum2sv(meta.optType));
        if(meta.optMode != OptMode::ENERGY)
            throw except::runtime_error("Wrong optimization mode [{}]. Expected [{}]", enum2sv(meta.optMode), enum2sv(OptMode::ENERGY));
        if(meta.optRitz == OptRitz::SM and not tensors.model->is_shifted())
            throw std::runtime_error("eig_optimize_energy with ritz [SM] requires energy-shifted MPO ");

        auto                 t_gs = tid::tic_scope("eig-ene");
        std::vector<opt_mps> results;
        switch(meta.optType) {
            case OptType::REAL: eig_executor<real>(tensors, initial_mps, results, meta); break;
            case OptType::CPLX: eig_executor<cplx>(tensors, initial_mps, results, meta); break;
        }

        if(results.size() >= 2) std::sort(results.begin(), results.end(), Comparator(meta, initial_mps.get_energy()));
        for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::info);
        reports::print_eigs_report();
        if(results.empty())
            return initial_mps; // Solver failed
        else
            return results.front();
    }

}
