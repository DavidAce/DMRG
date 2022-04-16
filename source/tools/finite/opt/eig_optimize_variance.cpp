#include "algorithms/AlgorithmStatus.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "math/eig.h"
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
    void eig_variance_executor(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results, const OptMeta &meta) {
        const auto problem_size = tensors.active_problem_size();
        if(problem_size <= settings::precision::max_size_full_diag) {
            tools::log->trace("Full diagonalization of (H-E)²");
            auto        matrix = tensors.get_effective_hamiltonian_squared<Scalar>();
            eig::solver solver;
            solver.eig(matrix.data(), matrix.dimension(0), 'I', 1, 1, 0.0, 1.0, 1);
            eigs_extract_results(tensors, initial_mps, meta, solver, results, false);
        } else {
            throw except::logic_error("eig_variance_executor: Hamiltonian too big for eig solver: {}", problem_size);
        }
    }

    opt_mps eig_optimize_variance(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status, OptMeta &meta) {
        initial_mps.validate_basis_vector();
        if(not tensors.model->is_shifted()) throw std::runtime_error("eig_optimize_variance requires energy-shifted MPO²");
        reports::eigs_add_entry(initial_mps, spdlog::level::info);
        auto                 t_var = tid::tic_scope("eig-var");
        std::vector<opt_mps> results;
        switch(meta.optType) {
            case OptType::REAL: eig_variance_executor<real>(tensors, initial_mps, results, meta); break;
            case OptType::CPLX: eig_variance_executor<cplx>(tensors, initial_mps, results, meta); break;
        }
        auto t_post = tid::tic_scope("post");
        if(results.empty()) {
            meta.optExit = OptExit::FAIL_ERROR;
            return initial_mps; // The solver failed
        }

        if(results.size() >= 2) {
            std::sort(results.begin(), results.end(), Comparator(meta)); // Smallest eigenvalue (i.e. variance) wins
        }

        for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::info);
        return results.front();
    }
}
