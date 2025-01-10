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
#include <queue>
namespace tools::finite::opt::internal {
    template<typename VecType>
    std::vector<long> get_k_largest(const VecType &vec, size_t k) {
        using Scalar   = typename VecType::Scalar;
        using idx_pair = std::pair<Scalar, long>;
        std::priority_queue<idx_pair, std::vector<idx_pair>, std::greater<idx_pair>> q;
        for(long i = 0; i < vec.size(); ++i) {
            if(q.size() < k)
                q.emplace(vec[i], i);
            else if(q.top().first < vec[i]) {
                q.pop();
                q.emplace(vec[i], i);
            }
        }
        k = q.size();
        std::vector<long> res(k);
        for(size_t i = 0; i < k; ++i) {
            res[k - i - 1] = q.top().second;
            q.pop();
        }
        return res;
    }
    template<typename VecType>
    std::vector<long> get_k_smallest(const VecType &vec, size_t k) {
        using Scalar = typename VecType::Scalar;
        std::priority_queue<Scalar> pq;
        for(auto d : vec) {
            if(pq.size() >= k && pq.top() > d) {
                pq.push(d);
                pq.pop();
            } else if(pq.size() < k) {
                pq.push(d);
            }
        }
        Scalar            kth_element = pq.top();
        std::vector<long> result;
        for(long i = 0; i < vec.size(); i++)
            if(vec[i] <= kth_element) { result.emplace_back(i); }
        return result;
    }

    template<typename Scalar>
    void optimize_generalized_shift_invert_eig_executor(const TensorsFinite &tensors, const opt_mps &initial_mps, std::vector<opt_mps> &results,
                                                        OptMeta &meta) {
        // Solve the generalized problem Hx=(1/E)H²x

        if(meta.optRitz == OptRitz::NONE) return;
        eig::solver solver;
        auto        matrixA = tensors.get_effective_hamiltonian<Scalar>();
        auto        matrixB = tensors.get_effective_hamiltonian_squared<Scalar>();
        auto        nev     = std::min<int>(static_cast<int>(matrixA.dimension(0)), meta.eigs_nev.value_or(1));
        switch(meta.optRitz) {
            case OptRitz::SR: {
                auto il = 1;
                auto iu = nev;
                solver.eig(matrixA.data(), matrixB.data(), matrixA.dimension(0), 'I', il, iu, 0.0, 1.0);
                extract_results(tensors, initial_mps, meta, solver, results, true);
                break;
            }
            case OptRitz::LR: {
                auto il = static_cast<int>(matrixA.dimension(0) - (nev - 1));
                auto iu = static_cast<int>(matrixA.dimension(0));
                solver.eig(matrixA.data(), matrixB.data(), matrixA.dimension(0), 'I', il, iu, 0.0, 1.0);
                extract_results(tensors, initial_mps, meta, solver, results, true);
                break;
            }
            case OptRitz::SM: {
                solver.eig(matrixA.data(), matrixB.data(), matrixA.dimension(0));
                // Determine the nev largest eigenvalues
                auto eigvals = eig::view::get_eigvals<fp64>(solver.result, false);
                auto indices = std::vector<long>();
                if(nev == 1) {
                    long                  maxidx = 0;
                    [[maybe_unused]] auto maxval = eigvals.cwiseAbs().minCoeff(&maxidx);
                    indices                      = {maxidx};
                } else {
                    indices = get_k_smallest(eigvals.cwiseAbs(), safe_cast<size_t>(nev));
                }
                extract_results(tensors, initial_mps, meta, solver, results, true, indices);
                break;
            }

            case OptRitz::LM: {
                solver.eig(matrixA.data(), matrixB.data(), matrixA.dimension(0));
                // Determine the nev largest eigenvalues
                auto eigvals = eig::view::get_eigvals<fp64>(solver.result, false);
                auto indices = std::vector<long>();
                if(nev == 1) {
                    long                  maxidx = 0;
                    [[maybe_unused]] auto maxval = eigvals.cwiseAbs().maxCoeff(&maxidx);
                    indices                      = {maxidx};
                } else {
                    indices = get_k_largest(eigvals.cwiseAbs(), safe_cast<size_t>(nev));
                }
                extract_results(tensors, initial_mps, meta, solver, results, true, indices);
                break;
                //
                // eig::solver solver_u, solver_l;
                // auto        matrixA_u = tensors.get_effective_hamiltonian<Scalar>();
                // auto        matrixB_u = tensors.get_effective_hamiltonian_squared<Scalar>();
                // auto        matrixA_l = matrixA_u;
                // auto        matrixB_l = matrixB_u;
                // auto        nev       = std::min<int>(static_cast<int>(matrixA_u.dimension(0)), meta.eigs_nev.value_or(1));
                // auto        il        = static_cast<int>(matrixA_u.dimension(0) - (nev - 1));
                // auto        iu        = static_cast<int>(matrixA_u.dimension(0));
                // solver_u.eig(matrixA_u.data(), matrixB_u.data(), matrixA_u.dimension(0), 'I', il, iu, 0.0, 1.0);
                // solver_l.eig(matrixA_l.data(), matrixB_l.data(), matrixA_l.dimension(0), 'I', 1, nev, 0.0, 1.0);
                // auto eigvals_u = eig::view::get_eigvals<fp64>(solver_u.result, false);
                // auto eigvals_l = eig::view::get_eigvals<fp64>(solver_l.result, false);
                // auto ev_u      = eigvals_u(eigvals_u.size() - 1);
                // auto ev_l      = eigvals_l(0);
                // if(std::abs(ev_l) >= std::abs(ev_u)) {
                //     extract_results(tensors, initial_mps, meta, solver_l, results, true);
                // } else {
                //     extract_results(tensors, initial_mps, meta, solver_u, results, true);
                // }
                // break;
            }
            default: {
            }
        }
    }

    opt_mps optimize_generalized_shift_invert_eig(const TensorsFinite &tensors, const opt_mps &initial_mps, [[maybe_unused]] const AlgorithmStatus &status,
                                                  OptMeta &meta) {
        if(meta.optSolver == OptSolver::EIGS) return optimize_generalized_shift_invert(tensors, initial_mps, status, meta);
        initial_mps.validate_initial_mps();
        const auto problem_size = tensors.active_problem_size();
        if(problem_size > settings::precision::eig_max_size)
            throw except::logic_error("optimize_generalized_shift_invert_eig: the problem size is too large for eig: {}", problem_size);

        tools::log->debug("optimize_generalized_shift_invert_eig: ritz {} | type {} | algo {}", enum2sv(meta.optRitz), enum2sv(meta.optType),
                          enum2sv(meta.optAlgo));

        reports::eigs_add_entry(initial_mps, spdlog::level::debug);
        auto                 t_var = tid::tic_scope("eig-gdmrg", tid::level::higher);
        std::vector<opt_mps> results;
        switch(meta.optType) {
            case OptType::REAL: optimize_generalized_shift_invert_eig_executor<fp64>(tensors, initial_mps, results, meta); break;
            case OptType::CPLX: optimize_generalized_shift_invert_eig_executor<cx64>(tensors, initial_mps, results, meta); break;
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
