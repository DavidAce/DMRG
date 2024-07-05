#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/eig.h"
#include "math/eig/matvec/matvec_mpo.h"
#include "math/eig/matvec/matvec_mpos.h"
#include "math/linalg/matrix.h"
#include "math/tenx.h"
#include "tensors/edges/EdgesFinite.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/contraction.h"
#include "tools/common/log.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
#include <primme/primme.h>

namespace tools::finite::opt {

    template<typename Scalar>
    void preconditioner(void *x, int *ldx, void *y, int *ldy, int *blockSize, primme_params *primme, int *ierr) {
        if(x == nullptr) return;
        if(y == nullptr) return;
        if(primme == nullptr) return;
        const auto H_ptr = static_cast<MatVecMPO<Scalar> *>(primme->matrix);
        H_ptr->FactorOP();
        H_ptr->MultOPv(x, ldx, y, ldy, blockSize, primme, ierr);
    }

    template<typename Scalar>
    std::vector<opt_mps> eigs_energy_executor(const TensorsFinite &tensors, const opt_mps &initial_mps, const OptMeta &meta) {
        if(meta.optCost != OptCost::ENERGY)
            throw except::runtime_error("eigs_energy_executor: Expected OptCost [{}] | Got [{}]", enum2sv(OptCost::ENERGY), enum2sv(meta.optCost));

        eig::Ritz ritz = eig::stringToRitz(enum2sv(meta.optRitz));
        tools::log->trace("Defining Hamiltonian matrix-vector product");
        tools::log->trace("Defining eigenvalue solver");
        eig::solver solver;

        const auto       &mpo  = tensors.get_multisite_mpo();
        const auto       &env  = tensors.get_multisite_env_ene_blk();
        const auto        size = initial_mps.get_tensor().size();
        MatVecMPO<Scalar> hamiltonian(env.L, env.R, mpo);
        hamiltonian.factorization = eig::Factorization::NONE; // No LU factorization by default

        // https://www.cs.wm.edu/~andreas/software/doc/appendix.html#c.primme_params.eps
        solver.config.tol             = meta.eigs_tol.value_or(settings::precision::eigs_tol_min);
        solver.config.compute_eigvecs = eig::Vecs::ON;
        solver.config.lib             = eig::Lib::PRIMME;
        solver.config.primme_method   = eig::PrimmeMethod::PRIMME_DYNAMIC;
        solver.config.maxIter         = meta.eigs_iter_max.value_or(settings::precision::eigs_iter_max);
        solver.config.maxTime         = 2 * 60 * 60;
        solver.config.maxNev          = meta.eigs_nev.value_or(1);
        solver.config.maxNcv          = meta.eigs_ncv.value_or(settings::precision::eigs_ncv); // arpack needs ncv ~512. Primme seems happy with 4-32.
        solver.config.primme_locking  = true;
        solver.config.loglevel        = 1;

        Eigen::Tensor<Scalar, 3> init;
        if constexpr(std::is_same_v<Scalar, double>) {
            init = initial_mps.get_tensor().real();
            if(ritz == eig::Ritz::SR) ritz = eig::Ritz::SA;
            if(ritz == eig::Ritz::LR) ritz = eig::Ritz::LA;
        } else {
            init = initial_mps.get_tensor();
        }
        solver.config.ritz = ritz;
        solver.config.initial_guess.push_back({init.data(), 0});
        tools::log->trace("Finding energy eigenstate {}", enum2sv(meta.optRitz));

        if(meta.optRitz == OptRitz::SM) {
            // This is used to find excited eigenstates of H
            if(size <= settings::precision::eigs_max_size_shift_invert and meta.optAlgo == OptAlgo::SHIFTINV) {
                hamiltonian.factorization       = eig::Factorization::LU;
                solver.config.shift_invert      = eig::Shinv::ON;
                solver.config.ritz              = eig::Ritz::primme_largest_abs;
                solver.config.sigma             = cplx(meta.eigv_target.value_or(0.0), 0.0);
                solver.config.primme_projection = "primme_proj_default";
                solver.config.primme_locking    = false;
            } else {
                if(meta.optAlgo == OptAlgo::SHIFTINV)
                    tools::log->debug("eigs_energy_executor: ignoring OptAlgo::SHIFTINV: "
                                      "problem size {} > settings::precision::eigs_max_size_shift_invert ({})",
                                      size, settings::precision::eigs_max_size_shift_invert);
                // The problem size is too big. Just find the eigenstate closest to zero
                // Remember: since the Hamiltonian is expected to be shifted, the target energy is near zero.
                hamiltonian.factorization         = eig::Factorization::NONE;
                solver.config.ritz                = eig::Ritz::primme_closest_abs;
                solver.config.primme_projection   = "primme_proj_refined";
                solver.config.primme_targetShifts = {meta.eigv_target.value_or(0.0)};
                // solver.config.primme_preconditioner = preconditioner<Scalar>;
            }
        }

        solver.eigs(hamiltonian);

        std::vector<opt_mps> results;
        internal::extract_results(tensors, initial_mps, meta, solver, results, false);
        return results;
    }

    template<typename Scalar>
    std::vector<opt_mps> eigs_energy_executor_large(const TensorsFinite &tensors, const opt_mps &initial_mps, const OptMeta &meta) {
        if(meta.optCost != OptCost::ENERGY)
            throw except::runtime_error("eigs_energy_executor: Expected OptCost [{}] | Got [{}]", enum2sv(OptCost::ENERGY), enum2sv(meta.optCost));

        eig::Ritz ritz = eig::stringToRitz(enum2sv(meta.optRitz));
        tools::log->trace("Defining Hamiltonian matrix-vector product");
        tools::log->trace("Defining eigenvalue solver");
        eig::solver solver;
        const auto &model       = *tensors.model;
        const auto &edges       = *tensors.edges;
        const auto &envs        = edges.get_ene_active();
        const auto &mpos        = model.get_mpo_active();
        auto        hamiltonian = MatVecMPOS<Scalar>(mpos, envs);
        auto        size        = hamiltonian.rows();
        // hamiltonian.factorization = eig::Factorization::NONE; // No LU factorization by default

        // https://www.cs.wm.edu/~andreas/software/doc/appendix.html#c.primme_params.eps
        solver.config.tol             = meta.eigs_tol.value_or(settings::precision::eigs_tol_min);
        solver.config.compute_eigvecs = eig::Vecs::ON;
        solver.config.lib             = eig::Lib::PRIMME;
        solver.config.primme_method   = eig::stringToMethod(meta.primme_method.value_or("PRIMME_DEFAULT_MIN_MATVECS"));
        solver.config.maxIter         = meta.eigs_iter_max.value_or(settings::precision::eigs_iter_max);
        solver.config.maxTime         = 2 * 60 * 60;
        solver.config.maxNev          = meta.eigs_nev.value_or(1);
        solver.config.maxNcv          = meta.eigs_ncv.value_or(settings::precision::eigs_ncv); // arpack needs ncv ~512. Primme seems happy with 4-32.
        solver.config.primme_locking  = true;
        solver.config.loglevel        = 1;

        Eigen::Tensor<Scalar, 3> init;
        if constexpr(std::is_same_v<Scalar, double>) {
            init = initial_mps.get_tensor().real();
            if(ritz == eig::Ritz::SR) ritz = eig::Ritz::SA;
            if(ritz == eig::Ritz::LR) ritz = eig::Ritz::LA;
        } else {
            init = initial_mps.get_tensor();
        }
        solver.config.ritz = ritz;
        solver.config.initial_guess.push_back({init.data(), 0});
        tools::log->trace("Finding energy eigenstate {}", enum2sv(meta.optRitz));

        if(meta.optRitz == OptRitz::SM) {
            // This is used to find excited eigenstates of H
            if(meta.optAlgo == OptAlgo::SHIFTINV)
                tools::log->debug("eigs_energy_executor: ignoring OptAlgo::SHIFTINV: "
                                  "problem size {} > settings::precision::eigs_max_size_shift_invert ({})",
                                  size, settings::precision::eigs_max_size_shift_invert);
            // The problem size is too big. Just find the eigenstate closest to zero
            // Remember: since the Hamiltonian is expected to be shifted, the target energy is near zero.
            // hamiltonian.factorization         = eig::Factorization::NONE;
            solver.config.ritz                = eig::Ritz::primme_closest_abs;
            solver.config.primme_projection   = "primme_proj_refined";
            solver.config.primme_targetShifts = {meta.eigv_target.value_or(0.0)};
        }

        solver.eigs(hamiltonian);

        std::vector<opt_mps> results;
        internal::extract_results(tensors, initial_mps, meta, solver, results, false);
        return results;
    }

    opt_mps internal::optimize_energy_eigs(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status, OptMeta &meta) {
        if(meta.optSolver == OptSolver::EIG) return optimize_energy_eig(tensors, initial_mps, status, meta);
        tools::log->debug("optimize_energy_eigs: ritz {} | type {} | func {} | algo {}", enum2sv(meta.optRitz), enum2sv(meta.optType), enum2sv(meta.optCost),
                          enum2sv(meta.optAlgo));

        initial_mps.validate_initial_mps();
        // if(not tensors.model->is_shifted()) throw std::runtime_error("optimize_variance_eigs requires energy-shifted MPO²");
        reports::eigs_add_entry(initial_mps, spdlog::level::debug);

        auto                 t_eigs = tid::tic_scope("eigs-ene", tid::higher);
        std::vector<opt_mps> results;
        if(meta.optType == OptType::REAL) results = eigs_energy_executor_large<real>(tensors, initial_mps, meta);
        if(meta.optType == OptType::CPLX) results = eigs_energy_executor_large<cplx>(tensors, initial_mps, meta);

        auto t_post = tid::tic_scope("post");
        if(results.empty()) {
            meta.optExit = OptExit::FAIL_ERROR;
            return initial_mps; // The solver failed
        }

        auto comparator = [&meta, &initial_mps](const opt_mps &lhs, const opt_mps &rhs) {
            auto diff = std::abs(lhs.get_eshift_eigval() - rhs.get_eshift_eigval());
            if(diff < settings::precision::eigs_tol_min) return lhs.get_overlap() > rhs.get_overlap();
            switch(meta.optRitz) {
                case OptRitz::SR: return lhs.get_energy() < rhs.get_energy();
                case OptRitz::LR: return lhs.get_energy() > rhs.get_energy();
                case OptRitz::SM: {
                    // return std::abs(lhs.get_eigs_eigval()) < std::abs(rhs.get_eigs_eigval());
                    auto diff_energy_lhs = std::abs(lhs.get_energy() - initial_mps.get_energy());
                    auto diff_energy_rhs = std::abs(rhs.get_energy() - initial_mps.get_energy());
                    return diff_energy_lhs < diff_energy_rhs;
                }
                default: throw except::runtime_error("Ground state optimization with ritz {} is not implemented", enum2sv(meta.optRitz));
            }
        };

        if(results.size() >= 2) std::sort(results.begin(), results.end(), comparator);
        for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::debug);
        return results.front();
    }

}
