#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/eig.h"
#include "math/eig/matvec/matvec_mpo.h"
#include "math/linalg/matrix.h"
#include "math/tenx.h"
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
        if(meta.optMode != OptMode::ENERGY and meta.optMode != OptMode::SIMPS)
            throw except::runtime_error("Wrong optimization mode [{}]. Expected [{}|{}]", enum2sv(meta.optMode), enum2sv(OptMode::ENERGY),
                                        enum2sv(OptMode::SIMPS));
        if(meta.optRitz == OptRitz::SM and not tensors.model->is_shifted())
            throw std::runtime_error("eigs_optimize_energy with ritz [SM] requires energy-shifted MPO ");

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
        solver.config.tol             = 1e-12;
        solver.config.compress        = settings::precision::use_compressed_mpo_squared_otf;
        solver.config.compute_eigvecs = eig::Vecs::ON;
        solver.config.lib             = eig::Lib::PRIMME;
        solver.config.primme_method   = eig::PrimmeMethod::PRIMME_DYNAMIC;
        solver.config.maxIter         = 10000;
        solver.config.maxTime         = 60 * 60;
        solver.config.maxNev          = static_cast<eig::size_type>(1);
        solver.config.maxNcv          = static_cast<eig::size_type>(4); // arpack needs ncv ~512 to handle all cases. Primme seems content with 4.
        solver.config.primme_locking  = false;
        solver.config.loglevel        = 2;

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
            if(size <= settings::solver::max_size_shift_invert) {
                hamiltonian.set_readyCompress(tensors.model->is_compressed_mpo_squared());
                hamiltonian.factorization       = eig::Factorization::LU;
                solver.config.shift_invert      = eig::Shinv::ON;
                solver.config.ritz              = eig::Ritz::primme_largest_abs;
                solver.config.sigma             = initial_mps.get_eigval();
                solver.config.maxIter           = 500;
                solver.config.maxNev            = static_cast<eig::size_type>(1);
                solver.config.maxNcv            = static_cast<eig::size_type>(4);
                solver.config.primme_projection = "primme_proj_default";
                solver.config.primme_locking    = false;
            } else {
                hamiltonian.factorization          = eig::Factorization::NONE;
                solver.config.maxIter              = 100000;
                solver.config.maxNev               = static_cast<eig::size_type>(1);
                solver.config.maxNcv               = static_cast<eig::size_type>(4);
                solver.config.primme_projection    = "primme_proj_default";
                solver.config.primme_locking       = false;
                solver.config.primme_target_shifts = {initial_mps.get_eigval()};
                solver.config.ritz                 = eig::Ritz::primme_closest_abs;
                if(meta.optMode == OptMode::SIMPS) solver.config.primme_preconditioner = preconditioner<Scalar>;
            }
        }

        if(meta.eigs_tol) solver.config.tol = meta.eigs_tol.value();
        if(meta.eigs_ncv) solver.config.maxNcv = meta.eigs_ncv.value();
        if(meta.eigs_iter_max) solver.config.maxIter = meta.eigs_iter_max.value();
        solver.eigs(hamiltonian);

        std::vector<opt_mps> results;
        internal::eigs_extract_results(tensors, initial_mps, meta, solver, results, false);

        auto comparator = [&ritz, &meta, &initial_mps](const opt_mps &lhs, const opt_mps &rhs) {
            auto diff = std::abs(lhs.get_eigval() - rhs.get_eigval());
            if(diff < settings::solver::eigs_tol_min) return lhs.get_overlap() > rhs.get_overlap();
            switch(ritz) {
                case eig::Ritz::SA:
                case eig::Ritz::SR: return lhs.get_energy() < rhs.get_energy();
                case eig::Ritz::LA:
                case eig::Ritz::LR: return lhs.get_energy() > rhs.get_energy();
                case eig::Ritz::SM:
                case eig::Ritz::primme_closest_abs: {
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
        return results;
    }

    opt_mps internal::eigs_optimize_energy(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status, OptMeta &meta) {
        if(tensors.active_problem_size() <= settings::solver::max_size_full_eigs) return eig_optimize_energy(tensors, initial_mps, status, meta);

        tools::log->debug("Energy optimization with ritz {} | type {}", enum2sv(meta.optRitz), enum2sv(meta.optType));
        auto                 t_eigs = tid::tic_scope("eigs-ene");
        std::vector<opt_mps> results;
        if(meta.optType == OptType::REAL) results = eigs_energy_executor<real>(tensors, initial_mps, meta);
        if(meta.optType == OptType::CPLX) results = eigs_energy_executor<cplx>(tensors, initial_mps, meta);

        reports::print_eigs_report();
        if(results.empty())
            return initial_mps; // Solver failed
        else
            return results.front();
    }

}
