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

    template<typename MatrixProductType>
    void simps_preconditioner(void *x, int *ldx, void *y, int *ldy, int *blockSize, [[maybe_unused]] primme_params *primme, int *ierr) {
        if(x == nullptr) return;
        if(y == nullptr) return;
        if(primme == nullptr) return;
        using T               = typename MatrixProductType::Scalar;
        const auto  H_ptr     = static_cast<MatrixProductType *>(primme->matrix);
        const auto  shape_mps = H_ptr->get_shape_mps();
        const auto &mpo       = H_ptr->get_mpo();
        const auto &envL      = H_ptr->get_envL();
        const auto &envR      = H_ptr->get_envR();
        for(int i = 0; i < *blockSize; i++) {
            auto mps_in  = Eigen::TensorMap<const Eigen::Tensor<T, 3>>(static_cast<T *>(x) + *ldx * i, shape_mps);
            auto mps_out = Eigen::TensorMap<Eigen::Tensor<T, 3>>(static_cast<T *>(y) + *ldy * i, shape_mps);
            tools::common::contraction::matrix_inverse_vector_product(mps_out, mps_in, mpo, envL, envR);
        }
        *ierr = 0;
    }

    template<typename Scalar>
    std::vector<opt_mps> eigs_energy_executor(const TensorsFinite &tensors, const opt_mps &initial_mps, const OptMeta &meta) {
        if(meta.optMode != OptMode::ENERGY)
            throw except::runtime_error("Wrong optimization mode [{}]. Expected [{}]", enum2sv(meta.optMode), enum2sv(OptMode::ENERGY));
        if(meta.optRitz == OptRitz::SM and not tensors.model->is_shifted())
            throw std::runtime_error("eigs_optimize_energy with ritz [SM] requires energy-shifted MPO ");

        eig::Ritz ritz = eig::stringToRitz(enum2sv(meta.optRitz));

        tools::log->trace("Defining Hamiltonian matrix-vector product");
        tools::log->trace("Defining eigenvalue solver");
        eig::solver solver;

        const auto       &mpo = tensors.get_multisite_mpo();
        const auto       &env = tensors.get_multisite_env_ene_blk();
        MatVecMPO<Scalar> hamiltonian(env.L, env.R, mpo);
        // https://www.cs.wm.edu/~andreas/software/doc/appendix.html#c.primme_params.eps
        solver.config.tol             = 1e-10; // This is the target residual_norm norm. 1e-10 seems to be sufficient:
        solver.config.compress        = settings::precision::use_compressed_mpo_squared_otf;
        solver.config.compute_eigvecs = eig::Vecs::ON;
        solver.config.lib             = eig::Lib::PRIMME;
        solver.config.primme_method   = eig::PrimmeMethod::PRIMME_GD_Olsen_plusK;
        solver.config.maxIter         = 4000;
        solver.config.maxTime         = 60 * 60;
        solver.config.maxNev          = static_cast<eig::size_type>(1);
        solver.config.maxNcv          = static_cast<eig::size_type>(16); // arpack needs ncv ~512 to handle all cases. Primme seems content with 16.
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
            solver.config.maxIter              = 100000;
            solver.config.maxNev               = static_cast<eig::size_type>(8);
            solver.config.maxNcv               = static_cast<eig::size_type>(32);
            solver.config.primme_projection    = "primme_proj_refined";
            solver.config.primme_locking       = false;
            solver.config.primme_target_shifts = {initial_mps.get_eigval()};
            solver.config.ritz                 = eig::Ritz::primme_closest_abs;
            if(meta.optMode == OptMode::SIMPS) solver.config.primme_preconditioner = simps_preconditioner<MatVecMPO<Scalar>>;

            solver.eigs(hamiltonian);
        } else {
            // Since we use energy-shifted mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
            // would otherwise cause trouble for the eigenvalue solver. This equates to subtracting sigma * identity from the bottom corner of the mpo.
            // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
            //        solver.config.sigma = 1.0;
            solver.eigs(hamiltonian);
        }

        std::vector<opt_mps> results;
        internal::eigs_extract_results(tensors, initial_mps, meta, solver, results, false);

        auto comparator = [&ritz, &meta, &initial_mps](const opt_mps &lhs, const opt_mps &rhs) {
            auto diff = std::abs(lhs.get_eigval() - rhs.get_eigval());
            if(diff < settings::precision::eigs_tolerance) return lhs.get_overlap() > rhs.get_overlap();
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
        for(const auto &mps : results) reports::eigs_add_entry(mps, spdlog::level::info);
        return results;
    }

    opt_mps internal::eigs_optimize_energy(const TensorsFinite &tensors, const opt_mps &initial_mps, const AlgorithmStatus &status, OptMeta &meta) {
        if(tensors.active_problem_size() <= settings::precision::max_size_full_diag) return eig_optimize_energy(tensors, initial_mps, status, meta);

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
