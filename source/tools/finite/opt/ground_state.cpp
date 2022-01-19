#include "../opt_meta.h"
#include "../opt_mps.h"
#include "debug/exceptions.h"
#include "math/linalg/matrix.h"
#include "math/linalg/tensor.h"
#include "math/tenx.h"
#include "opt-internal.h"
#include "report.h"
#include <config/settings.h>
#include <general/iter.h>
#include <math/eig.h>
#include <math/eig/matvec/matvec_mps.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

template<typename Scalar>
std::vector<tools::finite::opt::opt_mps> solve(const tools::finite::opt::opt_mps &initial_mps, const TensorsFinite &tensors,
                                               const tools::finite::opt::OptMeta &meta) {
    using namespace tools::finite::opt;
    const auto &mpo = tensors.get_multisite_mpo();
    const auto &env = tensors.get_multisite_env_ene_blk();

    eig::Ritz ritz = eig::stringToRitz(enum2sv(meta.optRitz));

    tools::log->trace("Defining Hamiltonian matrix-vector product");
    MatVecMps<Scalar> matrix(env.L, env.R, mpo);
    tools::log->trace("Defining eigenvalue solver");
    eig::solver solver;
    if(matrix.rows() < settings::precision::max_size_full_diag) {
        tools::log->trace("Finding ground state");
        solver.config.tag = std::is_same_v<double, Scalar> ? "dsyevd" : "zheevd";
        solver.eig(matrix.get_tensor().data(), matrix.rows());
    } else {
        solver.config.tol             = 1e-14;
        solver.config.compress        = settings::precision::use_compressed_mpo_squared_otf;
        solver.config.compute_eigvecs = eig::Vecs::ON;
        solver.config.lib             = eig::Lib::PRIMME;
        solver.config.tag             = "primme";
        solver.config.primme_method   = eig::PrimmeMethod::PRIMME_GD_Olsen_plusK;
        solver.config.maxIter         = 80000;
        solver.config.maxTime         = 60 * 60;
        solver.config.maxNev          = static_cast<eig::size_type>(std::min(matrix.rows() - 1, 2));
        solver.config.maxNcv          = static_cast<eig::size_type>(512);
        solver.config.primme_locking  = false;
        // Since we use reduced energy mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
        // would otherwise cause trouble for the eigenvalue solver. This equates to subtracting sigma * identity from the bottom corner of the mpo.
        // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
        solver.config.sigma = 1.0;
        solver.setLogLevel(0);
        Eigen::Tensor<Scalar, 3> init;
        if constexpr(std::is_same_v<Scalar, double>) {
            init               = initial_mps.get_tensor().real();
            solver.config.ritz = ritz == eig::Ritz::SR ? eig::Ritz::SA : eig::Ritz::LA; // Adapt to real version
        } else {
            init               = initial_mps.get_tensor();
            solver.config.ritz = ritz;
        }

        solver.config.initial_guess.push_back({init.data(), 0});

        tools::log->trace("Finding ground state");
        solver.eigs(matrix);
    }

    std::vector<opt_mps> results;
    internal::krylov_extract_solutions(tensors, initial_mps, solver, results, meta, true);

    auto comparator = [&ritz, &meta](const opt_mps &lhs, const opt_mps &rhs) {
        auto diff = std::abs(lhs.get_eigval() - rhs.get_eigval());
        if(diff < settings::precision::eig_tolerance) return lhs.get_overlap() > rhs.get_overlap();
        switch(ritz) {
            case eig::Ritz::SA:
            case eig::Ritz::SR: return lhs.get_energy() < rhs.get_energy();
            case eig::Ritz::LA:
            case eig::Ritz::LR: return lhs.get_energy() > rhs.get_energy();
            default: throw std::runtime_error(fmt::format("Ground state optimization with ritz {} is not implemented", enum2sv(meta.optRitz)));
        }
    };

    if(results.size() >= 2) std::sort(results.begin(), results.end(), comparator);
    return results;
}

tools::finite::opt::opt_mps tools::finite::opt::internal::ground_state_optimization(const TensorsFinite &tensors, const AlgorithmStatus &status,
                                                                                    OptMeta &meta) {
    double  energy_reduced = 0.0;
    opt_mps initial_mps("current state", tensors.get_multisite_mps(), tensors.active_sites,
                        tools::finite::measure::energy(tensors) - energy_reduced, // Eigval
                        energy_reduced,                                           // Energy reduced for full system
                        tools::finite::measure::energy_variance(tensors),
                        1.0, // Overlap
                        tensors.get_length());

    return ground_state_optimization(initial_mps, tensors, status, meta);
}

tools::finite::opt::opt_mps tools::finite::opt::internal::ground_state_optimization(const opt_mps &initial_mps, const TensorsFinite &tensors,
                                                                                    [[maybe_unused]] const AlgorithmStatus &status, OptMeta &meta) {
    tools::log->debug("Ground state optimization with ritz {} | type {}", enum2sv(meta.optRitz), enum2sv(meta.optType));
    auto                 t_gs = tid::tic_scope("gs");
    std::vector<opt_mps> results;
    if(meta.optType == OptType::REAL) results = solve<real>(initial_mps, tensors, meta);
    if(meta.optType == OptType::CPLX) results = solve<cplx>(initial_mps, tensors, meta);

    constexpr size_t max_print = settings::debug ? 32 : 4;
    for(const auto &[num, mps] : iter::enumerate(results)) {
        if(num >= max_print) break;
        reports::krylov_add_entry(mps);
    }
    reports::print_krylov_report();

    if(results.empty())
        return initial_mps; // Solver failed
    else
        return results.front();
}
