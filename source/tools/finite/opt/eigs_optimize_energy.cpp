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
#include <math/eig/matvec/matvec_mpo.h>
#include <tensors/TensorsFinite.h>
#include <tid/tid.h>
#include <tools/common/log.h>
#include <tools/finite/measure.h>

template<typename Scalar>
std::vector<tools::finite::opt::opt_mps> solve(const TensorsFinite &tensors, const tools::finite::opt::opt_mps &initial_mps,
                                               const tools::finite::opt::OptMeta &meta) {
    using namespace tools::finite::opt;
    auto      problem_size = initial_mps.get_tensor().size();
    eig::Ritz ritz         = eig::stringToRitz(enum2sv(meta.optRitz));

    tools::log->trace("Defining Hamiltonian matrix-vector product");
    tools::log->trace("Defining eigenvalue solver");
    eig::solver solver;
    if(problem_size < settings::precision::max_size_full_diag) {
        tools::log->trace("Finding ground state");
        const auto &matrix = tensors.template get_effective_hamiltonian<Scalar>();
        solver.config.tag  = std::is_same_v<double, Scalar> ? "dsyevd" : "zheevd";
        solver.eig<eig::Form::SYMM>(matrix.data(), matrix.dimension(0));
    } else {
        const auto       &mpo = tensors.get_multisite_mpo();
        const auto       &env = tensors.get_multisite_env_ene_blk();
        MatVecMPO<Scalar> matrix(env.L, env.R, mpo);
        // https://www.cs.wm.edu/~andreas/software/doc/appendix.html#c.primme_params.eps
        solver.config.tol             = 1e-10; // This is the target residual norm. 1e-10 seems to be sufficient:
        solver.config.compress        = settings::precision::use_compressed_mpo_squared_otf;
        solver.config.compute_eigvecs = eig::Vecs::ON;
        solver.config.lib             = eig::Lib::PRIMME;
        solver.config.tag             = "primme";
        solver.config.primme_method   = eig::PrimmeMethod::PRIMME_GD_Olsen_plusK;
        solver.config.maxIter         = 80000;
        solver.config.maxTime         = 60 * 60;
        solver.config.maxNev          = static_cast<eig::size_type>(1);
        solver.config.maxNcv          = static_cast<eig::size_type>(16); // arpack needs ncv ~512 to handle all cases. Primme seems content with 16.
        solver.config.primme_locking  = false;
        // Since we use energy-shifted mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
        // would otherwise cause trouble for the eigenvalue solver. This equates to subtracting sigma * identity from the bottom corner of the mpo.
        // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
        solver.config.sigma    = 1.0;
        solver.config.loglevel = 2;

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
    tools::finite::opt::internal::eigs_extract_results(tensors, initial_mps, meta, solver, results, true);

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
    for(const auto &mps : results) reports::eigs_add_entry(mps);
    return results;
}

tools::finite::opt::opt_mps tools::finite::opt::internal::eigs_optimize_energy(const TensorsFinite &tensors, const opt_mps &initial_mps,
                                                                               [[maybe_unused]] const AlgorithmStatus &status, OptMeta &meta) {
    tools::log->debug("Ground state optimization with ritz {} | type {}", enum2sv(meta.optRitz), enum2sv(meta.optType));
    auto                 t_gs = tid::tic_scope("gs");
    std::vector<opt_mps> results;
    if(meta.optType == OptType::REAL) results = solve<real>(tensors, initial_mps, meta);
    if(meta.optType == OptType::CPLX) results = solve<cplx>(tensors, initial_mps, meta);

    reports::print_eigs_report();
    if(results.empty())
        return initial_mps; // Solver failed
    else
        return results.front();
}
