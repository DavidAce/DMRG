#include "../opt_meta.h"
#include "../opt_mps.h"
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
                                                                                    [[maybe_unused]] const AlgorithmStatus &status,
                                                                                    OptMeta &meta) {
    tools::log->debug("Ground state optimization with ritz {} ...", enum2sv(meta.optRitz));
    using namespace internal;
    using namespace settings::precision;
    auto t_gs = tid::tic_scope("gs");

    eig::Ritz   ritz      = eig::stringToRitz(enum2sv(meta.optRitz));
    const auto &mpo       = tensors.get_multisite_mpo();
    const auto &env       = tensors.get_multisite_env_ene_blk();
    //    int nev = std::min(4l,(long)(tensors.state->active_problem_size()/2));
    tools::log->trace("Defining Hamiltonian matrix-vector product");
    MatVecMps<cplx> matrix(env.L, env.R, mpo);
    tools::log->trace("Defining eigenvalue solver");
    eig::solver solver;
    auto init = initial_mps.get_tensor();

    //    if(tensors.measurements.energy_variance)
    //        solver.config.tol = std::clamp( 1e-2 * tensors.measurements.energy_variance.value(), 1e-16 , settings::precision::eig_tolerance);
    solver.config.tol = settings::precision::eig_tolerance;
    solver.config.compress = settings::precision::use_compressed_mpo_squared_otf;
    solver.config.initial_guess.push_back({init.data(),0});
    solver.config.compute_eigvecs = eig::Vecs::ON;
    solver.config.sigma = 1.0;
    solver.config.maxNev = static_cast<eig::size_type>(1);;
    solver.config.maxNcv = static_cast<eig::size_type>(settings::precision::eig_default_ncv);;
    solver.config.ritz = ritz;
    solver.setLogLevel(2);

    // Since we use reduced energy mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
    // would otherwise cause trouble for Arpack. This equates to subtracting sigma * identity from the bottom corner of the mpo.
    // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
    tools::log->trace("Finding ground state");
    solver.eigs(matrix);

    std::vector<opt_mps> results;
    internal::krylov_extract_solutions(tensors, initial_mps, solver, results, meta, true);
    auto comp_energy = [&ritz, &meta](const opt_mps &lhs, const opt_mps &rhs) {
        switch(ritz) {
            case eig::Ritz::SR: return lhs.get_energy() <= rhs.get_energy();
            case eig::Ritz::LR: return lhs.get_energy() >= rhs.get_energy();
            default: throw std::runtime_error(fmt::format("Ground state optimization with ritz {} is not implemented", enum2sv(meta.optRitz)));
        }
    };

    if(results.size() >= 2) std::sort(results.begin(), results.end(), comp_energy);

    constexpr size_t max_print = settings::debug ? 32 : 4;
    for(const auto &[num, mps] : iter::enumerate(results)) {
        if(num >= max_print) break;
        reports::krylov_add_entry(mps);
    }
    reports::print_krylov_report();
    if(not results.empty())
        return results.front();
    else
        return initial_mps; // Solver failed

    //    try{
    //        solver.eigs(matrix, nev, ncv, ritz, eig::Form::SYMM, eig::Side::R, 1.0, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF);
    //        return tenx::asNormalized(eig::view::get_eigvec<Scalar>(solver.result, shape_mps, 0));
    //    }catch(const std::exception & ex){
    //        tools::log->warn("Arpack failed. shape_mps {} | shape_mpo {} | sites {}", shape_mps, shape_mpo, tensors.active_sites);
    //        return tenx::asNormalized(tensors.get_multisite_mps());
    //    }
}
