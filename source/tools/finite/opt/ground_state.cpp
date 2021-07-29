#include "../opt_mps.h"
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

tools::finite::opt::opt_mps tools::finite::opt::internal::ground_state_optimization(const TensorsFinite &tensors, const AlgorithmStatus &status,
                                                                                    StateRitz ritz) {
    double  energy_reduced = 0.0;
    opt_mps initial_mps("current state", tensors.get_multisite_mps(), tensors.active_sites,
                        tools::finite::measure::energy(tensors) - energy_reduced, // Eigval
                        energy_reduced,                                           // Energy reduced for full system
                        tools::finite::measure::energy_variance(tensors),
                        1.0, // Overlap
                        tensors.get_length());

    return ground_state_optimization(initial_mps, tensors, status, enum2str(ritz));
}

tools::finite::opt::opt_mps tools::finite::opt::internal::ground_state_optimization(const opt_mps &initial_mps, const TensorsFinite &tensors,
                                                                                    [[maybe_unused]] const AlgorithmStatus &status,
                                                                                    std::string_view                        ritzstring) {
    tools::log->debug("Ground state optimization with ritz {} ...", ritzstring);
    using namespace internal;
    using namespace settings::precision;
    auto t_gs = tid::tic_scope("gs");

    eig::Ritz   ritz      = eig::stringToRitz(ritzstring);
    const auto &mpo       = tensors.get_multisite_mpo();
    const auto &env       = tensors.get_multisite_ene_blk();
    auto        shape_mps = tensors.active_problem_dims();
    auto        shape_mpo = mpo.dimensions();
    //    int nev = std::min(4l,(long)(tensors.state->active_problem_size()/2));
    auto nev = static_cast<eig::size_type>(1);
    auto ncv = static_cast<eig::size_type>(settings::precision::eig_default_ncv);
    tools::log->trace("Defining Hamiltonian matrix-vector product");
    MatVecMPO<cplx> matrix(env.L.data(), env.R.data(), mpo.data(), shape_mps, shape_mpo);
    tools::log->trace("Defining eigenvalue solver");
    eig::solver solver;

    //    if(tensors.measurements.energy_variance)
    //        solver.config.tol = std::clamp( 1e-2 * tensors.measurements.energy_variance.value(), 1e-16 , settings::precision::eig_tolerance);
    solver.config.tol = settings::precision::eig_tolerance;
    solver.setLogLevel(2);

    // Since we use reduced energy mpo's, we set a sigma shift = 1.0 to move the smallest eigenvalue away from 0, which
    // would otherwise cause trouble for Arpack. This equates to subtracting sigma * identity from the bottom corner of the mpo.
    // The resulting eigenvalue will be shifted by the same amount, but the eigenvector will be the same, and that's what we keep.
    tools::log->trace("Finding ground state");
    solver.eigs(matrix, nev, ncv, ritz, eig::Form::SYMM, eig::Side::R, 1.0, eig::Shinv::OFF, eig::Vecs::ON, eig::Dephase::OFF);

    std::vector<opt_mps> eigvecs_mps;
    internal::krylov_extract_solutions(initial_mps, tensors, solver, eigvecs_mps);
    auto comp_energy = [&ritz, &ritzstring](const opt_mps &lhs, const opt_mps &rhs) {
        switch(ritz) {
            case eig::Ritz::SR: return lhs.get_energy() <= rhs.get_energy();
            case eig::Ritz::LR: return lhs.get_energy() >= rhs.get_energy();
            default: throw std::runtime_error(fmt::format("Ground state optimization with ritz {} is not implemented", ritzstring));
        }
    };

    if(eigvecs_mps.size() >= 2) std::sort(eigvecs_mps.begin(), eigvecs_mps.end(), comp_energy);

    constexpr size_t max_print = settings::debug ? 32 : 4;
    for(const auto &[num, mps] : iter::enumerate(eigvecs_mps)) {
        if(num >= max_print) break;
        internal::reports::krylov_add_entry(mps);
    }
    internal::reports::print_krylov_report();
    if(not eigvecs_mps.empty())
        return eigvecs_mps.front();
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
