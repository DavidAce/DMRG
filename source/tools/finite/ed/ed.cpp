#include "../ed.h"
#include "algorithms/AlgorithmStatus.h"
#include "general/iter.h"
#include "math/eig.h"
#include "math/num.h"
#include "math/tenx.h"
#include "tensors/model/ModelFinite.h"
#include "tensors/state/StateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tid/tid.h"
#include "tools/common/split.h"
#include "tools/finite/measure.h"
#include "tools/finite/opt/opt-internal.h"
#include "tools/finite/opt/report.h"
#include "tools/finite/opt_meta.h"
#include "tools/finite/opt_mps.h"
namespace tools::finite::ed {
    StateFinite find_exact_state(const TensorsFinite &tensors, const AlgorithmStatus &status) {
        auto sites      = num::range<size_t>(0ul, tensors.get_length<size_t>());
        auto tensors_ed = tensors;
        tensors_ed.clear_cache();
        tensors_ed.clear_measurements();
        tensors_ed.shift_mpo_energy(0);
        // The reduction clears our squared mpo's. So we have to rebuild.
        tensors_ed.rebuild_mpo_squared();
        tensors_ed.rebuild_edges();
        tensors_ed.activate_sites(sites);
        auto t_ham       = tid::tic_scope("get_ham");
        auto hamiltonian = tensors_ed.model->get_multisite_ham();
        t_ham.toc();
        auto t_tgt = tid::tic_scope("tgt_mps");

        tools::finite::opt::opt_mps target_mps("target_mps", tensors_ed.state->get_multisite_mps(), sites,
                                               //                        tools::finite::measure::energy(tensors_ed) - energy_shift, // Eigval
                                               tools::finite::measure::energy_minus_energy_shift(tensors_ed), // Eigval
                                               tools::finite::measure::energy_shift(tensors_ed),              // Energy shift for full system
                                               tools::finite::measure::energy_variance(tensors_ed),
                                               1.0, // Overlap
                                               tensors_ed.get_length());
        t_tgt.toc();
        tools::finite::opt::OptMeta meta;
        meta.optType      = OptType::CPLX;
        meta.problem_dims = target_mps.get_tensor().dimensions();
        meta.problem_size = target_mps.get_tensor().size();

        auto solver = eig::solver();
        solver.eig(hamiltonian.data(), hamiltonian.dimension(0));
        std::vector<tools::finite::opt::opt_mps> results;
        auto                                     t_ext = tid::tic_scope("extract");

        tools::finite::opt::internal::eigs_extract_results(tensors_ed, target_mps, meta, solver, results, true, 1.5);
        t_ext.toc();
        for(const auto &[num, mps] : iter::enumerate(results)) { tools::finite::opt::reports::eigs_add_entry(mps); }
        tools::finite::opt::reports::print_eigs_report();

        auto comparator = [](const tools::finite::opt::opt_mps &lhs, const tools::finite::opt::opt_mps &rhs) {
            return lhs.get_overlap() > rhs.get_overlap();
        };
        if(results.size() >= 2) std::sort(results.begin(), results.end(), comparator);

        auto  spin_dims     = std::vector<long>(tensors_ed.get_length<size_t>(), 2l);
        auto  positions     = num::range<size_t>(0, tensors_ed.get_length<size_t>());
        auto  centerpos     = tensors_ed.get_position<long>();
        auto &state_ed      = *tensors_ed.state;
        auto &multisite_mps = results.front().get_tensor();
        auto  mps_list      = tools::common::split::split_mps(multisite_mps, spin_dims, positions, centerpos, std::nullopt); // No svd truncation
        state_ed.set_name("state_ed");
        state_ed.set_mps_sites(mps_list);
        tools::finite::measure::do_all_measurements(state_ed);

        tools::log->info("Bond dimensions χ                  = {}", tools::finite::measure::bond_dimensions(state_ed));
        tools::log->info("Bond dimension  χ (mid)            = {}", tools::finite::measure::bond_dimension_midchain(state_ed));
        tools::log->info("Entanglement entropies Sₑ          = {:8.2e}", fmt::join(tools::finite::measure::entanglement_entropies(state_ed), ", "));
        tools::log->info("Entanglement entropy   Sₑ (mid)    = {:8.2e}", tools::finite::measure::entanglement_entropy_midchain(state_ed), ", ");
        if(status.algo_type == AlgorithmType::fLBIT) {
            tools::log->info("Number entropies Sₙ                = {:8.2e}", fmt::join(tools::finite::measure::number_entropies(state_ed), ", "));
            tools::log->info("Number entropy   Sₙ (mid)          = {:8.2e}", tools::finite::measure::number_entropy_midchain(state_ed), ", ");
        }
        tools::log->info("Spin components (global X,Y,Z)     = {:8.2e}", fmt::join(tools::finite::measure::spin_components(state_ed), ", "));

        auto expectation_values_xyz = tools::finite::measure::expectation_values_xyz(state_ed);
        auto structure_factor_xyz   = tools::finite::measure::structure_factor_xyz(state_ed);
        auto kvornings_marker       = tools::finite::measure::kvornings_marker(state_ed);

        tools::log->info("Expectation values ⟨σx⟩            = {:+9.6f}", fmt::join(tenx::span(expectation_values_xyz[0]), ", "));
        tools::log->info("Expectation values ⟨σy⟩            = {:+9.6f}", fmt::join(tenx::span(expectation_values_xyz[1]), ", "));
        tools::log->info("Expectation values ⟨σz⟩            = {:+9.6f}", fmt::join(tenx::span(expectation_values_xyz[2]), ", "));
        tools::log->info("Structure f. L⁻¹ ∑_ij ⟨σx_i σx_j⟩² = {:+.16f}", structure_factor_xyz[0]);
        tools::log->info("Structure f. L⁻¹ ∑_ij ⟨σy_i σy_j⟩² = {:+.16f}", structure_factor_xyz[1]);
        tools::log->info("Structure f. L⁻¹ ∑_ij ⟨σz_i σz_j⟩² = {:+.16f}", structure_factor_xyz[2]);
        tools::log->info("Kvornings marker                   = {:+9.6f}", fmt::join(tenx::span(kvornings_marker), ", "));
        tools::log->info("Truncation Errors ε                = {:8.2e}", fmt::join(state_ed.get_truncation_errors(), ", "));

        return state_ed;
    }

}
