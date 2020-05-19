//
// Created by david on 2019-06-25.
//
#include <math/nmspc_math.h>
#include <math/nmspc_random.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/model/class_mpo_factory.h>
#include <tensors/model/class_mpo_site.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/debug.h>
#include <tools/finite/measure.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>

void tools::finite::mpo::randomize(class_model_finite &model) {
    tools::log->trace("Setting random fields in MPO's");
    std::vector<class_mpo_site::TableMap> all_params;
    for(auto &mpo : model.MPO) {
        mpo->randomize_hamiltonian();
        all_params.emplace_back(mpo->get_parameters());
    }
    for(auto &mpo : model.MPO) mpo->set_averages(all_params, false);

//    tools::finite::mps::rebuild_edges(model);
    std::cerr << "MUST REBUILD ENVIRONMENTS AFTER RANDOMIZING HAMILTONIAN!!" << std::endl;
    tools::finite::debug::check_integrity(model);
}

void tools::finite::mpo::perturb_hamiltonian(class_model_finite &model, double coupling_ptb, double field_ptb, PerturbMode perturbMode) {
    std::vector<class_mpo_site::TableMap> all_params;
    for(auto &mpo : model.MPO) {
        mpo->set_perturbation(coupling_ptb, field_ptb, perturbMode);
        all_params.push_back(mpo->get_parameters());
    }
    for(auto &mpo : model.MPO) mpo->set_averages(all_params, false);
    model.clear_cache();
    std::cerr << "MUST REBUILD ENVIRONMENTS AFTER PERTURBING HAMILTONIAN!!" << std::endl;
//    tools::finite::mps::rebuild_edges(state);
    if(coupling_ptb == 0.0 and field_ptb == 0.0 and model.is_perturbed()) throw std::runtime_error("Model: Should have unperturbed!");
}

void tools::finite::mpo::damp_hamiltonian(class_model_finite &model, double coupling_damp, double field_damp) {
    for(auto &mpo : model.MPO) {
        mpo->set_coupling_damping(coupling_damp);
        mpo->set_field_damping(field_damp);
    }
    model.clear_cache();
    std::cerr << "MUST REBUILD ENVIRONMENTS AFTER DAMPING HAMILTONIAN!!" << std::endl;
//    tools::finite::mps::rebuild_edges(state);
    if(coupling_damp == 0.0 and field_damp == 0.0 and model.is_damped()) throw std::runtime_error("State: Should have undamped!");
}

void tools::finite::mpo::reduce_mpo_energy(class_model_finite & model, double site_energy) {
//    if(model.active_sites.empty() and state.active_sites.empty()){}
//    if(not math::all_equal(state.active_sites, model.active_sites, edges.active_sites)) throw std::runtime_error(fmt::format("Active sites not equal: state {} | model {} | edges {}",state.active_sites, model.active_sites, edges.active_sites));
//    const auto mps                              = state.get_multisite_tensor();
//    double     energy_per_site_before           = tools::finite::measure::energy_per_site(mps,model,edges);
    double     energy_per_site_reduced_before   = model.get_energy_per_site_reduced();
//    double energy_per_site_minus_reduced_before = tools::finite::measure::energy_minus_energy_reduced(mps, model, edges) / static_cast<double>(state.get_length());
//    double energy_variance_per_site_before      = tools::finite::measure::energy_variance_per_site(mps,model,edges);
//    log->debug("Variance check before reduce          : {:.16f}", std::log10(tools::finite::measure::energy_variance_per_site(mps,model,edges)));
//    log->debug("Status before reduce (multi)          : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
//               energy_per_site_before, energy_per_site_reduced_before, energy_per_site_minus_reduced_before, std::log10(energy_variance_per_site_before));

//    tools::log->trace("Reducing MPO energy by: {:<20.16f}", energy_per_site_before);
    model.set_reduced_energy_per_site(site_energy);
//    double energy_per_site_after               = tools::finite::measure::energy_per_site(mps,model,edges);
    double energy_per_site_reduced_after       = model.get_energy_per_site_reduced();
//    double energy_per_site_minus_reduced_after = tools::finite::measure::energy_minus_energy_reduced(mps,model,edges) / static_cast<double>(state.get_length());
//    double energy_variance_per_site_after      = tools::finite::measure::energy_variance_per_site(mps,model,edges);
//    log->debug("Variance check after reduce           : {:.16f}", std::log10(measure::energy_variance_per_site(mps,model,edges)));
//    log->debug("Status after reduce (multi)           : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
//               energy_per_site_after, energy_per_site_reduced_after, energy_per_site_minus_reduced_after, std::log10(energy_variance_per_site_after));

//    if(std::abs(energy_per_site_before - energy_per_site_after) > 0.1) throw std::logic_error("Energy before and after mpo reduction differ");
//    if(std::abs(energy_per_site_minus_reduced_after) > 0.1) throw std::logic_error("Energy reduction failed");
    std::cerr << "MUST REBUILD ENVIRONMENTS AFTER REDUCING MPO ENERGY!!" << std::endl;
}

