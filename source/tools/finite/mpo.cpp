//
// Created by david on 2019-06-25.
//
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>
#include <tools/finite/measure.h>
#include <tools/finite/debug.h>
#include <tools/common/log.h>
#include <state/class_state_finite.h>
#include <model/class_model_base.h>
#include <model/class_model_factory.h>
#include <math/nmspc_random.h>

void tools::finite::mpo::initialize(class_state_finite & state, const std::string &model_type, const size_t num_sites, size_t position){
    tools::log->trace("Initializing mpo with {} sites at position {}", num_sites,position);
    if (num_sites < 2    )      throw std::logic_error("Tried to initialize MPO with less than 2 sites");
    if (num_sites > 2048 )      throw std::logic_error("Tried to initialize MPO with more than 2048 sites");
    if (position >= num_sites ) throw std::logic_error("Tried to initialize MPO at a position larger than the number of sites");
    state.MPO_L.clear();
    state.MPO_R.clear();
    //Generate MPO
    for(size_t site = 0; site < num_sites; site++){
        if(site <= position)
            state.MPO_L.emplace_back(class_model_factory::create_mpo(site,model_type));
        else
            state.MPO_R.emplace_back(class_model_factory::create_mpo(site,model_type));
    }
    if(state.MPO_L.size() + state.MPO_R.size() != num_sites) throw std::logic_error("Initialized MPO with wrong size");
    if(state.MPO_L.size()-1 != position) throw std::logic_error("Initialized MPO at the wrong position");
    if(state.MPO_L.back()->get_position() != position) throw std::logic_error("Initialized MPO at the wrong position");

}


void tools::finite::mpo::randomize(class_state_finite &state) {
    tools::log->trace("Setting random fields in MPO's");
    std::vector<class_model_base::TableMap> all_params;

    for (auto &mpo : state.MPO_L){
        mpo->randomize_hamiltonian();
        all_params.emplace_back(mpo->get_parameters());
    }
    for (auto &mpo : state.MPO_R){
        mpo->randomize_hamiltonian();
        all_params.emplace_back(mpo->get_parameters());
    }

    for (auto &mpo : state.MPO_L){
        mpo->set_averages(all_params, false);
    }
    for (auto &mpo : state.MPO_R){
        mpo->set_averages(all_params, false);
    }

    for (size_t pos = 0; pos < state.get_length(); pos++){
        if(state.get_MPO(pos).hasNaN()) throw std::runtime_error("MPO " + std::to_string(pos) + " has initialized with NAN's");
    }

    tools::finite::mps::rebuild_environments(state);
    tools::finite::debug::check_integrity(state);
}


void tools::finite::mpo::perturb_hamiltonian(class_state_finite &state, double coupling_ptb, double field_ptb, PerturbMode perturbMode){
    std::vector<class_model_base::TableMap> all_params;
    for (auto &mpo : state.MPO_L){
        mpo->set_perturbation(coupling_ptb, field_ptb, perturbMode);
        all_params.push_back(mpo->get_parameters());
    }
    for (auto &mpo : state.MPO_R){
        mpo->set_perturbation(coupling_ptb, field_ptb, perturbMode);
        all_params.push_back(mpo->get_parameters());
    }
    for (auto &mpo : state.MPO_L){
        mpo->set_averages(all_params, false);
    }
    for (auto &mpo : state.MPO_R){
        mpo->set_averages(all_params, false);
    }

    state.clear_cache();
    state.clear_measurements();
    tools::finite::mps::rebuild_environments(state);
    if (coupling_ptb == 0.0 and field_ptb == 0.0 and state.is_perturbed())
        throw std::runtime_error("State: Should have unperturbed!");
}

void tools::finite::mpo::damp_hamiltonian(class_state_finite &state, double coupling_damp, double field_damp){
    for (auto &mpo : state.MPO_L){
        mpo->set_coupling_damping(coupling_damp);
        mpo->set_field_damping(field_damp);
    }
    for (auto &mpo : state.MPO_R){
        mpo->set_coupling_damping(coupling_damp);
        mpo->set_field_damping(field_damp);
    }

    state.clear_cache();
    state.clear_measurements();
    tools::finite::mps::rebuild_environments(state);
    if (coupling_damp == 0.0 and field_damp == 0.0 and state.is_damped())
        throw std::runtime_error("State: Should have undamped!");
}


void tools::finite::mpo::reduce_mpo_energy(class_state_finite &state){
    state.clear_measurements();
    state.clear_cache();
    if (state.active_sites.empty())
        reduce_mpo_energy_2site(state);
    else
        reduce_mpo_energy_multi(state);
//    reduce_mpo_energy(state, state.get_multitheta());
//
}

void tools::finite::mpo::reduce_mpo_energy_multi(class_state_finite &state){
    const auto theta                                    = state.get_multitheta();
    double energy_per_site_before                       = tools::finite::measure::energy_per_site(state,theta);
    double energy_per_site_reduced_before               = state.get_energy_per_site_reduced();
    double energy_per_site_minus_reduced_before         = tools::finite::measure::energy_minus_energy_reduced(state,theta)/state.get_length();
    double energy_variance_per_site_before              = tools::finite::measure::energy_variance_per_site(state,theta);
    log->debug("Variance check before reduce          : {:.16f}", std::log10(measure::energy_variance_per_site(state,theta)));
    log->debug("Status before reduce (multi)          : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
               energy_per_site_before,
               energy_per_site_reduced_before,
               energy_per_site_minus_reduced_before,
               std::log10(energy_variance_per_site_before));

    tools::log->trace("Reducing MPO energy by: {:<20.16f}",energy_per_site_before);
    state.set_reduced_energy_per_site(energy_per_site_before);
    double energy_per_site_after                        = tools::finite::measure::energy_per_site(state,theta);
    double energy_per_site_reduced_after                = state.get_energy_per_site_reduced();
    double energy_per_site_minus_reduced_after          = tools::finite::measure::energy_minus_energy_reduced(state,theta)/state.get_length();
    double energy_variance_per_site_after               = tools::finite::measure::energy_variance_per_site(state,theta);
    log->debug("Variance check after reduce           : {:.16f}", std::log10(measure::energy_variance_per_site(state,theta)));
    log->debug("Status after reduce (multi)           : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
               energy_per_site_after,
               energy_per_site_reduced_after,
               energy_per_site_minus_reduced_after,
               std::log10(energy_variance_per_site_after));

    if (std::abs(energy_per_site_before - energy_per_site_after) > 0.1  )
        throw std::logic_error("Energy before and after mpo reduction differ");
    if (std::abs(energy_per_site_minus_reduced_after) > 0.1  )
        throw std::logic_error("Energy reduction failed");


}


void tools::finite::mpo::reduce_mpo_energy_2site(class_state_finite &state){
    state.clear_measurements();
    state.clear_cache();
    auto   theta = state.get_theta();
    double energy_per_site_before                       = tools::finite::measure::energy_per_site(state,theta);
    double energy_per_site_reduced_before               = state.get_energy_per_site_reduced();
    double energy_per_site_minus_reduced_before         = tools::finite::measure::energy_minus_energy_reduced(state,theta)/state.get_length();
    double energy_variance_per_site_before              = tools::finite::measure::energy_variance_per_site(state,theta);
//    log->debug("Variance check before reduce          : {:.16f}", std::log10(measure::energy_variance_per_site(state)));
    log->debug("Status before reduce (2site)          : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
               energy_per_site_before,
               energy_per_site_reduced_before,
               energy_per_site_minus_reduced_before,
               std::log10(energy_variance_per_site_before));


    tools::log->trace("Reducing MPO energy by: {}",energy_per_site_before);
    state.set_reduced_energy_per_site(energy_per_site_before);
    state.clear_measurements();
    state.clear_cache();
    theta = state.get_theta();
    double energy_per_site_after                        = tools::finite::measure::energy_per_site(state,theta);
    double energy_per_site_reduced_after                = state.get_energy_per_site_reduced();
    double energy_per_site_minus_reduced_after          = tools::finite::measure::energy_minus_energy_reduced(state,theta)/state.get_length();
    double energy_variance_per_site_after               = tools::finite::measure::energy_variance_per_site(state,theta);
    //    log->debug("Variance check after reduce          : {:.16f}", std::log10(measure::energy_variance_per_site(state)));
    log->debug("Status after reduce (2site)           : E = {:<20.16f} | E_red = {:<20.16f} | E - E_red = {:<20.16f} | Var E = {:<20.16f}",
               energy_per_site_after,
               energy_per_site_reduced_after,
               energy_per_site_minus_reduced_after,
               std::log10(energy_variance_per_site_after));
    if (std::abs(energy_per_site_before - energy_per_site_after) > 0.1  )
        throw std::logic_error("Energy before and after differ");
    if (std::abs(energy_per_site_minus_reduced_after) > 0.1  )
        throw std::logic_error("Energy reduction failed");
}