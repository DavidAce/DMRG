//
// Created by david on 2020-05-11.
//

#include "class_model_finite.h"
#include <tensors/model/class_mpo_base.h>
#include <tools/finite/mpo.h>
#include <tools/finite/multisite.h>
#include <tools/finite/views.h>

// We need to make a destructor manually for the enclosing class "class_model_finite"
// that encloses "class_model_base". Otherwise unique_ptr will forcibly inline its
// own default deleter.
// This allows us to forward declare the abstract base class "class_model_base"
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_model_finite::~class_model_finite() = default;

class_model_finite::class_model_finite(const class_model_finite &other) { *this = other; }

class_model_finite &class_model_finite::operator=(const class_model_finite &other) {
    // check for self-assignment
    if(&other == this) return *this;
    // The MPO's are special and the whole point of doing this manually
    this->MPO.clear();
    for(auto &mpo : other.MPO) this->MPO.emplace_back(mpo->clone());
    this->cache = other.cache;
    return *this;
}

const class_mpo_base &class_model_finite::get_mpo(size_t pos) const {
    if(pos >= MPO.size()) throw std::range_error(fmt::format("get_mpo(pos) pos out of range: {}", pos));
    return **std::next(MPO.begin(), static_cast<long>(pos));
}

class_mpo_base &class_model_finite::get_mpo(size_t pos) { return const_cast<class_mpo_base &>(static_cast<const class_model_finite &>(*this).get_mpo(pos)); }

size_t class_model_finite::get_length() const { return MPO.size(); }

bool class_model_finite::is_real() const {
    for(auto &mpo : MPO)
        if(not mpo->isReal()) return false;
    ;
    return true;
}

bool class_model_finite::has_nan() const {
    for(auto &mpo : MPO)
        if(mpo->hasNaN()) return true;
    return false;
}

void class_model_finite::assert_validity() const {
    for(auto &mpo : MPO) mpo->assertValidity();
}

// For reduced energy MPO's

bool class_model_finite::is_reduced() const {
    bool reduced = MPO.front()->is_reduced();
    for(auto &mpo : MPO)
        if(reduced != mpo->is_reduced())
            throw std::runtime_error(
                fmt::format("First MPO has is_reduced: {}, but MPO at pos {} has is_reduced: {}", reduced, mpo->get_position(), mpo->is_reduced()));
    return reduced;
}

double class_model_finite::get_energy_reduced() const { return get_energy_per_site_reduced() * static_cast<double>(get_length()); }

double class_model_finite::get_energy_per_site_reduced() const {
    // Check that all energies are the same
    double e_reduced = MPO.front()->get_reduced_energy();
    for(auto &mpo : MPO)
        if(mpo->get_reduced_energy() != e_reduced) throw std::runtime_error("Reduced energy mismatch!");
    return e_reduced;
}

void class_model_finite::set_reduced_energy(double total_energy) { set_reduced_energy_per_site(total_energy / static_cast<double>(get_length())); }

void class_model_finite::set_reduced_energy_per_site(double site_energy) {
    if(get_energy_per_site_reduced() == site_energy) return;
    clear_cache();
    for(auto &mpo : MPO) mpo->set_reduced_energy(site_energy);
    std::cerr << "MUST REBUILD ENVIRONMENTS AFTER SETTING REDUCED ENERGY ON MPO'S!!" << std::endl;
    //    tools::finite::mps::rebuild_edges(*this);
}

void class_model_finite::perturb_hamiltonian(double coupling_ptb, double field_ptb, PerturbMode perturbMode) {
    tools::finite::mpo::perturb_hamiltonian(*this, coupling_ptb, field_ptb, perturbMode);
}

void class_model_finite::damp_hamiltonian(double coupling_damp, double field_damp) { tools::finite::mpo::damp_hamiltonian(*this, coupling_damp, field_damp); }

bool class_model_finite::is_perturbed() const {
    for(size_t pos = 0; pos < get_length(); pos++) {
        if(get_mpo(pos).is_perturbed()) return true;
    }
    return false;
}

bool class_model_finite::is_damped() const {
    for(size_t pos = 0; pos < get_length(); pos++) {
        if(get_mpo(pos).is_damped()) return true;
    }
    return false;
}

Eigen::DSizes<long, 4> class_model_finite::active_dimensions() const { return tools::finite::multisite::get_dimensions(*this); }

const Eigen::Tensor<class_model_finite::Scalar, 4> &class_model_finite::get_multisite_mpo() const { return tools::finite::views::get_multisite_mpo(*this); }

void class_model_finite::clear_cache() const { cache = Cache(); }
