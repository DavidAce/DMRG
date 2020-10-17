//
// Created by david on 2020-05-11.
//
#include <general/nmspc_tensor_extra.h>
// -- (textra first)
#include "class_model_finite.h"
#include "class_mpo_factory.h"
#include <tools/common/prof.h>
#include <tools/finite/multisite.h>
#include <tools/finite/views.h>

class_model_finite::class_model_finite() = default; // Can't initialize lists since we don't know the model size yet

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_model_finite::~class_model_finite()                                   = default;            // default dtor
class_model_finite::class_model_finite(class_model_finite &&other) = default;            // default move ctor
class_model_finite &class_model_finite::operator=(class_model_finite &&other) = default; // default move assign

/* clang-format off */
class_model_finite::class_model_finite(const class_model_finite &other) :
    cache(other.cache),
    active_sites(other.active_sites),
    model_type(other.model_type)
{
    MPO.clear();
    for(const auto &other_mpo : other.MPO) MPO.emplace_back(other_mpo->clone());
}
/* clang-format on */

class_model_finite &class_model_finite::operator=(const class_model_finite &other) {
    // check for self-assignment
    if(this != &other) {
        cache = other.cache;
        MPO.clear();
        for(auto &other_mpo : other.MPO) MPO.emplace_back(other_mpo->clone());
        active_sites = other.active_sites;
        model_type   = other.model_type;
    }
    return *this;
}

void class_model_finite::initialize(ModelType model_type_, size_t model_size) {
    tools::log->info("Initializing model with {} sites", model_size);
    if(model_size < 2) throw std::logic_error("Tried to initialize model with less than 2 sites");
    if(model_size > 2048) throw std::logic_error("Tried to initialize model with more than 2048 sites");
    if(not MPO.empty()) throw std::logic_error("Tried to initialize over an existing model. This is usually not what you want!");
    // Generate MPO
    model_type = model_type_;
    for(size_t site = 0; site < model_size; site++) {
        MPO.emplace_back(class_mpo_factory::create_mpo(site, model_type));
    }
    if(MPO.size() != model_size) throw std::logic_error("Initialized MPO with wrong size");
}

void class_model_finite::randomize(){
    tools::log->info("Randomizing model");
    std::vector<class_mpo_site::TableMap> all_params;
    for(auto &mpo : MPO) {
        mpo->randomize_hamiltonian();
        all_params.emplace_back(mpo->get_parameters());
    }
    for(auto &mpo : MPO) mpo->set_averages(all_params, false);
}


const class_mpo_site &class_model_finite::get_mpo(size_t pos) const {
    if(pos >= MPO.size()) throw std::range_error(fmt::format("get_2site_tensor(pos) pos out of range: {}", pos));
    return **std::next(MPO.begin(), static_cast<long>(pos));
}

class_mpo_site &class_model_finite::get_mpo(size_t pos) { return const_cast<class_mpo_site &>(static_cast<const class_model_finite &>(*this).get_mpo(pos)); }

size_t class_model_finite::get_length() const { return MPO.size(); }

bool class_model_finite::is_real() const {
    for(auto &mpo : MPO)
        if(not mpo->is_real()) return false;
    ;
    return true;
}

bool class_model_finite::has_nan() const {
    for(auto &mpo : MPO)
        if(mpo->has_nan()) return true;
    return false;
}

void class_model_finite::assert_validity() const {
    for(auto &mpo : MPO) mpo->assert_validity();
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

}

void class_model_finite::perturb_hamiltonian(double coupling_ptb, double field_ptb, PerturbMode perturbMode) {
    std::vector<class_mpo_site::TableMap> all_params;
    for(auto &mpo : MPO) {
        mpo->set_perturbation(coupling_ptb, field_ptb, perturbMode);
        all_params.push_back(mpo->get_parameters());
    }
    for(auto &mpo : MPO) mpo->set_averages(all_params, false);
    clear_cache();
    if(coupling_ptb == 0.0 and field_ptb == 0.0 and is_perturbed()) throw std::runtime_error("Model: Should have unperturbed!");
}

void class_model_finite::damp_hamiltonian(double coupling_damp, double field_damp) {
    for(auto &mpo : MPO) {
        mpo->set_coupling_damping(coupling_damp);
        mpo->set_field_damping(field_damp);
    }
    clear_cache();
    if(coupling_damp == 0.0 and field_damp == 0.0 and is_damped()) throw std::runtime_error("Model: Should have undamped!");
}

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

Eigen::Tensor<class_model_finite::Scalar, 4> class_model_finite::get_multisite_mpo(const std::vector<size_t> &sites) const {
    if(sites.empty()) throw std::runtime_error("No active sites on which to build a multisite mpo tensor");
    if(sites == active_sites and cache.multisite_mpo) return cache.multisite_mpo.value();
    tools::log->trace("Contracting multisite mpo tensor with {} sites", sites.size());
    tools::common::profile::get_default_prof()["t_mpo"]->tic();
    Eigen::Tensor<Scalar, 4> multisite_mpo;
    constexpr auto           shuffle_idx  = Textra::array6{0, 3, 1, 4, 2, 5};
    constexpr auto           contract_idx = Textra::idx({1}, {0});
    Textra::array4           new_dims;
    Eigen::Tensor<Scalar, 4> temp;
    bool                     first = true;
    for(auto &site : sites) {
        if(first) {
            multisite_mpo    = get_mpo(site).MPO();
            first            = false;
            continue;
        }
        const auto &M    = get_mpo(site).MPO();
        long        dim0 = multisite_mpo.dimension(0);
        long        dim1 = M.dimension(1);
        long        dim2 = multisite_mpo.dimension(2) * M.dimension(2);
        long        dim3 = multisite_mpo.dimension(3) * M.dimension(3);
        new_dims         = {dim0, dim1, dim2, dim3};
        temp.resize(new_dims);
        temp.device(*Textra::omp::dev) = multisite_mpo.contract(M, contract_idx).shuffle(shuffle_idx).reshape(new_dims);
        multisite_mpo                  = temp;
    }
    tools::common::profile::get_default_prof()["t_mpo"]->toc();
    return multisite_mpo;
}

const Eigen::Tensor<class_model_finite::Scalar, 4> &class_model_finite::get_multisite_tensor() const {
    if(cache.multisite_tensor) return cache.multisite_tensor.value();
    cache.multisite_tensor = get_multisite_tensor(active_sites);
    return cache.multisite_tensor.value();
}

void class_model_finite::clear_cache() const {
    tools::log->trace("Clearing model cache");
    cache = Cache();
}
