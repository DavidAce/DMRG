
#include "class_tensors_finite.h"
#include <math/num.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/common/log.h>
#include <tools/finite/env.h>
#include <tools/finite/measure.h>
#include <tools/finite/mps.h>
#include <tools/finite/multisite.h>
#include <tools/finite/ops.h>

class_tensors_finite::class_tensors_finite()
    : state(std::make_unique<class_state_finite>()), model(std::make_unique<class_model_finite>()), edges(std::make_unique<class_edges_finite>()) {
    tools::log->trace("Constructing class_tensors_finite");
}

// We need to define the destructor and other special functions
// because we enclose data in unique_ptr for this pimpl idiom.
// Otherwise unique_ptr will forcibly inline its own default deleter.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_tensors_finite::~class_tensors_finite()                                        = default; // default dtor
class_tensors_finite::class_tensors_finite(class_tensors_finite &&other)             = default; // default move ctor
class_tensors_finite &class_tensors_finite::operator=(class_tensors_finite &&other)  = default; // default move assign

class_tensors_finite::class_tensors_finite(const class_tensors_finite &other)
    : state(std::make_unique<class_state_finite>(*other.state)), model(std::make_unique<class_model_finite>(*other.model)),
      edges(std::make_unique<class_edges_finite>(*other.edges)), active_sites(other.active_sites), measurements(other.measurements) {}
class_tensors_finite &class_tensors_finite::operator=(const class_tensors_finite &other) {
    // check for self-assignment
    if(this != &other) {
        state        = std::make_unique<class_state_finite>(*other.state);
        model        = std::make_unique<class_model_finite>(*other.model);
        edges        = std::make_unique<class_edges_finite>(*other.edges);
        active_sites = other.active_sites;
        measurements = other.measurements;
    }
    return *this;
}

void class_tensors_finite::initialize(ModelType model_type, size_t model_size, size_t position) {
    state->initialize(model_type, model_size, position);
    model->initialize(model_type, model_size);
    edges->initialize(model_size);
}

void class_tensors_finite::randomize_model() {
    eject_all_edges();
    model->randomize();
    rebuild_edges();
}

void class_tensors_finite::randomize_state(StateType state_type,const std::string &sector, long chi_lim, bool use_eigenspinors, std::optional<long> bitfield, std::optional<double> svd_threshold) {
    clear_measurements();
    tools::finite::mps::randomize_state(*state, sector,state_type, chi_lim, use_eigenspinors, bitfield);
    if(state_type == StateType::RANDOMIZE_PREVIOUS_STATE or state_type == StateType::RANDOM_ENTANGLED_STATE)
        project_to_nearest_sector(sector,chi_lim,svd_threshold);
    else
        normalize_state(chi_lim,svd_threshold,NormPolicy::ALWAYS);
}

void class_tensors_finite::normalize_state(long chi_lim, std::optional<double> svd_threshold, NormPolicy norm_policy) {
    // Normalize if unity was lost for some reason (numerical error buildup)
    auto has_normalized = tools::finite::mps::normalize_state(*state, chi_lim, svd_threshold, norm_policy);
    if(has_normalized) {
        clear_measurements();
        eject_all_edges();
        rebuild_edges();
        assert_validity();
    }
}


void class_tensors_finite::project_to_nearest_sector(const std::string &sector, std::optional<long> chi_lim, std::optional<double> svd_threshold) {
    clear_measurements();
    if (not chi_lim) chi_lim =  state->find_largest_chi();
    tools::finite::mps::normalize_state(*state, chi_lim.value(), svd_threshold, NormPolicy::IFNEEDED);
    tools::finite::ops::project_to_nearest_sector(*state, sector);
    normalize_state(chi_lim.value(),svd_threshold, NormPolicy::ALWAYS);
}

void class_tensors_finite::perturb_model_params(double coupling_ptb, double field_ptb, PerturbMode perturbMode) {
    measurements = tensors_measure_finite(); // State measurements can remain
    model->perturb_hamiltonian(coupling_ptb, field_ptb, perturbMode);
    eject_all_edges();
    rebuild_edges();
    model->assert_validity();
}

void class_tensors_finite::reduce_mpo_energy(std::optional<double> site_energy) {
    if(not site_energy)
        site_energy = tools::finite::measure::energy_per_site(*this);
    measurements = tensors_measure_finite(); // State measurements can remain
    model->set_reduced_energy_per_site(site_energy.value());
    clear_cache();
    clear_measurements();
    eject_all_edges();
    rebuild_edges();
    model->assert_validity();

}

void class_tensors_finite::damp_model_disorder(double coupling_damp, double field_damp) {
    measurements = tensors_measure_finite(); // State measurements can remain
    model->damp_hamiltonian(coupling_damp, field_damp);
    eject_all_edges();
    rebuild_edges();
    model->assert_validity();
}

// Active sites
void class_tensors_finite::sync_active_sites() {
    active_sites        = state->active_sites;
    model->active_sites = state->active_sites;
    edges->active_sites = state->active_sites;
    rebuild_edges();
}

void class_tensors_finite::activate_sites(long threshold, size_t max_sites, size_t min_sites) {
    active_sites        = tools::finite::multisite::generate_site_list(*state, threshold, max_sites, min_sites);
    state->active_sites = active_sites;
    model->active_sites = active_sites;
    edges->active_sites = active_sites;
    rebuild_edges();
    clear_cache();
    clear_measurements();
}

void class_tensors_finite::activate_truncated_sites(long threshold, long chi_lim, size_t max_sites, size_t min_sites) {
    active_sites        = tools::finite::multisite::generate_truncated_site_list(*state, threshold, chi_lim, max_sites, min_sites);
    state->active_sites = active_sites;
    model->active_sites = active_sites;
    edges->active_sites = active_sites;
    rebuild_edges(); // will only produce missing edges
    clear_cache();   // clear multisite tensors and mpos which become invalid
    // TODO: Should measurements be cleared here? Probably not
    //    clear_measurements();
}

long class_tensors_finite::active_problem_size() const { return tools::finite::multisite::get_problem_size(*state, active_sites); }
bool class_tensors_finite::is_real() const { return state->is_real() and model->is_real() and edges->is_real(); }
bool class_tensors_finite::has_nan() const { return state->has_nan() or model->has_nan() or edges->has_nan(); }
void class_tensors_finite::assert_validity() const {
    state->assert_validity();
    model->assert_validity();
    edges->assert_validity();
}

size_t class_tensors_finite::get_length() const {
    if(not num::all_equal(state->get_length(), model->get_length(), edges->get_length()))
        throw std::runtime_error("All lengths are not equal: state {} | model {} | edges {}");
    return state->get_length();
}

size_t class_tensors_finite::get_position() const { return state->get_position(); }

bool class_tensors_finite::position_is_the_middle() const { return state->position_is_the_middle(); }
bool class_tensors_finite::position_is_the_middle_any_direction() const { return state->position_is_the_middle_any_direction(); }
bool class_tensors_finite::position_is_left_edge() const { return state->position_is_left_edge(); }
bool class_tensors_finite::position_is_right_edge() const { return state->position_is_right_edge(); }
bool class_tensors_finite::position_is_any_edge() const { return state->position_is_any_edge(); }
bool class_tensors_finite::position_is_at(size_t pos) const { return state->position_is_at(pos); }
void class_tensors_finite::move_center_point(long chi_lim,std::optional<double> svd_threshold) { tools::finite::mps::move_center_point(*state, chi_lim, svd_threshold); }
void class_tensors_finite::merge_multisite_tensor(const Eigen::Tensor<Scalar, 3> &multisite_tensor, long chi_lim, std::optional<double> svd_threshold) {
    // Make sure the active sites are the same everywhere
    if(not num::all_equal(active_sites, state->active_sites, model->active_sites, edges->active_sites))
        throw std::runtime_error("All active sites are not equal: tensors {} | state {} | model {} | edges {}");
    clear_measurements();
    tools::finite::mps::merge_multisite_tensor(*state, multisite_tensor, active_sites, get_position(), chi_lim);
    normalize_state(chi_lim, svd_threshold,NormPolicy::IFNEEDED);
}

void class_tensors_finite::rebuild_edges() { tools::finite::env::rebuild_edges(*state, *model, *edges); }
void class_tensors_finite::eject_all_edges() { edges->eject_all_edges(); }
void class_tensors_finite::eject_inactive_edges() { edges->eject_inactive_edges(); }
void class_tensors_finite::do_all_measurements() const { tools::finite::measure::do_all_measurements(*this); }
void class_tensors_finite::clear_measurements() const {
    tools::log->trace("Clearing tensor measurements");
    measurements = tensors_measure_finite();
    state->clear_measurements();
}

void class_tensors_finite::clear_cache() const {
    model->clear_cache();
    state->clear_cache();
}