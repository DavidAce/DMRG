
#include "class_tensors_finite.h"
#include <math/nmspc_math.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/edges/class_env_ene.h>
#include <tensors/edges/class_env_var.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/finite/measure.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>

class_tensors_finite::class_tensors_finite() {
    state = std::make_unique<class_state_finite>();
    model = std::make_unique<class_model_finite>();
    edges = std::make_unique<class_edges_finite>();
}

void class_tensors_finite::initialize(ModelType model_type, size_t model_size, size_t position) {
    state->set_chi_lim(2); // Can't call chi_init() <-- it's a pure virtual function
    if(state->has_nan()) throw std::runtime_error("State has NAN's before initializing it");

    tools::finite::mps::initialize(*state, model_type, model_size, position);
    tools::finite::mpo::initialize(*model, model_type, model_size);
    tools::finite::mpo::randomize(*model);
}

// Active sites
void class_tensors_finite::sync_active_sites() {
    active_sites = state->active_sites;
    model->active_sites = state->active_sites;
    edges->active_sites = state->active_sites;
}

void class_tensors_finite::activate_sites(const size_t threshold, const size_t max_sites, const size_t min_sites){
    active_sites = state->activate_sites(threshold, max_sites, min_sites);
    model->active_sites = state->active_sites;
    edges->active_sites = state->active_sites;
}

bool class_tensors_finite::is_real() const { return state->is_real() and model->is_real() and edges->is_real(); }
bool class_tensors_finite::has_nan() const { return state->has_nan() or model->has_nan() or edges->has_nan(); }
void class_tensors_finite::assert_validity() const {
    state->assert_validity();
    model->assert_validity();
    edges->assert_validity();
}

size_t class_tensors_finite::get_length() const {
    if(not math::all_equal(state->get_length(), model->get_length(), edges->get_length()))
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


void class_tensors_finite::move_center_point() const {


}


void class_tensors_finite::do_all_measurements() const {
    tools::finite::measure::do_all_measurements(*this);
}


void class_tensors_finite::clear_measurements() const {
    measurements = tensors_measure_finite();
    state->clear_measurements();
}
void class_tensors_finite::clear_cache() const {
    state->clear_cache();
    model->clear_cache();
}