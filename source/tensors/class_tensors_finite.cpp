
#include "class_tensors_finite.h"
#include <math/nmspc_math.h>
#include <tensors/edges/class_edges_finite.h>
#include <tensors/model/class_model_finite.h>
#include <tensors/state/class_state_finite.h>
#include <tools/finite/env.h>
#include <tools/finite/measure.h>
#include <tools/finite/mpo.h>
#include <tools/finite/mps.h>

class_tensors_finite::class_tensors_finite() :
    state(std::make_unique<class_state_finite>()),
    model(std::make_unique<class_model_finite>()),
    edges(std::make_unique<class_edges_finite>())
{
    tools::log->trace("Constructing tensors");
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
class_tensors_finite::~class_tensors_finite() = default;                                                    // default dtor
class_tensors_finite::class_tensors_finite(class_tensors_finite &&other)  noexcept = default;               // default move ctor
class_tensors_finite &class_tensors_finite::operator=(class_tensors_finite &&other) noexcept = default;     // default move assign

class_tensors_finite::class_tensors_finite(const class_tensors_finite &other):
    state(std::make_unique<class_state_finite>(*other.state)),
    model(std::make_unique<class_model_finite>(*other.model)),
    edges(std::make_unique<class_edges_finite>(*other.edges)),
    active_sites(other.active_sites),
    measurements(other.measurements)
{

}
class_tensors_finite &class_tensors_finite::operator=(const class_tensors_finite &other)
{
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
    state->initialize(model_type,model_size, position);
    model->initialize(model_type,model_size);
    edges->initialize(model_size);
    tools::finite::env::rebuild_edges(*state,*model,*edges);
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


void class_tensors_finite::move_center_point() {


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