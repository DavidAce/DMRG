//
// Created by david on 2020-05-14.
//

#include "class_tensors_infinite.h"
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/state/class_state_infinite.h>
#include <tensors/state/class_mps_site.h>
#include <tools/infinite/env.h>
#include <tools/infinite/measure.h>
#include <tools/common/log.h>

class_tensors_infinite::class_tensors_infinite():
    state(std::make_unique<class_state_infinite>()),
    model(std::make_unique<class_model_infinite>()),
    edges(std::make_unique<class_edges_infinite>())
{
    tools::log->trace("Constructing tensors");
}


// We need to make a destructor manually because we enclose
// the classes with unique_ptr. Otherwise unique_ptr will
// forcibly inline its own default deleter.
// This is a classic pimpl idiom.
// Here we follow "rule of five", so we must also define
// our own copy/move ctor and copy/move assignments
// This has the side effect that we must define our own
// operator= and copy assignment constructor.
// Read more: https://stackoverflow.com/questions/33212686/how-to-use-unique-ptr-with-forward-declared-type
// And here:  https://stackoverflow.com/questions/6012157/is-stdunique-ptrt-required-to-know-the-full-definition-of-t
class_tensors_infinite::~class_tensors_infinite() = default;                                                      // default dtor
class_tensors_infinite::class_tensors_infinite(class_tensors_infinite &&other)  noexcept = default;               // default move ctor
class_tensors_infinite &class_tensors_infinite::operator=(class_tensors_infinite &&other) noexcept = default;     // default move assign

class_tensors_infinite::class_tensors_infinite(const class_tensors_infinite &other):
    state(std::make_unique<class_state_infinite>(*other.state)),
    model(std::make_unique<class_model_infinite>(*other.model)),
    edges(std::make_unique<class_edges_infinite>(*other.edges)),
    measurements(other.measurements)
{}

class_tensors_infinite &class_tensors_infinite::operator=(const class_tensors_infinite &other) {
    // check for self-assignment
    if(this != &other) {
        state = std::make_unique<class_state_infinite>(*other.state);
        model = std::make_unique<class_model_infinite>(*other.model);
        edges = std::make_unique<class_edges_infinite>(*other.edges);
        measurements = other.measurements;
    }
    return *this;
}



void class_tensors_infinite::initialize(ModelType model_type_) {
    state->initialize(model_type_);
    model->initialize(model_type_);
    edges->initialize();
}

bool class_tensors_infinite::is_real() const { return state->is_real() and model->is_real() and edges->is_real(); }
bool class_tensors_infinite::has_nan() const { return state->has_nan() or model->has_nan() or edges->has_nan(); }
void class_tensors_infinite::assert_validity() const {
    state->assert_validity();
    model->assert_validity();
    edges->assert_validity();
}

void class_tensors_infinite::update_mps(const Eigen::Tensor<Scalar,3> & twosite_tensor) {
    state->set_mps(twosite_tensor);
    clear_cache();
    clear_measurements();
}

void class_tensors_infinite::enlarge() {
    tools::infinite::env::enlarge_edges(*state,*model,*edges);
    state->swap_AB();
}


void class_tensors_infinite::do_all_measurements() const {
    state->do_all_measurements();
    tools::infinite::measure::do_all_measurements(*this);
}

void class_tensors_infinite::clear_measurements() const {
    state->clear_measurements();
    measurements = tensors_measure_infinite();
}

void class_tensors_infinite::clear_cache() const {
    state->clear_cache();
    model->clear_cache();
}