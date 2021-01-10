//
// Created by david on 2020-05-14.
//

#include "class_tensors_infinite.h"
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/common/log.h>
#include <tools/infinite/env.h>
#include <tools/infinite/measure.h>
#include <tools/infinite/mps.h>

class_tensors_infinite::class_tensors_infinite()
    : state(std::make_unique<class_state_infinite>()), model(std::make_unique<class_model_infinite>()), edges(std::make_unique<class_edges_infinite>()) {
    tools::log->trace("Constructing class_tensors_infinite");
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
class_tensors_infinite::~class_tensors_infinite()                              = default;            // default dtor
class_tensors_infinite::class_tensors_infinite(class_tensors_infinite &&other) = default;            // default move ctor
class_tensors_infinite &class_tensors_infinite::operator=(class_tensors_infinite &&other) = default; // default move assign

class_tensors_infinite::class_tensors_infinite(const class_tensors_infinite &other)
    : state(std::make_unique<class_state_infinite>(*other.state)), model(std::make_unique<class_model_infinite>(*other.model)),
      edges(std::make_unique<class_edges_infinite>(*other.edges)), measurements(other.measurements) {}

class_tensors_infinite &class_tensors_infinite::operator=(const class_tensors_infinite &other) {
    // check for self-assignment
    if(this != &other) {
        state        = std::make_unique<class_state_infinite>(*other.state);
        model        = std::make_unique<class_model_infinite>(*other.model);
        edges        = std::make_unique<class_edges_infinite>(*other.edges);
        measurements = other.measurements;
    }
    return *this;
}

void class_tensors_infinite::initialize(ModelType model_type_) {
    state->initialize(model_type_);
    model->initialize(model_type_);
    edges->initialize();
}

void class_tensors_infinite::randomize_model() {
    model->randomize();
    model->rebuild_mpo_squared();
    reset_edges();
}

void class_tensors_infinite::reset_to_random_product_state(const std::string &sector, long bitfield, bool use_eigenspinors) {
    eject_edges();
    state->clear_cache(); // Other caches can remain intact
    tools::infinite::mps::random_product_state(*state, sector, bitfield, use_eigenspinors);
    reset_edges();
}

bool class_tensors_infinite::is_real() const { return state->is_real() and model->is_real() and edges->is_real(); }
bool class_tensors_infinite::has_nan() const { return state->has_nan() or model->has_nan() or edges->has_nan(); }
void class_tensors_infinite::assert_validity() const {
    state->assert_validity();
    model->assert_validity();
    edges->assert_validity();
}

size_t class_tensors_infinite::get_length() const { return edges->get_length(); }

void class_tensors_infinite::reset_edges() {
    tools::infinite::env::reset_edges(*state, *model, *edges);
    clear_measurements();
}

void class_tensors_infinite::eject_edges() {
    edges->eject_edges();
    clear_measurements();
}

void class_tensors_infinite::merge_twosite_tensor(const Eigen::Tensor<Scalar, 3> &twosite_tensor, long chi_lim, std::optional<double> svd_threshold) {
    state->clear_cache();
    clear_measurements();
    tools::infinite::mps::merge_twosite_tensor(*state, twosite_tensor, chi_lim, svd_threshold);
    //    normalize_state(chi_lim, svd_threshold, NormPolicy::IFNEEDED);
}

void class_tensors_infinite::enlarge() {
    tools::infinite::env::enlarge_edges(*state, *model, *edges);
    state->swap_AB();
    clear_cache();
    clear_measurements();
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