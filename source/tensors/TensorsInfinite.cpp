#include "TensorsInfinite.h"
#include <tensors/edges/EdgesInfinite.h>
#include <tensors/model/ModelInfinite.h>
#include <tensors/site/mps/MpsSite.h>
#include <tensors/state/StateInfinite.h>
#include <tools/common/log.h>
#include <tools/infinite/env.h>
#include <tools/infinite/measure.h>
#include <tools/infinite/mps.h>

TensorsInfinite::TensorsInfinite()
    : state(std::make_unique<StateInfinite>()), model(std::make_unique<ModelInfinite>()), edges(std::make_unique<EdgesInfinite>()) {
    tools::log->trace("Constructing TensorsInfinite");
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
TensorsInfinite::~TensorsInfinite()                       = default;            // default dtor
TensorsInfinite::TensorsInfinite(TensorsInfinite &&other) = default;            // default move ctor
TensorsInfinite &TensorsInfinite::operator=(TensorsInfinite &&other) = default; // default move assign

TensorsInfinite::TensorsInfinite(const TensorsInfinite &other)
    : state(std::make_unique<StateInfinite>(*other.state)), model(std::make_unique<ModelInfinite>(*other.model)),
      edges(std::make_unique<EdgesInfinite>(*other.edges)), measurements(other.measurements) {}

TensorsInfinite &TensorsInfinite::operator=(const TensorsInfinite &other) {
    // check for self-assignment
    if(this != &other) {
        state        = std::make_unique<StateInfinite>(*other.state);
        model        = std::make_unique<ModelInfinite>(*other.model);
        edges        = std::make_unique<EdgesInfinite>(*other.edges);
        measurements = other.measurements;
    }
    return *this;
}

void TensorsInfinite::initialize(ModelType model_type_) {
    state->initialize(model_type_);
    model->initialize(model_type_);
    edges->initialize();
}

void TensorsInfinite::randomize_model() {
    model->randomize();
    model->rebuild_mpo_squared();
    reset_edges();
}

void TensorsInfinite::reset_to_random_product_state(std::string_view sector, long bitfield, bool use_eigenspinors) {
    eject_edges();
    state->clear_cache(); // Other caches can remain intact
    tools::infinite::mps::random_product_state(*state, sector, bitfield, use_eigenspinors);
    reset_edges();
}

bool TensorsInfinite::is_real() const { return state->is_real() and model->is_real() and edges->is_real(); }
bool TensorsInfinite::has_nan() const { return state->has_nan() or model->has_nan() or edges->has_nan(); }
void TensorsInfinite::assert_validity() const {
    state->assert_validity();
    model->assert_validity();
    edges->assert_validity();
}

size_t TensorsInfinite::get_length() const { return edges->get_length(); }

void TensorsInfinite::reset_edges() {
    tools::infinite::env::reset_edges(*state, *model, *edges);
    clear_measurements();
}

void TensorsInfinite::eject_edges() {
    edges->eject_edges();
    clear_measurements();
}

void TensorsInfinite::merge_twosite_tensor(const Eigen::Tensor<Scalar, 3> &twosite_tensor, long chi_lim, std::optional<svd::settings> svd_settings) {
    state->clear_cache();
    clear_measurements();
    tools::infinite::mps::merge_twosite_tensor(*state, twosite_tensor, chi_lim, svd_settings);
    //    normalize_state(chi_lim, svd_threshold, NormPolicy::IFNEEDED);
}

void TensorsInfinite::enlarge() {
    tools::infinite::env::enlarge_edges(*state, *model, *edges);
    state->swap_AB();
    clear_cache();
    clear_measurements();
}

void TensorsInfinite::do_all_measurements() const {
    state->do_all_measurements();
    tools::infinite::measure::do_all_measurements(*this);
}

void TensorsInfinite::clear_measurements() const {
    state->clear_measurements();
    measurements = tensors_measure_infinite();
}

void TensorsInfinite::clear_cache() const {
    state->clear_cache();
    model->clear_cache();
}