//
// Created by david on 2020-05-14.
//

#include "class_tensors_infinite.h"
#include <math/nmspc_math.h>
#include <tensors/edges/class_edges_infinite.h>
#include <tensors/edges/class_env_ene.h>
#include <tensors/edges/class_env_var.h>
#include <tensors/model/class_model_infinite.h>
#include <tensors/model/class_mpo_base.h>
#include <tensors/state/class_mps_2site.h>
#include <tensors/state/class_mps_site.h>
#include <tensors/state/class_state_infinite.h>
#include <tools/infinite/env.h>
#include <tools/infinite/measure.h>
#include <tools/infinite/mpo.h>
#include <tools/infinite/mps.h>

class_tensors_infinite::class_tensors_infinite() {
    state = std::make_unique<class_state_infinite>();
    model = std::make_unique<class_model_infinite>();
    edges = std::make_unique<class_edges_infinite>();
}

void class_tensors_infinite::initialize(ModelType model_type_) {
    state->set_chi_lim(2); // Can't call chi_init() <-- it's a pure virtual function
    if(state->has_nan()) throw std::runtime_error("State has NAN's before initializing it");
    tools::infinite::mps::initialize(*state, model_type_);
    tools::infinite::mpo::initialize(*model, model_type_);
    tools::infinite::mpo::randomize(*model);
    tools::infinite::env::initialize(*edges);
}

bool class_tensors_infinite::is_real() const { return state->is_real() and model->is_real() and edges->is_real(); }
bool class_tensors_infinite::has_nan() const { return state->has_nan() or model->has_nan() or edges->has_nan(); }
void class_tensors_infinite::assert_validity() const {
    state->assert_validity();
    model->assert_validity();
    edges->assert_validity();
}

void class_tensors_infinite::do_all_measurements() const {
    tools::infinite::measure::do_all_measurements(*state);
    tools::infinite::measure::do_all_measurements(*this);
}

void class_tensors_infinite::clear_measurements() const {
    measurements = tensors_measure_infinite();
    state->clear_measurements();
}

void class_tensors_infinite::clear_cache() const {
    state->clear_cache();
    model->clear_cache();
}