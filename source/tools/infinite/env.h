#pragma once
class class_state_infinite;
class class_model_infinite;
class class_edges_infinite;

namespace tools::infinite::env {
    extern void reset_edges(const class_state_infinite &state, const class_model_infinite &model, class_edges_infinite &edges);
    extern void enlarge_edges(const class_state_infinite &state, const class_model_infinite &model, class_edges_infinite &edges);
}
