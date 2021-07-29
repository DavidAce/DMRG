#pragma once
class StateInfinite;
class ModelInfinite;
class EdgesInfinite;

namespace tools::infinite::env {
    extern void reset_edges(const StateInfinite &state, const ModelInfinite &model, EdgesInfinite &edges);
    extern void enlarge_edges(const StateInfinite &state, const ModelInfinite &model, EdgesInfinite &edges);
}
