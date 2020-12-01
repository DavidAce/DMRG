#pragma once

class class_state_finite;
class class_model_finite;
class class_edges_finite;
class class_tensors_finite;
namespace tools::finite::print {
    extern void dimensions(const class_state_finite & state);
    extern void dimensions(const class_edges_finite & edges);
    extern void dimensions(const class_tensors_finite & tensors);
    extern void model(const class_model_finite & model);

}