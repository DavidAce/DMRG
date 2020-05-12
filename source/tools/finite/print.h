#pragma once

class class_state_finite;
class class_model_finite;
class class_edges_finite;
namespace tools::finite::print {
    extern void print_system        (const class_state_finite & state,const class_model_finite & model,const class_edges_finite & edges);
    extern void print_model         (const class_state_finite & state);
}