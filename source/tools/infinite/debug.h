#pragma once
class class_state_infinite;
class class_model_infinite;
class class_edges_infinite;
class class_tensors_infinite;
namespace tools::infinite::debug {
    extern void check_integrity(const class_tensors_infinite &tensors);
    extern void check_integrity(const class_state_infinite &state);
    extern void check_integrity(const class_model_infinite &model);
    extern void check_integrity(const class_edges_infinite &edges);
}