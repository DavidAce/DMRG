#pragma once
class class_state_infinite;
class class_model_infinite;
class class_tensors_infinite;
namespace tools::infinite::print {
    extern void dimensions(const class_tensors_infinite &tensors);         /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    extern void print_state_compact(const class_state_infinite &state); /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    extern void print_hamiltonians(const class_model_infinite &model);
}