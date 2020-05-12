#pragma once
#include <simulation/enums.h>
#include <string>

/* clang-format off */
class class_state_finite;
class class_model_finite;
namespace tools::finite::mpo {
    extern void initialize                 (class_model_finite & model, ModelType model_type, size_t num_sites, size_t position);
    extern void randomize                  (class_model_finite & model);
    extern void perturb_hamiltonian        (class_model_finite & model, double coupling_ptb, double field_ptb, PerturbMode PerturbMode);
    extern void damp_hamiltonian           (class_model_finite & model, double coupling_damp, double field_damp);
    extern void reduce_mpo_energy          (class_model_finite & model, const class_state_finite & state);
    extern void reduce_mpo_energy_multi    (class_model_finite & model, const class_state_finite & state);
    extern void reduce_mpo_energy_2site    (class_model_finite & model, const class_state_finite & state);
}
/* clang-format on */
