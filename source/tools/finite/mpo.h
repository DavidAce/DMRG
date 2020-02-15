#pragma once
#include <simulation/enums.h>
#include <string>

class class_state_finite;
namespace tools::finite::mpo {
    extern void initialize                 (class_state_finite & state, size_t length, std::string model_type);
    extern void randomize                  (class_state_finite & state);
    extern void perturb_hamiltonian        (class_state_finite & state, double coupling_ptb, double field_ptb, PerturbMode PerturbMode);
    extern void damp_hamiltonian           (class_state_finite & state, double coupling_damp, double field_damp);
    extern void reduce_mpo_energy          (class_state_finite & state);
    extern void reduce_mpo_energy_multi    (class_state_finite & state);
    extern void reduce_mpo_energy_2site    (class_state_finite & state);
}
