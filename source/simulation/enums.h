//
// Created by david on 2019-10-16.
//

#pragma once
enum class SimulationType { iDMRG, fDMRG, xDMRG, iTEBD };
enum class StorageLevel : size_t { NONE, LIGHT, NORMAL, FULL };
enum class PerturbMode {
    PERCENTAGE,                // J_ptb = couplingPtb * J_rnd
    ABSOLUTE,                  // J_ptb = couplingPtb
    UNIFORM_RANDOM_PERCENTAGE, // J_ptb = std::random_uniform(-couplingPtb, couplingPtb) * J_rnd
    UNIFORM_RANDOM_ABSOLUTE,   // J_ptb = std::random_uniform(-couplingPtb, couplingPtb)
};