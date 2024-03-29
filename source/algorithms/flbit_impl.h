#pragma once
#include "AlgorithmFinite.h"
#include "qm/gate.h"
#include "qm/lbit.h"
#include <deque>

/*!
// * \brief Class that runs the finite LBIT algorithm.
 */

class StateFinite;

namespace flbit_impl {
    std::pair<StateFinite, AlgorithmStatus> update_state(const size_t time_index, cplx_t time_point, const StateFinite &state_lbit_init,
                                                         const std::vector<std::vector<qm::SwapGate>> &gates_tevo,
                                                         const std::vector<std::vector<qm::Gate>> &unitary_circuit, const AlgorithmStatus &status_init);
    std::vector<std::vector<qm::SwapGate>>  get_time_evolution_gates(const cplx_t &time_point, const std::vector<std::vector<qm::SwapGate>> &ham_swap_gates);
    StateFinite                             time_evolve_lbit_state(const StateFinite &state_lbit_init, const std::vector<std::vector<qm::SwapGate>> &gates_tevo,
                                                                   const AlgorithmStatus &status);
    StateFinite                             transform_to_real_basis(const StateFinite &state_lbit, const std::vector<std::vector<qm::Gate>> &unitary_circuit,
                                                                    const AlgorithmStatus &status);
    AlgorithmStatus check_convergence (const AlgorithmStatus & status);
    void print_status(const StateFinite &state_real, const AlgorithmStatus &status);

}