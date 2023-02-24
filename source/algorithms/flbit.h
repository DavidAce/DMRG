#pragma once
#include "AlgorithmFinite.h"
#include <deque>
#include <qm/gate.h>
/*!
// * \brief Class that runs the finite LBIT algorithm.
 */

class StateFinite;
class flbit : public AlgorithmFinite {
    public:
    bool                                                ready_hamiltonian_gates      = false;
    bool                                                ready_hamiltonian_swap_gates = false;
    std::unique_ptr<StateFinite>                        state_lbit, state_lbit_init;
    std::vector<qm::Gate>                               ham_gates_1body, time_gates_1site;
    std::vector<qm::Gate>                               ham_gates_2body, time_gates_2site;
    std::vector<qm::Gate>                               ham_gates_3body, time_gates_3site;
    std::vector<qm::Gate>                               ham_gates_Lsite, time_gates_Lsite;
    std::vector<qm::SwapGate>                           ham_swap_gates_1body, time_swap_gates_1site;
    std::vector<qm::SwapGate>                           ham_swap_gates_2body, time_swap_gates_2site;
    std::vector<qm::SwapGate>                           ham_swap_gates_3body, time_swap_gates_3site;
    std::vector<std::vector<qm::Gate>>                  unitary_gates_2site_layers;
    std::vector<Eigen::Tensor<std::complex<double>, 4>> unitary_gates_mpo_layers;
    Eigen::Tensor<std::complex<double>, 1>              ledge, redge;
    std::vector<std::complex<double>>                   time_points;
    Eigen::Tensor<Scalar, 2>                            lbit_overlap; // The real-space support of the l-bits
    Eigen::Tensor<Scalar, 1>                            Upsi_ed;
    // Inherit the constructor of class_algorithm_base
    using AlgorithmFinite::AlgorithmFinite;
    explicit flbit(std::shared_ptr<h5pp::File> h5file_);
    void update_time_step();
    void resume() final;
    void run_task_list(std::deque<flbit_task> &task_list);
    void run_default_task_list() final;
    void run_preprocessing() final;
    void run_algorithm() final;
    void run_fes_analysis() final;
    void update_state() final;
    void check_convergence() final;
    void create_time_points();
    void create_hamiltonian_gates();
    void create_hamiltonian_swap_gates();
    void update_time_evolution_gates();
    void update_time_evolution_swap_gates();
    void create_unitary_circuit_gates();
    void transform_to_real_basis();
    void transform_to_lbit_basis();
    void write_to_file(StorageEvent storage_event = StorageEvent::ITER_STATE, CopyPolicy copy_policy = CopyPolicy::TRY) final;
    void print_status() final;
};
