//
// Created by david on 2018-01-31.
//

#pragma once
#include "class_algorithm_finite.h"
#include <physics/class_quantum_gates.h>
#include <deque>
/*!
// * \brief Class that runs the finite LBIT algorithm.
 */

class class_state_finite;
class class_flbit : public class_algorithm_finite {
    public:
    std::unique_ptr<class_state_finite> state_lbit;
    std::vector<qm::Gate>               ham_gates_1body, time_gates_1site;
    std::vector<qm::Gate>               ham_gates_2body, time_gates_2site;
    std::vector<qm::Gate>               ham_gates_3body, time_gates_3site;
    std::vector<qm::Gate>               ham_gates_Lsite, time_gates_Lsite;
    std::vector<std::vector<qm::Gate>>  unitary_gates_2site_layers;
    std::vector<std::complex<double>>   time_points;
    Eigen::Tensor<Scalar,2>             lbit_overlap; // The real-space support of the l-bits
    Eigen::Tensor<Scalar, 1> Upsi_ed;
    // Inherit the constructor of class_algorithm_base
    using class_algorithm_finite::class_algorithm_finite;
    explicit class_flbit(std::shared_ptr<h5pp::File> h5pp_file_);
    void   single_flbit_step();
    void   update_time_step();
    void   resume() final;
    void   run_task_list(std::deque<flbit_task> &task_list);
    void   run_default_task_list() final;
    void   run_preprocessing() final;
    void   run_algorithm() final;
    void   check_convergence() final;
    void   create_time_points();
    void   create_hamiltonian_gates();
    void   create_time_evolution_gates();
    void   create_lbit_transform_gates();
    void   transform_to_real_basis();
    void   transform_to_lbit_basis();
    void   write_to_file(StorageReason storage_reason = StorageReason::CHECKPOINT, std::optional<CopyPolicy> copy_file = std::nullopt) final;
    bool   cfg_algorithm_is_on() final;
    long   cfg_chi_lim_max() final;
    size_t cfg_print_freq() final;
    bool   cfg_chi_lim_grow() final;
    long   cfg_chi_lim_init() final;
    bool   cfg_store_wave_function() final;
    void   print_status_update() final;
};
