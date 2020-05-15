//
// Created by david on 2019-06-24.
//

#pragma once

#include <algorithms/class_algorithm_base.h>
#include <tensors/class_tensors_finite.h>
class class_state_finite;
class class_model_finite;
class class_edges_finite;

//class h5pp_table_measurements_finite;
class class_algorithm_finite : public class_algorithm_base {
    public:
    // Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_algorithm_finite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type);
    ~class_algorithm_finite(); // Read comment on implementation
//    std::unique_ptr<class_state_finite> state;// The finite chain state
//    std::unique_ptr<class_model_finite> model;// The finite chain model
//    std::unique_ptr<class_edges_finite> edges;// The finite chain edges
    class_tensors_finite tensors; // State, model and edges

    size_t state_number = 0;

    // Control behavior when stuck
    size_t              min_stuck_iters      = 1;     /*!< If stuck for this many sweeps -> do subspace instead of direct */
    size_t              max_stuck_iters      = 3;     /*!< If stuck for this many sweeps -> try stuff, or stop. */
    size_t              min_saturation_iters = 1;     /*!< If both var and ent saturated  this long -> got_stuck: true */
    size_t              max_saturation_iters = 3;     /*!< If either var or ent saturated this long -> got_stuck: true */
    bool                has_projected        = false; /*!< True if projection has already been tried */
    bool                has_damped           = false; /*!< True if damping of hamiltonian parameters is ongoing */
    size_t              chi_quench_steps     = 0;     /*!< Number of steps left doing chi-quenching */
    size_t              num_chi_quenches     = 0;     /*!< Number of bond dimension quench trials that have occurred */
    size_t              max_chi_quenches     = 2;     /*!< Maximum number of bond dimension quench trials allowed */
    size_t              chi_lim_quench_ahead = 32;    /*!< Bond dimension during a quench */
    size_t              chi_lim_quench_trail = 32;    /*!< Bond dimension during a quench */
    size_t              num_perturbations    = 0;     /*!< Number of perturbation trials done */
    size_t              max_perturbations    = 2;     /*!< Maximum number of perturbation trials allowed */
    size_t              perturbation_steps   = 0;     /*!< Number of steps left doing perturbation of MPOs */
    size_t              damping_steps        = 0;     /*!< Number of steps left doing disorder damping of MPOs */
    size_t              num_dampings         = 0;     /*!< Number of damping trials done */
    size_t              max_dampings         = 2;     /*!< Maximum number of damping trials allowed */
    std::vector<double> damping_exponents;            /*!< Exponents for for the damping trials */

    public:
    virtual void run_simulation() = 0;
    virtual void run_preprocessing();
    virtual void run_postprocessing();
    virtual bool store_wave_function() = 0;
    void         try_projection(class_state_finite & state);
    void         try_bond_dimension_quench(class_state_finite & state);
    void         try_hamiltonian_perturbation(class_state_finite & state);
    void         try_disorder_damping(class_model_finite & model);
    void         move_center_point(std::optional<size_t> num_moves = std::nullopt);
    void         update_truncation_limit() final;
    void         update_bond_dimension_limit(std::optional<long> tmp_bond_limit = std::nullopt) final;
    void         run() final;
    void         run_old() final;
    void         clear_saturation_status() override;
    void         reset_to_random_product_state(const std::string &parity_sector = "random") final;
    void         reset_to_random_current_state(std::optional<double> chi_lim = std::nullopt) final;
    void         reset_to_initial_state() final;
    void         write_to_file(StorageReason storage_reason = StorageReason::CHECKPOINT) final;
    void         copy_from_tmp(StorageReason storage_reason = StorageReason::CHECKPOINT) final;
    void         print_status_update() final;
    void         print_status_full() final;
    void         check_convergence_variance(double threshold = quietNaN, double slope_threshold = quietNaN);
    void         check_convergence_entg_entropy(double slope_threshold = quietNaN);
    void         write_to_file(StorageReason storage_reason, const class_state_finite &state);
//
    std::list<double> V_mpo_vec;    // History of variances
    std::list<size_t> X_mpo_vec;    // History of moves numbers
    std::list<double> V_mpo_slopes; // History of variance slopes

    std::vector<std::list<double>> S_mat;
    std::vector<std::list<size_t>> X_mat;
    std::list<double>              S_slopes; // History of variance slopes
};
