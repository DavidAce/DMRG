//
// Created by david on 2019-06-24.
//

#pragma once

#include <algorithms/class_algorithm_base.h>
#include <tensors/class_tensors_finite.h>
class class_state_finite;
class class_model_finite;
class class_edges_finite;

// class h5pp_table_measurements_finite;
class class_algorithm_finite : public class_algorithm_base {
    public:
    // Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_algorithm_finite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type);
    class_tensors_finite tensors; // State, model and edges

    size_t excited_state_number = 0; /*!< Keeps track of found excited states. Convention: "0" is the biased seed,
                                      * and others become unbiased by randomizing from state 0 (or any other).
                                      * Therefore, this number is only incremented when randomizing the current state*/

    // Control behavior when stuck
    size_t max_stuck_iters      = 20;  //  5;  /*!< If stuck for this many sweeps -> stop. */
    size_t max_saturation_iters = 100;  // 10;  /*!< If either var or ent saturated this long -> algorithm saturated: true */
    size_t min_saturation_iters = 1;   // 1;   /*!< Saturated at least this many iters before stopping */
    size_t min_converged_iters  = 2;           /*!< Converged at least this many iters before success */

    bool                    has_projected        = false;        /*!< True if projection has already been tried */
    bool                    has_expanded         = false;        /*!< True if full expansion has already been tried */
    bool                    has_damped           = false;        /*!< True if damping of hamiltonian parameters is ongoing */
    size_t                  chi_quench_steps     = 0;            /*!< Number of steps left doing chi-quenching */
    size_t                  num_chi_quenches     = 0;            /*!< Number of bond dimension quench trials that have occurred */
    size_t                  max_chi_quenches     = 2;            /*!< Maximum number of bond dimension quench trials allowed */
    long                    chi_lim_quench_ahead = 32;           /*!< Bond dimension during a quench */
    long                    chi_lim_quench_trail = 32;           /*!< Bond dimension during a quench */
    size_t                  num_perturbations    = 0;            /*!< Number of perturbation trials done */
    size_t                  max_perturbations    = 2;            /*!< Maximum number of perturbation trials allowed */
    size_t                  perturbation_steps   = 0;            /*!< Number of steps left doing perturbation of MPOs */
    size_t                  damping_steps        = 0;            /*!< Number of steps left doing disorder damping of MPOs */
    size_t                  num_dampings         = 0;            /*!< Number of damping trials done */
    size_t                  max_dampings         = 2;            /*!< Maximum number of damping trials allowed */
    size_t                  iter_discard         = 0;            /*!< Iteration when last discard occurred */
    size_t                  num_discards         = 0;            /*!< Counter for number of times discarding the smallest schmidt values */
    size_t                  max_discards         = 2;            /*!< Maximum number of times to discard the smallest schmidt values */
    size_t                  num_expansion_iters  = 0;            /*!< Counter for number of iterations of subspace expansion */
    size_t                  max_expansion_iters  = 16;           /*!< Maximum number of iterations doing subspace expansion */
    double                  max_expansion_alpha  = 0.2;          /*!< Maximum value of subspace expansion factor */
    std::optional<double>   sub_expansion_alpha = std::nullopt;  /*!< The current value of the subspace expansion factor */
    std::vector<double>     damping_exponents;                   /*!< Exponents for for the damping trials */
    std::optional<OptMode>  last_optmode  = std::nullopt;
    std::optional<OptSpace> last_optspace = std::nullopt;

    public:
    virtual void run_algorithm()     = 0;
    virtual void run_preprocessing() = 0; // Specific for each algorithm type
    virtual void run_postprocessing();
    virtual bool cfg_store_wave_function() = 0;
    virtual void resume()                  = 0;
    virtual void run_default_task_list()   = 0;
    void         try_projection(std::optional<std::string> target_sector = std::nullopt);
    void         try_full_expansion();
    void         try_discard_small_schmidt();
    void         try_bond_dimension_quench();
    void         try_hamiltonian_perturbation();
    void         try_disorder_damping();
    void         move_center_point(std::optional<long> num_moves = std::nullopt);
    void         reduce_mpo_energy();
    void         rebuild_mpo_squared();
    void         update_variance_max_digits(std::optional<double> energy = std::nullopt) final;
    void         update_bond_dimension_limit(std::optional<long> tmp_bond_limit = std::nullopt) final;
    void         randomize_model();
    void         run() final;
    void         clear_convergence_status() override;
    void         randomize_state(ResetReason reason, StateInit state_init, std::optional<StateInitType> state_type = std::nullopt,
                                 std::optional<std::string> sector = std::nullopt, std::optional<long> chi_lim = std::nullopt,
                                 std::optional<bool> use_eigenspinors = std::nullopt, std::optional<long> bitfield = std::nullopt,
                                 std::optional<double> svd_threshold = std::nullopt);

    void write_to_file(StorageReason storage_reason = StorageReason::CHECKPOINT, std::optional<CopyPolicy> copy_file = std::nullopt) override;
    void print_status_update() override;
    void print_status_full() final;
    void check_convergence_variance(std::optional<double> threshold = std::nullopt, std::optional<double> saturation_sensitivity = std::nullopt);
    void check_convergence_entg_entropy(std::optional<double> saturation_sensitivity = std::nullopt);
    void check_convergence_spin_parity_sector(const std::string & target_sector, double threshold = 1e-12);
    void setup_prefix(const StorageReason &storage_reason, StorageLevel &storage_level, const std::string &state_name, std::string &state_prefix,
                      std::string &model_prefix, std::vector<std::string> &table_prefxs);
    void write_to_file(StorageReason storage_reason, const class_state_finite &state, std::optional<CopyPolicy> copy_policy = std::nullopt);
    template<typename T>
    void write_to_file(StorageReason storage_reason, const T &data, const std::string &name, std::optional<CopyPolicy> copy_policy = std::nullopt);
    std::vector<double> var_mpo_step;               // History of energy variances (from mpo) at each step
    std::vector<double> var_mpo_iter;               // History of energy variances (from mpo) at each iteration
    std::vector<std::vector<double>> entropy_iter;  // History of entanglement entropies at each iteration

};
