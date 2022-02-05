#pragma once

#include <algorithms/AlgorithmBase.h>
#include <measure/MeasurementsStateFinite.h>
#include <tensors/TensorsFinite.h>
class StateFinite;
class ModelFinite;
class EdgesFinite;

// class h5pp_table_measurements_finite;
class AlgorithmFinite : public AlgorithmBase {
    public:
    // Inherit the constructor of class_algorithm_base
    using AlgorithmBase::AlgorithmBase;
    explicit AlgorithmFinite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type);
    TensorsFinite tensors; // State, model and edges

    size_t excited_state_number = 0; /*!< Keeps track of found excited states. */

    size_t                   projected_iter        = 0;  /*!< The last iteration when projection was tried */
    size_t                   expanded_iter         = 0;  /*!< The last iteration when expansion was tried */
    size_t                   bond_quench_steps     = 0;  /*!< Number of steps left doing bond dimension quenching */
    size_t                   num_bond_quenches     = 0;  /*!< Number of bond dimension quench trials that have occurred */
    size_t                   max_bond_quenches     = 2;  /*!< Maximum number of bond dimension quench trials allowed */
    long                     bond_lim_quench_ahead = 32; /*!< Bond dimension during a quench */
    long                     bond_lim_quench_trail = 32; /*!< Bond dimension during a quench */
    std::optional<OptMode>   last_optmode          = std::nullopt;
    std::optional<OptSolver> last_optspace         = std::nullopt;

    public:
    virtual void run_algorithm()     = 0;
    virtual void run_fes_analysis()  = 0;
    virtual void run_preprocessing() = 0; // Specific for each algorithm type
    virtual void run_postprocessing();
    virtual void resume()                = 0;
    virtual void run_default_task_list() = 0;
    void         try_projection(std::optional<std::string> target_sector = std::nullopt);
    void         try_full_expansion();
    void         try_bond_dimension_quench();
    void         move_center_point(std::optional<long> num_moves = std::nullopt);
    void         shift_mpo_energy();
    void         rebuild_mpo_squared();
    void         update_variance_max_digits(std::optional<double> energy = std::nullopt) final;
    void         update_bond_dimension_limit() final;
    void         reduce_bond_dimension_limit();
    void         update_expansion_factor_alpha();
    void         randomize_model();
    void         run() final;
    void         clear_convergence_status() override;
    void         randomize_state(ResetReason reason, StateInit state_init, std::optional<StateInitType> state_type = std::nullopt,
                                 std::optional<std::string> sector = std::nullopt, std::optional<long> bond_limit = std::nullopt,
                                 std::optional<bool> use_eigenspinors = std::nullopt, std::optional<long> bitfield = std::nullopt,
                                 std::optional<double> svd_threshold = std::nullopt);

    void write_to_file(StorageReason storage_reason = StorageReason::CHECKPOINT, std::optional<CopyPolicy> copy_file = std::nullopt) override;
    void print_status_update() override;
    void print_status_full() final;
    void check_convergence_variance(std::optional<double> threshold = std::nullopt, std::optional<double> saturation_sensitivity = std::nullopt);
    void check_convergence_entg_entropy(std::optional<double> saturation_sensitivity = std::nullopt);
    void check_convergence_spin_parity_sector(std::string_view target_sector, double threshold = 1e-8);
    void write_to_file(StorageReason storage_reason, const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges,
                       std::optional<CopyPolicy> copy_policy = std::nullopt);
    template<typename T>
    void write_to_file(StorageReason storage_reason, const T &data, std::string_view name, std::optional<CopyPolicy> copy_policy = std::nullopt);

    struct log_entry {
        AlgorithmStatus     status;
        double              variance;
        std::vector<double> entropies;
        log_entry(const AlgorithmStatus &s, const TensorsFinite &t);
    };

    std::vector<log_entry> algorithm_history;
    std::vector<double>    var_mpo_step; // History of energy variances (from mpo) at each step
};
