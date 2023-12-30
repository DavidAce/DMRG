#pragma once

#include "algorithms/AlgorithmBase.h"
#include "measure/MeasurementsStateFinite.h"
#include "tensors/TensorsFinite.h"
class StateFinite;
class ModelFinite;
class EdgesFinite;

// class h5pp_table_measurements_finite;
class AlgorithmFinite : public AlgorithmBase {
    private:
    size_t                             iter_last_bond_reduce = 0;
    std::optional<std::vector<size_t>> sites_mps, sites_mpo; // Used when moving sites

    public:
    // Inherit the constructor of class_algorithm_base
    using AlgorithmBase::AlgorithmBase;
    explicit AlgorithmFinite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type);
    TensorsFinite tensors; // State, model and edges

    size_t excited_state_number = 0; /*!< Keeps track of found excited states. */

    size_t                   projected_iter = 0; /*!< The last iteration when projection was tried */
    size_t                   expanded_iter  = 0; /*!< The last iteration when expansion was tried */
    std::optional<OptMode>   last_optmode   = std::nullopt;
    std::optional<OptSolver> last_optspace  = std::nullopt;

    public:
    virtual void run_fes_analysis()      = 0;
    virtual void resume()                = 0;
    virtual void run_default_task_list() = 0;
    void         try_projection(std::optional<std::string> target_sector = std::nullopt);
    void         try_parity_shift();
    void         try_parity_sep();
    void         try_moving_sites();
    void         move_center_point(std::optional<long> num_moves = std::nullopt);
    void         shift_mpo_energy();
    void         update_variance_max_digits(std::optional<double> energy = std::nullopt) final;
    void         update_bond_dimension_limit() final;
    void         reduce_bond_dimension_limit(double rate, UpdateWhen when, StorageEvent storage_event);
    void         update_truncation_error_limit() final;
    void         update_expansion_factor_alpha();
    void         initialize_model();
    void         run() final;
    void         run_postprocessing() override;
    void         clear_convergence_status() override;
    void         initialize_state(ResetReason reason, StateInit state_init, std::optional<StateInitType> state_type = std::nullopt,
                                 std::optional<std::string> axis = std::nullopt, std::optional<bool> use_eigenspinors = std::nullopt,
                                 std::optional<std::string> pattern = std::nullopt, std::optional<long> bond_lim = std::nullopt,
                                 std::optional<double> trnc_lim = std::nullopt);

    void write_to_file(StorageEvent storage_event = StorageEvent::ITERATION, CopyPolicy copy_policy = CopyPolicy::TRY) override;
    void print_status() override;
    void print_status_full() final;
    void check_convergence_variance(std::optional<double> threshold = std::nullopt, std::optional<double> saturation_sensitivity = std::nullopt);
    void check_convergence_entg_entropy(std::optional<double> saturation_sensitivity = std::nullopt);
    void check_convergence_spin_parity_sector(std::string_view target_sector, double threshold = 1e-8);
    void write_to_file(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges, StorageEvent storage_event,
                       CopyPolicy copy_policy = CopyPolicy::TRY);
    template<typename T>
    void write_to_file(const T &data, std::string_view name, StorageEvent storage_event, CopyPolicy copy_policy = CopyPolicy::TRY);

    struct log_entry {
        AlgorithmStatus     status;
        double              variance;
        std::vector<double> entropies;
        log_entry(const AlgorithmStatus &s, const TensorsFinite &t);
    };

    std::vector<log_entry> algorithm_history;
    std::vector<double>    var_mpo_step; // History of energy variances (from mpo) at each step
};
