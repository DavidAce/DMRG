#pragma once

#include "algorithms/AlgorithmBase.h"
#include "math/svd/config.h"
#include "measure/MeasurementsStateFinite.h"
#include "tensors/TensorsFinite.h"
#include "tools/finite/opt_meta.h"

class StateFinite;
class ModelFinite;
class EdgesFinite;
namespace tools::finite::opt {
    class opt_mps;
    struct OptMeta;
}

// class h5pp_table_measurements_finite;
class AlgorithmFinite : public AlgorithmBase {
    using OptMeta = tools::finite::opt::OptMeta;

    private:
    size_t                             iter_last_bond_reduce = 0;
    std::optional<std::vector<size_t>> sites_mps, sites_mpo; // Used when moving sites

    public:
    // Inherit the constructor of class_algorithm_base
    using AlgorithmBase::AlgorithmBase;
    explicit      AlgorithmFinite(OptRitz opt_ritz_, AlgorithmType algo_type);
    explicit      AlgorithmFinite(std::shared_ptr<h5pp::File> h5ppFile_, OptRitz opt_ritz_, AlgorithmType algo_type);
    TensorsFinite tensors; // State, model and edges

    size_t                   projected_iter = 0; /*!< The last iteration when projection was tried */
    std::optional<OptCost>   last_optcost   = std::nullopt;
    std::optional<OptAlgo>   last_optalgo   = std::nullopt;
    std::optional<OptSolver> last_optsolver = std::nullopt;

    public:
    virtual void          resume()                = 0;
    virtual void          run_default_task_list() = 0;
    void                  try_projection(std::optional<std::string> target_axis = std::nullopt);
    void                  set_parity_shift_mpo(std::optional<std::string> target_axis = std::nullopt);
    void                  set_parity_shift_mpo_squared(std::optional<std::string> target_axis = std::nullopt);
    void                  try_moving_sites();
    void                  expand_environment(EnvExpandMode envexpMode, EnvExpandSide envexpSide = EnvExpandSide::FORWARD,
                                             std::optional<double> alpha = std::nullopt, std::optional<svd::config> svd_cfg = std::nullopt);
    void                  move_center_point(std::optional<long> num_moves = std::nullopt);
    virtual void          set_energy_shift_mpo(); // We override this in xdmrg
    void                  rebuild_tensors();
    void                  update_precision_limit(std::optional<double> energy_upper_bound = std::nullopt) final;
    void                  update_bond_dimension_limit() final;
    void                  reduce_bond_dimension_limit(double rate, UpdatePolicy when, StorageEvent storage_event);
    void                  update_truncation_error_limit() final;
    void                  update_environment_expansion_alpha();
    void                  initialize_model();
    void                  run() final;
    void                  run_rbds_analysis();
    void                  run_rtes_analysis();
    void                  run_postprocessing() override;
    [[nodiscard]] OptMeta get_opt_meta();
    void                  clear_convergence_status() override;
    void                  initialize_state(ResetReason reason, StateInit state_init, std::optional<StateInitType> state_type = std::nullopt,
                                           std::optional<std::string> axis = std::nullopt, std::optional<bool> use_eigenspinors = std::nullopt,
                                           std::optional<std::string> pattern = std::nullopt, std::optional<long> bond_lim = std::nullopt,
                                           std::optional<double> trnc_lim = std::nullopt);

    void write_to_file(StorageEvent storage_event = StorageEvent::ITERATION, CopyPolicy copy_policy = CopyPolicy::TRY) override;
    void print_status() override;
    void print_status_full() final;
    void check_convergence() override;
    void check_convergence_variance(std::optional<double> threshold = std::nullopt, std::optional<double> saturation_sensitivity = std::nullopt);
    void check_convergence_entg_entropy(std::optional<double> saturation_sensitivity = std::nullopt);
    void check_convergence_spin_parity_sector(std::string_view target_axis, double threshold = 1e-8);
    void write_to_file(const StateFinite &state, const ModelFinite &model, const EdgesFinite &edges, StorageEvent storage_event,
                       CopyPolicy copy_policy = CopyPolicy::TRY);
    template<typename T>
    void write_to_file(const T &data, std::string_view name, StorageEvent storage_event, CopyPolicy copy_policy = CopyPolicy::TRY);

    struct log_entry {
        AlgorithmStatus     status;
        double              energy;
        double              variance;
        std::vector<double> entropies;
                            log_entry(const AlgorithmStatus &s, const TensorsFinite &t);
    };

    std::vector<log_entry> algorithm_history;
    std::vector<double>    var_mpo_step; // History of energy variances (from mpo) at each step
    double                 ene_latest = 0.0;
    double                 var_latest = 1.0;
    double                 var_envexp = 1.0;
    double                 ene_delta  = 0.0;
    double                 var_delta  = 0.0;
    double                 var_change = 0.0; // Variance change from normal optimization
};
