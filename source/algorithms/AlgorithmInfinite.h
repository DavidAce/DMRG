#pragma once
#include <algorithms/AlgorithmBase.h>
#include <tensors/TensorsInfinite.h>
class StateInfinite;
class ModelInfinite;
class EdgesInfinite;

class AlgorithmInfinite : public AlgorithmBase {
    public:
    // Inherit the constructor of class_algorithm_base
    using AlgorithmBase::AlgorithmBase;
    explicit AlgorithmInfinite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type);
    TensorsInfinite tensors;

    virtual void run_simulation() = 0;
    virtual void run_preprocessing();
    virtual void run_postprocessing();

    /* clang-format off */
    void run()                                                                                               final;
    void update_variance_max_digits(std::optional<double> energy = std::nullopt)                             final;
    void update_bond_dimension_limit()                                                                       final;
    void update_truncation_error_limit()                                                                     final;
    void randomize_model();
    void randomize_state(ResetReason reason,
                         std::optional<std::string> sector = std::nullopt,
                         std::optional<bool> use_eigenspinors = std::nullopt, std::optional<size_t> bitfield = std::nullopt);

    void clear_convergence_status()                                                                          override;
    void write_to_file(StorageReason storage_reason = StorageReason::CHECKPOINT,
                       std::optional<CopyPolicy> copy_policy = std::nullopt)                                 final;
    void print_status_update()                                                                               final;
    void print_status_full()                                                                                 final;
    /* clang-format on */

    void check_convergence_variance_mpo(std::optional<double> threshold = std::nullopt, std::optional<double> sensitivity = std::nullopt);
    void check_convergence_variance_ham(std::optional<double> threshold = std::nullopt, std::optional<double> sensitivity = std::nullopt);
    void check_convergence_variance_mom(std::optional<double> threshold = std::nullopt, std::optional<double> sensitivity = std::nullopt);
    void check_convergence_entg_entropy(std::optional<double> sensitivity = std::nullopt);

    std::vector<double> var_mpo_iter; // History of energy variances (from mpo) at each iteration
    std::vector<double> var_ham_iter;
    std::vector<double> var_mom_iter;
    std::vector<double> entropy_iter;
};
