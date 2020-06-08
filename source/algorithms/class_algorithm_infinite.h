//
// Created by david on 2019-06-24.
//

#pragma once

#include <algorithms/class_algorithm_base.h>
#include <tensors/class_tensors_infinite.h>
class class_state_infinite;
class class_model_infinite;
class class_edges_infinite;

class class_algorithm_infinite : public class_algorithm_base {
    public:
    // Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_algorithm_infinite(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType sim_type);
    class_tensors_infinite tensors;
    //    std::shared_ptr<class_state_infinite> state;
    // Tables
    //    std::shared_ptr<class_h5table_buffer<class_h5table_measurements_infinite>>  h5tbuf_measurements; //Written every iteration

    virtual void run_simulation() = 0;
    virtual void run_preprocessing();
    virtual void run_postprocessing();

    /* clang-format off */
    void run()                                                                                               final;
    void update_truncation_limit()                                                                           final;
    void update_bond_dimension_limit(std::optional<long> max_bond_dim = std::nullopt)                        final;
    void randomize_state(ResetReason reason, std::optional<std::string> sector = std::nullopt,
                                       std::optional<long> bitfield = std::nullopt, std::optional<bool> use_eigenspinors = std::nullopt);
    void clear_convergence_status()                                                                           override;
    void write_to_file(StorageReason storage_reason = StorageReason::CHECKPOINT)                             final;
    void copy_from_tmp(StorageReason storage_reason = StorageReason::CHECKPOINT)                             final;
    void print_status_update()                                                                               final;
    void print_status_full()                                                                                 final;
    /* clang-format on */

    void check_convergence_variance_mpo(double threshold = quietNaN, double slope_threshold = quietNaN);
    void check_convergence_variance_ham(double threshold = quietNaN, double slope_threshold = quietNaN);
    void check_convergence_variance_mom(double threshold = quietNaN, double slope_threshold = quietNaN);
    void check_convergence_entg_entropy(double slope_threshold = quietNaN);

    std::vector<bool>   B_mpo_vec; // History of saturation true/false
    std::vector<double> V_mpo_vec; // History of variances
    std::vector<size_t> X_mpo_vec; // History of moves numbers
    double            V_mpo_slope = 0;

    std::vector<bool>   B_ham_vec; // History of saturation true/false
    std::vector<double> V_ham_vec;
    std::vector<size_t> X_ham_vec;
    double            V_ham_slope = 0;

    std::vector<bool>   B_mom_vec; // History of saturation true/false
    std::vector<double> V_mom_vec;
    std::vector<size_t> X_mom_vec;
    double            V_mom_slope = 0;

    std::vector<bool>   BS_vec; // History of saturation true/false
    std::vector<double> S_vec;
    std::vector<size_t> XS_vec;
    double            S_slope = 0;
};
