//
// Created by david on 2019-06-24.
//

#pragma once

#include <algorithms/class_algorithm_base.h>
class table_measurements_finite;
class class_state_finite;


class class_algorithm_finite: public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_algorithm_finite(
            std::shared_ptr<h5pp::File> h5ppFile_,
            std::string sim_name,
            SimulationType sim_type,
            size_t num_sites
            );
    // Logs
    std::shared_ptr<class_hdf5_log<table_measurements_finite>>  measurements_logger; //Written every sweep
    std::shared_ptr<class_hdf5_log<table_measurements_finite>>  measurements_result; //Written when we update bond dimension or finish


    //MPS
    std::unique_ptr<class_state_finite>    state,state_backup;


//    size_t min_stuck_iters               = 2;
    size_t max_stuck_iters               = 4;
    size_t min_stuck_iters               = 2;
    size_t min_saturation_iters          = 1;
    size_t max_saturation_iters          = 3;
    bool   has_projected  = false;

private:
    void write_log_sim_status ();
    void write_log_measurement();
    void write_log_profiling();
public:

    virtual void run_simulation()                    = 0;
    virtual void run_preprocessing();
    virtual void run_postprocessing();
    virtual void single_DMRG_step(std::string ritz);
    virtual bool store_wave_function()               = 0;
    void move_center_point();
    void update_bond_dimension_limit(std::optional<long> max_bond_dim = std::nullopt)           final;
    void run()                                                                                  final;
    void compute_observables()                                                                  final;
    void clear_saturation_status()                                                              override;
    void reset_to_random_state(const std::string parity_sector = "random", int seed_state = -1) final;
    void backup_best_state(const class_state_finite &state);
    void check_convergence_variance(double threshold = quietNaN, double slope_threshold = quietNaN);
    void check_convergence_entg_entropy(double slope_threshold = quietNaN);

    void write_measurements(bool force = false)                                                 final;
    void write_state(bool force = false)                                                        final;
    void write_status(bool force = false)                                                       final;
    void write_logs(bool force = false)                                                         final;
    void write_results()                                                                        final;
    void print_status_update()                                                                  final;
    void print_status_full()                                                                    final;

    std::list<bool>   B_mpo_vec; //History of saturation true/false
    std::list<double> V_mpo_vec; //History of variances
    std::list<int>    X_mpo_vec; //History of moves numbers
    double V_mpo_slope = 0;

    std::vector<std::list<bool>  > BS_mat;
    std::vector<std::list<double>> S_mat;
    std::vector<std::list<int>>    XS_mat;
    double S_slope = 0;

};

