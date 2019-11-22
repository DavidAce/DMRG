//
// Created by david on 2019-06-24.
//

#pragma once

#include <algorithms/class_algorithm_base.h>

class class_h5table_measurements_finite;
class class_state_finite;
class class_h5table_buffer_dynamic;


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

    // Tables
    std::shared_ptr<class_h5table_buffer<class_h5table_measurements_finite>>  h5tbuf_measurements; //Written every sweep


    //MPS
    std::unique_ptr<class_state_finite>    state;
//    std::list<std::unique_ptr<class_state_finite>> state_champions; // We keep the best from each sweep
//    std::list<class_state_finite> state_champions; // We keep the best from each sweep
    // What happens when stuck this many iterations:
    // 1: direct, 2: subspace, 2: subspace, 4: update bond dim if possile, else stop
    size_t max_stuck_iters               = 4; // If stuck for this many sweeps -> stop simulation
    size_t min_stuck_iters               = 2; // If stuck for this many sweeps -> do subspace instead of direct
    //What happens when saturated this many iterations
    size_t min_saturation_iters          = 1; // If both var and ent saturated  this long -> got_stuck: true
    size_t max_saturation_iters          = 3; // If either var or ent saturated this long -> got_stuck: true
    bool   has_projected  = false;

public:

    virtual void run_simulation()                    = 0;
    virtual void run_preprocessing();
    virtual void run_postprocessing();
    virtual void single_DMRG_step(std::string ritz);
    virtual bool store_wave_function()               = 0;
    void move_center_point();
    void update_bond_dimension_limit(std::optional<long> tmp_bond_limit = std::nullopt)         final;
    void run()                                                                                  final;
    void clear_saturation_status()                                                              override;
    void reset_to_random_state(const std::string parity_sector = "random", int seed_state = -1) final;

    void write_state        (bool result = false)                                               final;
    void write_measurements (bool result = false)                                               final;
    void write_sim_status   (bool result = false)                                               final;
    void write_profiling    (bool result = false)                                               final;
    void copy_from_tmp      (bool result = false)                                               final;
    void print_status_update()                                                                  final;
    void print_status_full()                                                                    final;
    void check_convergence_variance(double threshold = quietNaN, double slope_threshold = quietNaN);
    void check_convergence_entg_entropy(double slope_threshold = quietNaN);
    void write_projection(const class_state_finite & state_projected, std::string parity_sector);


//    std::list<bool>   B_mpo_vec; //History of saturation true/false
    std::list<double> V_mpo_vec; //History of variances
    std::list<int>    X_mpo_vec; //History of moves numbers
//    double V_mpo_slope = 0;
    std::list<double> V_mpo_slopes; //History of variance slopes

//    std::vector<std::list<bool>  > BS_mat;
    std::vector<std::list<double>> S_mat;
    std::vector<std::list<int>>    X_mat;
//    double S_slope = 0;
    std::list<double> S_slopes; //History of variance slopes

};

