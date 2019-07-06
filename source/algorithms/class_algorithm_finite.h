//
// Created by david on 2019-06-24.
//

#ifndef DMRG_CLASS_ALGORITHM_FINITE_H
#define DMRG_CLASS_ALGORITHM_FINITE_H

#include <algorithms/class_algorithm_base.h>
class class_finite_state;
class class_table_finite_chain;
class class_table_dmrg;


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

    //MPS
    std::unique_ptr<class_finite_state>    state;
    std::unique_ptr<class_hdf5_table<class_table_dmrg>> table_dmrg;
    std::unique_ptr<class_hdf5_table<class_table_finite_chain>> table_dmrg_chain;


    size_t min_saturation_length          = 0;
    size_t max_saturation_length          = 0;


    virtual void run_simulation()           = 0;
    virtual void run_preprocessing();
    virtual void run_postprocessing();
    virtual void single_DMRG_step(std::string ritz);
    virtual bool store_wave_function()               = 0;

    void move_center_point();
    void store_table_entry_site_state(bool force = false);
    void store_chain_entry_to_file(bool force = false);
    void initialize_chain();



    void run()                                                              final;
    void compute_observables()                                              final;
    void clear_saturation_status()                                          override;
    void reset_to_random_state(const std::string parity)                    final;
    void store_state_and_measurements_to_file(bool force = false)           final;
//    void store_table_entry_progress(bool force = false)                     final;
    void print_status_update()                                              final;
    void print_status_full()                                                final;
    void print_profiling()                                                  final;
    void print_profiling_sim(class_tic_toc &t_parent)                       final;


    void check_convergence_variance(double threshold = quietNaN, double slope_threshold = quietNaN);
    void check_convergence_entg_entropy(double slope_threshold = quietNaN);
    std::list<bool>   B_mpo_vec; //History of saturation true/false
    std::list<double> V_mpo_vec; //History of variances
    std::list<int>    X_mpo_vec; //History of step numbers
    double V_mpo_slope = 0;

    std::list<bool>   BS_vec; //History of saturation true/false
    std::list<double> S_vec;
    std::list<int>    XS_vec;
    double S_slope = 0;
};

#endif //DMRG_CLASS_ALGORITHM_FINITE_H
