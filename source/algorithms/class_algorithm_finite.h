//
// Created by david on 2019-06-24.
//

#ifndef DMRG_CLASS_ALGORITHM_FINITE_H
#define DMRG_CLASS_ALGORITHM_FINITE_H

#include <algorithms/class_algorithm_base.h>
class class_finite_state;


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


    size_t min_saturation_iters          = 0;
    size_t max_saturation_iters          = 0;


    virtual void run_simulation()           = 0;
    virtual void run_preprocessing();
    virtual void run_postprocessing();
    virtual void single_DMRG_step(std::string ritz);
    virtual bool store_wave_function()               = 0;

    void move_center_point();
    void run()                                                              final;
    void compute_observables()                                              final;
    void clear_saturation_status()                                          override;
    void reset_to_random_state(const std::string symmetry)                    final;
    void write_measurements(bool force = false)                             final;
    void write_state(bool force = false)                                    final;
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

    std::vector<std::list<bool>  > BS_mat;
    std::vector<std::list<double>> S_mat;
    std::vector<std::list<int>>    XS_mat;
    double S_slope = 0;

};

#endif //DMRG_CLASS_ALGORITHM_FINITE_H
