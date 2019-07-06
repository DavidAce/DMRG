//
// Created by david on 2019-06-24.
//

#ifndef DMRG_CLASS_ALGORITHM_INFINITE_H
#define DMRG_CLASS_ALGORITHM_INFINITE_H
#include <algorithms/class_algorithm_base.h>

class class_table_dmrg;
class class_infinite_state;

class class_algorithm_infinite: public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_algorithm_infinite(
            std::shared_ptr<h5pp::File> h5ppFile_,
            std::string sim_name,
            SimulationType sim_type
    );
    std::unique_ptr<class_hdf5_table<class_table_dmrg>> table_dmrg;
    std::shared_ptr<class_infinite_state> state;


    virtual void run_simulation()         = 0;
    virtual void run_preprocessing() ;
    virtual void run_postprocessing();


    void enlarge_environment();
    void swap();
    void run()                                                          override;
    void compute_observables()                                          final;
    void clear_saturation_status()                                      override;
    void reset_to_random_state(const std::string parity)                final;
    void store_state_and_measurements_to_file(bool force = false)       final;
//    void store_table_entry_progress(bool force = false)                 final;
    void print_status_update()                                          final;
    void print_status_full()                                            final;
    void print_profiling()                                              final;
    void print_profiling_sim(class_tic_toc &t_parent)                   final;



    void check_convergence_variance_mpo(double threshold = quietNaN,double slope_threshold = quietNaN);
    void check_convergence_variance_ham(double threshold = quietNaN,double slope_threshold = quietNaN);
    void check_convergence_variance_mom(double threshold = quietNaN,double slope_threshold = quietNaN);
    void check_convergence_entg_entropy(double slope_threshold = quietNaN);

    std::list<bool>   B_mpo_vec; //History of saturation true/false
    std::list<double> V_mpo_vec; //History of variances
    std::list<int>    X_mpo_vec; //History of step numbers
    double V_mpo_slope = 0;

    std::list<bool>   B_ham_vec; //History of saturation true/false
    std::list<double> V_ham_vec;
    std::list<int>    X_ham_vec;
    double V_ham_slope = 0;

    std::list<bool>   B_mom_vec; //History of saturation true/false
    std::list<double> V_mom_vec;
    std::list<int>    X_mom_vec;
    double V_mom_slope = 0;

    std::list<bool>   BS_vec; //History of saturation true/false
    std::list<double> S_vec;
    std::list<int>    XS_vec;
    double S_slope = 0;
};


#endif //DMRG_CLASS_ALGORITHM_INFINITE_H
