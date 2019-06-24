//
// Created by david on 2019-06-24.
//

#ifndef DMRG_CLASS_ALGORITHM_FINITE_H
#define DMRG_CLASS_ALGORITHM_FINITE_H

#include <algorithms/class_algorithm_base.h>

class class_finite_chain_state;
class class_table_finite_chain;
class class_table_dmrg;


class class_algorithm_finite: public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_algorithm_finite(
            std::shared_ptr<h5pp::File> h5ppFile_,
            std::string sim_name,
            SimulationType sim_type
            );

    //MPS
    std::shared_ptr<class_finite_chain_state>    state;
    std::unique_ptr<class_hdf5_table<class_table_dmrg>> table_dmrg;
    std::unique_ptr<class_hdf5_table<class_table_finite_chain>> table_dmrg_chain;


    size_t min_saturation_length          = 0;
    size_t max_saturation_length          = 0;


    virtual void single_DMRG_step()                     = 0;


    void store_table_entry_site_state(bool force = false);
    void store_chain_entry_to_file(bool force = false);
    void initialize_chain();

    void run()                                                              final;
    void run_postprocessing()                                               final;
    void compute_observables()                                              final;
    void check_convergence_entg_entropy(double slope_threshold = quietNaN)  final;
    void store_state_and_measurements_to_file(bool force = false)           final;
    void store_table_entry_progress(bool force = false)                     final;
    void print_profiling()                                                  final;
    void print_profiling_sim(class_tic_toc &t_parent)                       final;



};

#endif //DMRG_CLASS_ALGORITHM_FINITE_H
