//
// Created by david on 2018-02-09.
//

#ifndef DMRG_CLASS_EXITED_DMRG_H
#define DMRG_CLASS_EXITED_DMRG_H

#include "class_algorithm_base.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Core>
class class_table_finite_chain;
class class_table_dmrg;




/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_finite_chain_state;
class class_xDMRG : public class_algorithm_base {
private:

public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_);
    std::unique_ptr<class_hdf5_table<class_table_dmrg>>         table_xdmrg;
    std::unique_ptr<class_hdf5_table<class_table_finite_chain>> table_xdmrg_chain;

    enum class xDMRG_Mode {KEEP_BEST_OVERLAP,FULL_EIG_OPT,PARTIAL_EIG_OPT, DIRECT_OPT};
    int    min_saturation_length;
    int    max_saturation_length;
    bool   projected_during_saturation  = false;

    //Energy ranges


    void run()                                          override;
    void run_simulation()                               override;
    void run_preprocessing()                            override;
    void run_postprocessing()                           override;
    void check_convergence()                            override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_state_and_measurements_to_file(bool force = false)        override;
    void store_table_entry_progress(bool force = false)     override;
    void store_table_entry_site_state(bool force = false);
    void single_xDMRG_step();
    void initialize_chain();
    void find_energy_range();
    long   chi_max()                                    override;
    int    num_sites()                                  override;
    int    store_freq()                                 override;
    int    print_freq()                                 override;
    bool   chi_grow()                                   override;

};



#endif //DMRG_CLASS_EXITED_DMRG_H
