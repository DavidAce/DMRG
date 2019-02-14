//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_INFINITE_DMRG_H
#define DMRG_CLASS_INFINITE_DMRG_H


#include "class_algorithm_base.h"
class class_table_dmrg;

/*!
 * \brief Class that runs the infinite DMRG algorithm.
 */
class class_iDMRG : public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_iDMRG(std::shared_ptr<class_hdf5_file> hdf5_);

    std::unique_ptr<class_hdf5_table<class_table_dmrg>> table_idmrg;

    void run()                                          override;
    void run_simulation()                               override;
    void run_preprocessing()                            override;
    void run_postprocessing()                           override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_state_to_file(bool force = false)        override;
    void store_progress_to_file(bool force = false)     override;
    long   chi_max()                                    override;
    int    num_sites()                                  override;
    int    store_freq()                                 override;
    int    print_freq()                                 override;
    bool   chi_grow()                                   override;
};


#endif //DMRG_CLASS_INFINITE_DMRG_H
