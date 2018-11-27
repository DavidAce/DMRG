//
// Created by david on 2018-01-31.
//

#ifndef DMRG_CLASS_FINITE_DMRG_H
#define DMRG_CLASS_FINITE_DMRG_H

#include "class_algorithm_base.h"

class class_table_dmrg;

/*!
 * \brief Class that runs the finite DMRG algorithm.
 */

class class_finite_chain_sweeper;
class class_fDMRG : public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_fDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    std::unique_ptr<class_hdf5_table<class_table_dmrg>> table_fdmrg;
    std::unique_ptr<class_hdf5_table<class_table_finite_chain>> table_fdmrg_chain;

    int    max_length   ;
    int    max_sweeps   ;

    void run()                                          override;
    void check_convergence_all()                    override;
    void initialize_constants()                         override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_table_entry_to_file(bool force = false)  override;
    void store_chain_entry_to_file(bool force = false);
    void initialize_chain();

};



#endif //DMRG_CLASS_FINITE_DMRG_H
