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
    int    max_steps  ;
    void run()                                          override;
    void initialize_constants()                         override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_table_entry_to_file()                    override;

};


#endif //DMRG_CLASS_INFINITE_DMRG_H
