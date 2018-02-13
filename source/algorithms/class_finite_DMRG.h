//
// Created by david on 2018-01-31.
//

#ifndef DMRG_CLASS_FINITE_DMRG_H
#define DMRG_CLASS_FINITE_DMRG_H


#include "class_base_algorithm.h"
/*!
 * \brief Class that runs the finite DMRG algorithm.
 */
class class_environment_storage;
class class_finite_DMRG : public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_finite_DMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    long chi_max    = settings::fdmrg::chi_max;
    int  max_length = settings::fdmrg::max_length;
    int  print_freq = settings::fdmrg::print_freq;
    int  max_sweeps = settings::fdmrg::max_sweeps;
    int  direction  = 1;
    int  sweep      = 0;
    int  position   = 0;
    int  middle_of_chain = max_length/2;

    std::shared_ptr<class_environment_storage>  env_storage;
    void run() override;

    void env_storage_insert();
    void env_storage_load();
    void env_storage_overwrite_MPS();
    void env_storage_move();


    void print_status_full()   override;
    void print_status_update() override;

    void store_table_entry()   override;
    void print_profiling()     override;

};



#endif //DMRG_CLASS_FINITE_DMRG_H
