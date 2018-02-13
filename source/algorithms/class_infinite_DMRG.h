//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_INFINITE_DMRG_H
#define DMRG_CLASS_INFINITE_DMRG_H



#include "class_base_algorithm.h"
/*!
 * \brief Class that runs the infinite DMRG algorithm.
 */
class class_infinite_DMRG : public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_infinite_DMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    long chi_max    = settings::idmrg::chi_max;
    int  max_length = settings::idmrg::max_length;
    int  print_freq = settings::idmrg::print_freq;
    int  iteration  = 0;

    void run() override;

    void print_status_full()   override;
    void print_status_update() override;

    void store_table_entry()   override;
    void print_profiling()     override;
};


#endif //DMRG_CLASS_INFINITE_DMRG_H
