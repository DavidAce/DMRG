//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_IMAGINARY_TEBD_H
#define DMRG_CLASS_IMAGINARY_TEBD_H
#include "class_base_algorithm.h"
class class_table_tebd;

/*!
 * \brief Class that runs the imaginary TEBD algorithm.
 */
class class_iTEBD :public class_base_algorithm {
public:
    using class_base_algorithm::class_base_algorithm;
    explicit class_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_);

    std::unique_ptr<class_hdf5_table<class_table_tebd>> table_itebd;

    int    max_steps    ;
    double phys_time = 0;
    double delta_t   = 0; //Make sure this one gets initialized to delta_t0!
    double delta_t0     ;
    double delta_tmin   ;
    int    suzuki_order ;
    void run()                                          override;
    void initialize_constants()                         override;
    void update_chi()                                   override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_table_entry_to_file()                    override;

};


#endif //DMRG_CLASS_IMAGINARY_TEBD_H
