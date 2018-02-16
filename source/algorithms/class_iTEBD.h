//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_IMAGINARY_TEBD_H
#define DMRG_CLASS_IMAGINARY_TEBD_H
#include "class_base_algorithm.h"

/*!
 * \brief Class that runs the imaginary TEBD algorithm.
 */
class class_iTEBD :public class_algorithm_base {
public:
    using class_algorithm_base::class_algorithm_base;
    explicit class_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_);
    long   chi_max      = settings::itebd::chi_max;
    int    max_steps    = settings::itebd::max_steps;
    double delta_t0     = settings::itebd::delta_t0;
    double delta_tmin   = settings::itebd::delta_tmin;
    int    print_freq   = settings::itebd::print_freq;
    int    suzuki_order = settings::itebd::suzuki_order;
    double delta_t      = delta_t0;
    int    iteration    = 0;
    double phys_time    = 0;

    double old_entropy;
    double new_entropy;

    void reduce_timestep();


    void run() override;

    void print_status_full() override;
    void print_status_update() override;

    void store_table_entry() override;
    void print_profiling()   override;
};


#endif //DMRG_CLASS_IMAGINARY_TEBD_H
