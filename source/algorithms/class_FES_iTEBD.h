//
// Created by david on 2018-01-31.
//

#ifndef DMRG_FES_TEBD_H
#define DMRG_FES_TEBD_H


#include "class_base_algorithm.h"
/*!
 * \brief Class that studies Finite-entanglement scaling using the imaginary TEBD algorithm.
 */
class class_FES_iTEBD : public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_FES_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_);
    long   chi_min    = settings::fes_itebd::chi_min;
    long   chi_max    = settings::fes_itebd::chi_max;
    long   chi_num    = settings::fes_itebd::chi_num;
    int    max_steps  = settings::fes_itebd::max_steps;
    double delta_t0   = settings::fes_itebd::delta_t0;
    double delta_tmin = settings::fes_itebd::delta_tmin;
    int    print_freq = settings::fes_itebd::print_freq;
    int    suzuki_order = settings::fes_itebd::suzuki_order;
    double delta_t    = delta_t0;
    int    iteration  = 0;

    void run() override;
    void run2();

    void print_status_full()   override;
    void print_status_update() override;

    void store_table_entry()   override;
    void print_profiling()     override;
};
#endif //DMRG_FES_TEBD_H
