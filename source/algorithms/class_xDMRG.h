//
// Created by david on 2018-02-09.
//

#pragma once

#include "class_algorithm_finite.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Core>


/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_state_finite;
class class_xDMRG : public class_algorithm_finite {
private:
    double energy_window_growth_factor = 1.0;
    size_t force_overlap = 0;
    size_t num_chi_quenches = 0;

    public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_finite::class_algorithm_finite;
    explicit class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_);


    void find_energy_range();
    void inflate_initial_state();
    void single_xDMRG_step();
    void reset_to_random_state_in_energy_window(const std::string &parity_sector,bool inflate, std::string reason );
    void try_projection();
    void try_chi_quench();
    void run_preprocessing()                final;
    void run_simulation()                   final;
    void check_convergence()                final;
    bool   sim_on()                         final;
    long   chi_max()                        final;
    size_t num_sites()                      final;
    size_t write_freq()                     final;
    size_t print_freq()                     final;
    bool   chi_grow()                       final;
    long   chi_init()                       final;
    bool   store_wave_function()            final;

};



