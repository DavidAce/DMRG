//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_IMAGINARY_TEBD_H
#define DMRG_CLASS_IMAGINARY_TEBD_H
#include "class_algorithm_base.h"
class class_table_tebd;

/*!
 * \brief Class that runs the imaginary TEBD algorithm.
 */
class class_iTEBD :public class_algorithm_base {
public:
    using class_algorithm_base::class_algorithm_base;
    explicit class_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_);

    std::unique_ptr<class_hdf5_table<class_table_tebd>> table_itebd;
    std::vector<Eigen::Tensor<Scalar,4>> unitary_time_evolving_operators;
    Eigen::MatrixXcd h_evn;
    Eigen::MatrixXcd h_odd;

    int    max_steps    ;
    double phys_time = 0;
    double delta_t   = 0; //Make sure this one gets initialized to delta_t0!
    double delta_t0     ;
    double delta_tmin   ;
    int    suzuki_order ;
    bool   time_step_has_converged;

    void run()                                          override;
    void single_TEBD_step(long chi_max);
    void initialize_constants()                         override;
    void check_convergence_time_step();
    void check_convergence_all()                        override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_table_entry_to_file(bool force = false)  override;

};


#endif //DMRG_CLASS_IMAGINARY_TEBD_H
