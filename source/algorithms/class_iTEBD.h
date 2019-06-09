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
    explicit class_iTEBD(std::shared_ptr<h5pp::File> h5ppFile_);
//    explicit class_iTEBD(std::shared_ptr<class_hdf5_file> hdf5_);

    std::unique_ptr<class_hdf5_table<class_table_tebd>> table_itebd;
    std::vector<Eigen::Tensor<Scalar,4>> unitary_time_evolving_operators;
    Eigen::MatrixXcd h_evn;
    Eigen::MatrixXcd h_odd;



    void run()                                          override;
    void run_simulation()                               override;
    void run_preprocessing()                            override;
    void run_postprocessing()                           override;
    void single_TEBD_step(long chi_max);
//    void initialize_constants()                         override;
    void check_convergence_time_step();
    void check_convergence()                            override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_state_and_measurements_to_file(bool force = false)        override;
    void store_table_entry_progress(bool force = false)     override;
    long   chi_max()                                    override;
    size_t num_sites()                                  override;
    size_t store_freq()                                 override;
    size_t print_freq()                                 override;
    bool   chi_grow()                                   override;
};


#endif //DMRG_CLASS_IMAGINARY_TEBD_H
