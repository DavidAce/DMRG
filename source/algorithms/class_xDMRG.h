//
// Created by david on 2018-02-09.
//

#ifndef DMRG_CLASS_EXITED_DMRG_H
#define DMRG_CLASS_EXITED_DMRG_H

#include "class_base_algorithm.h"
class class_table_dmrg;

/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_finite_chain_sweeper;
class class_xDMRG : public class_base_algorithm {
private:

public:
    //Inherit the constructor of class_base_algorithm
    using class_base_algorithm::class_base_algorithm;
    explicit class_xDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    std::unique_ptr<class_hdf5_table<class_table_dmrg>> table_xdmrg;
    int    max_length   ;
    int    max_sweeps   ;
    double r_strength = 0 ;  //Randomness strength for the random field.

    void run()                                          override;
    void check_convergence_overall()                    override;
    void initialize_constants()                         override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_table_entry_to_file()                    override;
    void single_xDMRG_step(long chi_max);
    void initialize_random_chain();
    auto find_greatest_overlap(Eigen::Tensor<Scalar,4> &theta);

};



#endif //DMRG_CLASS_EXITED_DMRG_H
