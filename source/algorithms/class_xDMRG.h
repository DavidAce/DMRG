//
// Created by david on 2018-02-09.
//

#ifndef DMRG_CLASS_EXITED_DMRG_H
#define DMRG_CLASS_EXITED_DMRG_H

#include "class_algorithm_base.h"
class class_table_dmrg;

/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_finite_chain_sweeper;
class class_xDMRG : public class_algorithm_base {
private:
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_xDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    std::unique_ptr<class_hdf5_table<class_table_dmrg>> table_xdmrg;
    std::unique_ptr<class_hdf5_table<class_table_finite_chain>> table_xdmrg_chain;

    int    max_length   ;
    int    max_sweeps   ;
    double r_strength = 0 ;  //Randomness strength for the random field.

    //Energy ranges
    double energy_min = 0;
    double energy_max = 0;
    double energy_mid = 0;
    double energy_at_site = 0;

    void run()                                          override;
    void check_convergence_overall()                    override;
    void initialize_constants()                         override;
    void print_profiling()                              override;
    void print_profiling_sim(class_tic_toc &t_parent)   override;
    void store_table_entry_to_file()                    override;
    void store_chain_entry_to_file();
    void single_xDMRG_step(long chi_max, double energy_target = 0);
    void initialize_chain();
    void reset_chain_mps_to_random_product_state();
    void set_random_fields_in_chain_mpo();
    void find_energy_range();
    Eigen::Tensor<Scalar,4> find_state_with_greatest_overlap_part_diag(Eigen::Tensor<Scalar, 4> &theta, double energy_target = 0);
    Eigen::Tensor<Scalar,4> find_state_with_greatest_overlap_full_diag(Eigen::Tensor<Scalar, 4> &theta, double energy_target = 0);
};



#endif //DMRG_CLASS_EXITED_DMRG_H
