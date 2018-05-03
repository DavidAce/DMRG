//
// Created by david on 2018-02-09.
//

#ifndef DMRG_CLASS_EXITED_DMRG_H
#define DMRG_CLASS_EXITED_DMRG_H

#include "class_base_algorithm.h"

/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_finite_chain_storage;
class class_xDMRG : public class_base_algorithm {
private:
    Eigen::Tensor<Scalar,4> state;
    double randomness_strength = 0;
public:
    //Inherit the constructor of class_base_algorithm
    using class_base_algorithm::class_base_algorithm;
    explicit class_xDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    std::shared_ptr<class_finite_chain_storage>  env_storage;
    void run() override;
    void single_xDMRG_step(long chi_max);
    auto find_greatest_overlap();
    int  initialize_random_chain();
    int  env_storage_insert();
    void env_storage_overwrite_MPS();
    int  env_storage_move();
    void print_profiling()     override;
    void print_profiling_sim(class_tic_toc &t_parent) override;

};



#endif //DMRG_CLASS_EXITED_DMRG_H
