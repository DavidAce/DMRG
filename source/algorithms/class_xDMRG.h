//
// Created by david on 2018-02-09.
//

#ifndef DMRG_CLASS_EXITED_DMRG_H
#define DMRG_CLASS_EXITED_DMRG_H

#include "class_base_algorithm.h"

/*!
 * \brief Class that runs the excited-state DMRG algorithm.
 */

class class_environment_storage;
class class_xDMRG : public class_algorithm_base {
private:
    using Scalar = double;
    Eigen::Tensor<Scalar,4> state;
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_xDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    std::shared_ptr<class_environment_storage>  env_storage;
    void run() override;
    void single_xDMRG_step(long chi_max);
    void find_greatest_overlap();
    int  initialize_random_chain();
    int  env_storage_insert();
    void env_storage_overwrite_MPS();
    int  env_storage_move();
    void print_profiling()     override;
};



#endif //DMRG_CLASS_EXITED_DMRG_H
