//
// Created by david on 2018-01-31.
//

#ifndef DMRG_CLASS_FINITE_DMRG_H
#define DMRG_CLASS_FINITE_DMRG_H


#include "class_base_algorithm.h"
/*!
 * \brief Class that runs the finite DMRG algorithm.
 */
class class_environment_storage;
class class_fDMRG : public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_fDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    int  middle_of_chain = max_length/2;

    std::shared_ptr<class_environment_storage>  env_storage;
    void run() override;
    int  initialize_chain();
    int  env_storage_insert();
    void env_storage_overwrite_MPS();
    int  env_storage_move();
    void print_profiling()     override;

};



#endif //DMRG_CLASS_FINITE_DMRG_H
