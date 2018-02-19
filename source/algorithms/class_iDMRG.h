//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_INFINITE_DMRG_H
#define DMRG_CLASS_INFINITE_DMRG_H



#include "class_base_algorithm.h"
/*!
 * \brief Class that runs the infinite DMRG algorithm.
 */

class class_iDMRG : public class_algorithm_base {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_base::class_algorithm_base;
    explicit class_iDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    void run() override;
    void print_profiling() override;
};


#endif //DMRG_CLASS_INFINITE_DMRG_H
