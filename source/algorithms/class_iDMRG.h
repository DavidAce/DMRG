//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_INFINITE_DMRG_H
#define DMRG_CLASS_INFINITE_DMRG_H


#include "class_algorithm_infinite.h"
class class_table_dmrg;

/*!
 * \brief Class that runs the infinite DMRG algorithm.
 */
class class_iDMRG : public class_algorithm_infinite {
public:
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_infinite::class_algorithm_infinite;
//    explicit class_iDMRG(std::shared_ptr<class_hdf5_file> hdf5_);
    explicit class_iDMRG(std::shared_ptr<h5pp::File> h5ppFile_);

    void run_simulation()                                 override;
//    void run_simulation()                               override;
//    void run_preprocessing()                            override;
//    void run_postprocessing()                           override;

    void single_DMRG_step(std::string ritz);

    void check_convergence()                             override;
    bool    sim_on ()                                    override;
    long    chi_max()                                    override;
    size_t  num_sites()                                  override;
    size_t  store_freq()                                 override;
    size_t  print_freq()                                 override;
    bool    chi_grow()                                   override;
};


#endif //DMRG_CLASS_INFINITE_DMRG_H
