//
// Created by david on 2019-06-24.
//

#ifndef DMRG_CLASS_ALGORITHM_INFINITE_H
#define DMRG_CLASS_ALGORITHM_INFINITE_H
#include <algorithms/class_algorithm_base.h>

class class_algorithm_infinite: public class_algorithm_base {
    //Inherit the constructor of class_algorithm_base
    using class_algorithm_infinite::class_algclass_algorithm_infiniteorithm_finite;
    explicit class_iDMRG(std::shared_ptr<h5pp::File> h5ppFile_);

    void run()                                          override;
    void run_simulation()                               override;
    void run_preprocessing()                            override;
    void run_postprocessing()                           override;


};


#endif //DMRG_CLASS_ALGORITHM_INFINITE_H
