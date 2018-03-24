//
// Created by david on 7/30/17.
//

#ifndef DMRG_CLASS_ALGORITHMS_H
#define DMRG_CLASS_ALGORITHMS_H

#include <memory>
#include <IO/class_custom_cout.h>

class class_hdf5_file;

class class_algorithm_launcher  {
public:

    std::shared_ptr <class_hdf5_file> hdf5;
    class_custom_cout       ccout;
    class_algorithm_launcher(std::shared_ptr<class_hdf5_file> hdf5_);
    class_algorithm_launcher();

    void run_algorithms(){
        run_infinite_DMRG();
        run_finite_DMRG();
        run_excited_state_DMRG();
        run_imaginary_TEBD();
        run_FES_iDMRG();
        run_FES_iTEBD();
        std::cout << "All algorithms finished." << std::endl;
    };

    void run_infinite_DMRG();
    void run_finite_DMRG();
    void run_excited_state_DMRG();
    void run_imaginary_TEBD();
    void run_FES_iDMRG();
    void run_FES_iTEBD();
};


#endif //DMRG_CLASS_ALGORITHMS_H
