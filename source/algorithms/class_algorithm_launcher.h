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
        run_iDMRG();
        run_fDMRG();
        run_xDMRG();
        run_iTEBD();
        std::cout << "All simulations finished." << std::endl;
    };

    void run_iDMRG();
    void run_fDMRG();
    void run_xDMRG();
    void run_iTEBD();
};


#endif //DMRG_CLASS_ALGORITHMS_H
