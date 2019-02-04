//
// Created by david on 2018-11-15.
//

#ifndef DMRG_CLASS_RESUME_H
#define DMRG_CLASS_RESUME_H
#include <memory>
#include <iostream>
#include <sim_parameters/nmspc_sim_settings.h>

class class_hdf5_file;
class class_superblock;
class class_finite_chain_state;


class class_resume_from_hdf5 {
private:
    std::shared_ptr<class_hdf5_file>            hdf5;
    std::shared_ptr<class_superblock>           superblock;
    std::shared_ptr<class_finite_chain_state>   state;
    std::string                                 sim_name;
    SimulationType                              sim_type;

    void load_mps_chain ();

public:

    bool resume_success = false;

    class_resume_from_hdf5(
            std::shared_ptr<class_hdf5_file>            hdf5_,
            std::shared_ptr<class_superblock>           superblock_,
            std::shared_ptr<class_finite_chain_state>   state_,
            std::string                                 sim_name_,
            SimulationType                              sim_type_
            );

};


#endif //DMRG_CLASS_RESUME_H
