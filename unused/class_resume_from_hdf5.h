//
// Created by david on 2018-11-15.
//

#ifndef DMRG_CLASS_RESUME_H
#define DMRG_CLASS_RESUME_H
#include <memory>
#include <iostream>
#include <sim_parameters/nmspc_sim_settings.h>

class class_superblock;
class class_finite_chain_state;


class class_resume_from_hdf5 {
private:

public:

    static bool load_success;
    static void load_state_from_hdf5 (
            const class_hdf5_file & hdf5,
                  class_finite_chain_state & state,
                  class_superblock &superblock,
                  std::string sim_name,
                  SimulationType sim_type
            );


};







#endif //DMRG_CLASS_RESUME_H
