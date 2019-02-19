//
// Created by david on 2019-02-20.
//

#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/class_superblock.h>
#include <model/class_hamiltonian_base.h>

void MPS_Tools::Infinite::Print::print_hamiltonians(const class_superblock & superblock){
    superblock.HA->print_parameter_names();
    superblock.HA->print_parameter_values();
    superblock.HB->print_parameter_values();

}
