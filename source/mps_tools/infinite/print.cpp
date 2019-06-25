//
// Created by david on 2019-02-20.
//

#include <mps_tools/nmspc_mps_tools.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_mps_2site.h>
#include <model/class_hamiltonian_base.h>

void mpstools::infinite::print::print_hamiltonians(const class_superblock & superblock){
    superblock.HA->print_parameter_names();
    superblock.HA->print_parameter_values();
    superblock.HB->print_parameter_values();

}

void mpstools::infinite::print::print_state(const class_superblock & superblock){
    using namespace Textra;
    auto & MPS_A  = superblock.MPS->MPS_A;
    auto & MPS_B  = superblock.MPS->MPS_A;
    auto & MPS_C  = superblock.MPS->LC;
    auto & ENV_L  = superblock.Lblock;
    auto & ENV_R  = superblock.Rblock;

    std::cout << std::setprecision(10);

    std::cout << "State length              : "    << superblock.get_length()   << std::endl;
    std::cout << "State position            : "    << superblock.get_position() << std::endl;
}