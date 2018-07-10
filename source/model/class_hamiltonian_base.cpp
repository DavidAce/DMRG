//
// Created by david on 2018-07-04.
//

#include "class_hamiltonian_base.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_quantum_mechanics.h>
#include <general/nmspc_random_numbers.h>
#include <iomanip>

using namespace qm;
using Scalar = std::complex<double>;




void class_hamiltonian_base::set_position(int new_pos) {
    position = new_pos;
}

int class_hamiltonian_base::get_position() const{
    return position;
}

