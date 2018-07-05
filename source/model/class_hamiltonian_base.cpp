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

class_hamiltonian_base::class_hamiltonian_base() {
    extent4 = {1, 1, spin_dim, spin_dim};
    extent2 = {spin_dim, spin_dim};
    random_field = rn::uniform_double(-randomness_strength,randomness_strength);
};


void class_hamiltonian_base::set_position(int new_pos) {
    position = new_pos;
}

int class_hamiltonian_base::get_position() const{
    return position;
}


void class_hamiltonian_base::set_site_reduced_energy(double single_site_energy) {
    energy_reduced = single_site_energy;
}


double class_hamiltonian_base::get_site_energy()const {return energy_reduced;}
double class_hamiltonian_base::get_site_random_field() const {return random_field;}
double class_hamiltonian_base::get_site_randomness_strength() const {return randomness_strength;}


