//
// Created by david on 2018-04-17.
//

#ifndef CLASS_HAMILTONIAN_H
#define CLASS_HAMILTONIAN_H
#include "../../../../../../../../usr/include/c++/7/memory"
#include"../libs/eigen3/include/eigen3/unsupported/Eigen/CXX11/Tensor"
#include "../source/sim_parameters/nmspc_sim_settings.h"

class class_hamiltonian {
    using Scalar = std::complex<double>;

private:
    double J = 1.0;
    double g = 1.0;
    double e = 0.0;
    double r = 0.0;
    Eigen::array<long, 4> extent4 = {1, 1, 2, 2};     //Extent of pauli matrices in a rank-4 tensor.
    Eigen::array<long, 2> extent2 = {2, 2};           //Extent of pauli matrices in a rank-2 tensor.

    public:
    class_hamiltonian(){
        set_parameters(J, g, e, r);
    };

    class_hamiltonian(double coupling_J, double field_g, double energy_e, double field_r = 0.0){
        set_parameters(coupling_J, field_g, energy_e, field_r);
    };


    Eigen::Tensor<Scalar,4> MPO;
    Eigen::Tensor<Scalar,4> MPO_reduced();
    Eigen::Tensor<Scalar,4> MPO_reduced(double energy_e_temporary);

    void build_mpo();
    void set_parameters(double coupling_J, double field_g, double energy_e, double field_r = 0.0);
    void update_site_coupling(double coupling_J);
    void update_site_field(double field_r);
    void update_site_energy(double energy_e);
    void update_site_random_field(double field_r);

    void set_random_field(double randomness_strength);

    double get_site_coupling() const;
    double get_site_field() const;
    double get_site_energy() const;
    double get_site_random_field() const;


};


#endif //CLASS_HAMILTONIAN_H
