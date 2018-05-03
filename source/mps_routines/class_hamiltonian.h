//
// Created by david on 2018-04-17.
//

#ifndef CLASS_HAMILTONIAN_H
#define CLASS_HAMILTONIAN_H

#include<unsupported/Eigen/CXX11/Tensor>

class class_hamiltonian {
    using Scalar = std::complex<double>;

private:
    double J = 1.0;
    double g = 1.0;
    double e = 0.0;
    Eigen::array<long, 4> extent4 = {1, 1, 2, 2};     //Extent of pauli matrices in a rank-4 tensor.
    Eigen::array<long, 2> extent2 = {2, 2};           //Extent of pauli matrices in a rank-2 tensor.
public:
    class_hamiltonian(){
        set_parameters(J, g, e);
    };

    class_hamiltonian(double coupling_J, double field_g, double energy_e){
        set_parameters(coupling_J, field_g, energy_e);
    };

    Eigen::Tensor<Scalar,4> MPO;
    void build_mpo();
    void set_parameters(double coupling_J, double field_g, double energy_e);
    void update_site_coupling(double coupling_J);
    void update_site_energy(double energy_e);
    void update_site_field(double field_g);

    Eigen::Tensor<Scalar,4> MPO_zero_site_energy();

};


#endif //CLASS_HAMILTONIAN_H
