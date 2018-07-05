//
// Created by david on 2018-07-04.
//

#ifndef CLASS_TF_ISING_H
#define CLASS_TF_ISING_H

#include <iostream>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include "class_hamiltonian_base.h"

class class_tf_ising : public class_hamiltonian_base {
    using Scalar = std::complex<double>;
private:
    double J_coupling = settings::model::tf_ising::J;
    double g_field    = settings::model::tf_ising::g;

public:

    class_tf_ising(): class_hamiltonian_base(){
        build_mpo();
        std::cout << "random field  : " << random_field   << std::endl;
        std::cout << "random stren  : " << randomness_strength   << std::endl;
        std::cout << "J             : " << J_coupling     << std::endl;
        std::cout << "g             : " << g_field        << std::endl;
        std::cout << "energy_reduce : " << energy_reduced << std::endl;
        std::cout << "MPO: \n"  << std::setprecision(4) << MPO.reshape(Textra::array2{6,6}) << std::endl;
    };


    void build_mpo() override;
    void randomize_field() override;

    Eigen::Tensor<Scalar,4> MPO_reduced_view() const override;
    Eigen::Tensor<Scalar,4> MPO_reduced_view(double single_site_energy) const override;
    std::unique_ptr<class_hamiltonian_base> clone() const override;

    void print_parameter_names () const override;
    void print_parameter_values() const override;

    std::vector<double> get_all_parameters() const;


};




#endif //CLASS_TF_ISING_H
