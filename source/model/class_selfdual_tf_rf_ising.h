//
// Created by david on 2018-07-06.
//

#ifndef CLASS_SELFDUAL_TF_ISING_H
#define CLASS_SELFDUAL_TF_ISING_H

#include <iostream>
#include <general/nmspc_tensor_extra.h>
#include <iomanip>
#include "class_hamiltonian_base.h"


class class_selfdual_tf_rf_ising : public class_hamiltonian_base {
    using Scalar = std::complex<double>;
private:
    int    spin_dim            = settings::model::selfdual_tf_rf_ising::d;           /*!< Spin dimension */
    double J_rnd               = 0;
    double h_rnd               = 0;
    double J_avg               = settings::model::selfdual_tf_rf_ising::J_avg;
    double h_avg               = settings::model::selfdual_tf_rf_ising::h_avg;
    double J_std               = settings::model::selfdual_tf_rf_ising::J_std;
    double h_std               = settings::model::selfdual_tf_rf_ising::h_std;
    double lambda              = settings::model::selfdual_tf_rf_ising::lambda;
    double e_reduced           = 0;                            /*!< Energy offset for this mpo (to make "reduced" MPO views) */

public:

    class_selfdual_tf_rf_ising();
    void build_mpo()                                                            override;
    void randomize_hamiltonian()                                                override;
    Eigen::Tensor<Scalar,4> MPO_reduced_view()                            const override;
    Eigen::Tensor<Scalar,4> MPO_reduced_view(double site_energy)          const override;
    std::unique_ptr<class_hamiltonian_base> clone()                       const override;
    void   set_reduced_energy(double site_energy)                               override;
    int    get_spin_dimension()                                           const override;
//    double get_energy_reduced()                                           const override;
//    double get_random_field()                                             const override;
//    double get_randomness_strength()                                      const override;
    void   print_parameter_names ()                                       const override;
    void   print_parameter_values()                                       const override;
    std::vector<double> get_all_parameters()                              const;
};

#endif //DMRG_CLASS_SELFDUAL_TF_ISING_H
