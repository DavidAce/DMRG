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
    int    spin_dim            = settings::model::tf_ising::d;           /*!< Spin dimension */
    double J_coupling          = settings::model::tf_ising::J;
    double g_mag_field         = settings::model::tf_ising::g;
    double w_rnd_strength      = settings::model::tf_ising::w;           /*!< Randomness strength. The random field is uniformly distributed in (-w,w) */
    double r_rnd_field         = 0;                            /*!< Random field value */
    double e_reduced           = 0;                            /*!< Energy offset for this mpo (to make "reduced" MPO views) */

public:

    class_tf_ising();
    void build_mpo()                                                                  override;
    void randomize_hamiltonian()                                                      override;
    Eigen::Tensor<Scalar,4> MPO_reduced_view()                                  const override;
    Eigen::Tensor<Scalar,4> MPO_reduced_view(double site_energy)                const override;
    Eigen::MatrixXcd single_site_hamiltonian(
            int position,
            int sites,
            std::vector<Eigen::MatrixXcd> &SX,
            std::vector<Eigen::MatrixXcd> &SY,
            std::vector<Eigen::MatrixXcd> &SZ)                                  const override;

    std::unique_ptr<class_hamiltonian_base> clone()                             const override;
    void   set_reduced_energy(double site_energy)                                     override;
    int    get_spin_dimension()                                                 const override;
//    double get_energy_reduced()                                                 const override;
//    double get_random_field()                                                   const override;
//    double get_randomness_strength()                                            const override;
    void   print_parameter_names ()                                             const override;
    void   print_parameter_values()                                             const override;
    std::vector<double> get_all_parameters()                                    const;
};

#endif //CLASS_TF_ISING_H
