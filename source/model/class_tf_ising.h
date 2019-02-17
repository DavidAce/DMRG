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
    int    num_params          = 7;  //Number of parameters for this model excluding this one.

public:

    class_tf_ising();
    void set_hamiltonian(const Eigen::Tensor<Scalar,4> MPO, std::vector<double> parameters)  override;
    void set_hamiltonian(const std::vector<double> parameters)                               override;
    void set_hamiltonian(const Eigen::MatrixXd all_parameters, int position)                 override;
    void set_hamiltonian(const Eigen::VectorXd parameters)                                   override;
    void build_mpo()                                                                         override;
    void randomize_hamiltonian()                                                             override;
    Eigen::Tensor<Scalar,4> MPO_reduced_view()                                         const override;
    Eigen::Tensor<Scalar,4> MPO_reduced_view(double site_energy)                       const override;
    Eigen::MatrixXcd single_site_hamiltonian(
            int position,
            int sites,
            std::vector<Eigen::MatrixXcd> &SX,
            std::vector<Eigen::MatrixXcd> &SY,
            std::vector<Eigen::MatrixXcd> &SZ)                                          const override;
    std::unique_ptr<class_hamiltonian_base> clone()                                     const override;
    void   set_reduced_energy(double site_energy)                                             override;
    size_t get_spin_dimension()                                                         const override;
    void   print_parameter_names ()                                                     const override;
    void   print_parameter_values()                                                     const override;
    std::vector<std::string> get_parameter_names()                                      const override;
    std::vector<double>      get_parameter_values()                                     const override;
    void   set_non_local_parameters(const std::vector<std::vector<double>> chain_parameters)  override;
};

#endif //CLASS_TF_ISING_H
