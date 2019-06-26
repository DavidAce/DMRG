//
// Created by david on 2018-07-06.
//

#ifndef CLASS_SELFDUAL_TF_ISING_H
#define CLASS_SELFDUAL_TF_ISING_H

#include <iostream>
#include <general/nmspc_tensor_extra.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <iomanip>
#include "class_hamiltonian_base.h"

template<typename table_type> class class_hdf5_table;
class class_selfdual_tf_rf_ising_table;

class class_selfdual_tf_rf_ising : public class_hamiltonian_base {
    using Scalar = std::complex<double>;
private:
    int    spin_dim            = 0;           /*!< Spin dimension */
    double J_rnd               = 0;
    double h_rnd               = 0;
    double J_log_mean          = 0;
    double h_log_mean          = 0;
    double J_avg               = 0;
    double h_avg               = 0;
    double J_sigma             = 0;
    double h_sigma             = 0;
    double lambda              = 0;
    double delta               = 0;
    double e_reduced           = 0;                            /*!< Energy offset for this mpo (to make "reduced" MPO views) */

    int    num_params = 13;  //Number of parameters for this model excluding this one.


public:

    class_selfdual_tf_rf_ising(size_t position_,std::string logName = "SDUAL-ISING");

    // Functions that extend the base (no override)
    void set_realization_averages(double J_avg_,double h_avg_);
    auto get_J_rnd(){return J_rnd;};
    auto get_h_rnd(){return h_rnd;};


    // Functions that override the base
    void set_hamiltonian(const Eigen::Tensor<Scalar,4> MPO_, std::vector<double> parameters) override;
    void set_hamiltonian(const std::vector<double> parameters)                                      override;
    void set_hamiltonian(const Eigen::MatrixXd all_parameters, int position)                  override;
    void set_hamiltonian(const Eigen::VectorXd parameters)                                    override;


    void build_mpo()                                                                          override;
    void randomize_hamiltonian()                                                              override;
    Eigen::Tensor<Scalar,4> MPO_reduced_view()                                          const override;
    Eigen::Tensor<Scalar,4> MPO_reduced_view(double site_energy)                        const override;
    Eigen::MatrixXcd single_site_hamiltonian(
            int position,
            int sites,
            std::vector<Eigen::MatrixXcd> &SX,
            std::vector<Eigen::MatrixXcd> &SY,
            std::vector<Eigen::MatrixXcd> &SZ)                                          const override;
    std::unique_ptr<class_hamiltonian_base> clone()                                     const override;
    void   set_reduced_energy(double site_energy)                                             override;
    size_t get_spin_dimension()                                                         const override;

//    double get_energy_reduced()                                                        const override;
//    double get_random_field()                                                          const override;
//    double get_randomness_strength()                                                   const override;
    void   print_parameter_names ()                                                     const override;
    void   print_parameter_values()                                                     const override;
    std::vector<std::string> get_parameter_names()                                      const override;
    std::vector<double>      get_parameter_values()                                     const override;

    void   set_full_lattice_parameters(const std::vector<std::vector<double>> chain_parameters)  override;


//    void   write_to_hdf5_table()                                                override;




};

#endif //DMRG_CLASS_SELFDUAL_TF_ISING_H
