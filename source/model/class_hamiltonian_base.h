//
// Created by david on 2018-07-04.
//

#ifndef CLASS_HAMILTONIAN_BASE_H
#define CLASS_HAMILTONIAN_BASE_H

#include <memory>
#include <unsupported/Eigen/CXX11/Tensor>
#include <sim_parameters/nmspc_sim_settings.h>

class class_hamiltonian_base{
    using Scalar = std::complex<double>;
protected:
    Eigen::array<long, 4> extent4;                           /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2;                           /*!< Extent of pauli matrices in a rank-2 tensor */
    int position = 0;                                        /*!< Position on a finite chain */

public:

    class_hamiltonian_base() = default;
    virtual ~class_hamiltonian_base() = default;
    virtual std::unique_ptr<class_hamiltonian_base> clone()                     const = 0;
    Eigen::Tensor<Scalar,4> MPO;


    virtual Eigen::Tensor<Scalar,4> MPO_reduced_view()                          const = 0;
    virtual Eigen::Tensor<Scalar,4> MPO_reduced_view(double single_site_energy) const = 0;

    Eigen::MatrixXcd        MPO_matrix_view();    /*!< Matrix representation of full 2-site Hamiltonian */
    virtual Eigen::MatrixXcd single_site_hamiltonian(int position, int sites, std::vector<Eigen::MatrixXcd> &SX, std::vector<Eigen::MatrixXcd> &SY, std::vector<Eigen::MatrixXcd> &SZ) const = 0;
    virtual void   build_mpo()                                     = 0;
    virtual void   randomize_hamiltonian()                         = 0;
    virtual void   print_parameter_names ()                  const = 0;
    virtual void   print_parameter_values()                  const = 0;
    virtual void   set_reduced_energy(double site_energy)          = 0;
    virtual size_t get_spin_dimension()                      const = 0;
    virtual std::vector<std::string> get_parameter_names()   const = 0;
    virtual std::vector<double>      get_parameter_values()  const = 0;
//    virtual void   write_to_hdf5_table()                           = 0;
//    virtual double get_energy_reduced()                      const = 0;
//    virtual double get_random_field()                        const = 0;
//    virtual double get_randomness_strength()                 const = 0;

    void set_position(int new_pos);

    size_t    get_position()                 const;

};



#endif //CLASS_HAMILTONIAN_BASE_H
