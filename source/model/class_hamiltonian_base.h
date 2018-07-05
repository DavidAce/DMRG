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
    double energy_reduced      = 0;                          /*!< Energy offset for this mpo (to make "reduced" MPO views) */
    double random_field        = 0;                          /*!< Random field value */
    double randomness_strength = settings::model::w;         /*!< Randomness strength. The random field is uniformly distributed in (-w,w) */
    int    spin_dim            = settings::model::d;         /*!< Spin dimension */
    Eigen::array<long, 4> extent4;                           /*!< Extent of pauli matrices in a rank-4 tensor */
    Eigen::array<long, 2> extent2;                           /*!< Extent of pauli matrices in a rank-2 tensor */
    int position = 0;                                        /*!< Position on a finite chain */

public:

    class_hamiltonian_base();

    virtual ~class_hamiltonian_base() = default;
    virtual std::unique_ptr<class_hamiltonian_base> clone() const = 0;
    Eigen::Tensor<Scalar,4> MPO;



    virtual Eigen::Tensor<Scalar,4> MPO_reduced_view() const = 0;
    virtual Eigen::Tensor<Scalar,4> MPO_reduced_view(double single_site_energy) const = 0;

    virtual void build_mpo()                    = 0;
    virtual void randomize_field()              = 0;
    virtual void print_parameter_names () const = 0;
    virtual void print_parameter_values() const = 0;

    void set_position(int new_pos);
    void set_site_reduced_energy(double energy_e_reduced);

    int    get_position() const;
    double get_site_energy() const;
    double get_site_random_field() const;
    double get_site_randomness_strength() const;

};



#endif //CLASS_HAMILTONIAN_BASE_H
