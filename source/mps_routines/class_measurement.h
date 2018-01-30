//
// Created by david on 2017-11-12.
//

#ifndef DMRG_CLASS_OBSERVABLES_H
#define DMRG_CLASS_OBSERVABLES_H
#include <map>
#include <memory>
#include <IO/class_custom_cout.h>
#include <general/n_tensor_extra.h>
#include <mps_routines/class_mps.h>
#include <sim_parameters/n_sim_settings.h>
#include <mps_routines/algorithms/class_base.h>

class class_superblock;
class class_measurement {
public:
    using Scalar = class_mps::Scalar;
private:
    std::shared_ptr<class_superblock> superblock;
    class_custom_cout ccout;

public:
    SimulationType sim;

    double first_moment();
    explicit class_measurement(std::shared_ptr<class_superblock> superblock_, SimulationType sim_);

    double variance1,variance2,variance3;
    double get_expectationvalue(const Textra::Tensor<double,4> &MPO);
    double get_expectationvalue(const Textra::Tensor<std::complex<double>,4> &MPO);
    double get_expectationvalue_sq(const Textra::Tensor<double,4> &MPO);
    double get_expectationvalue_sq(const Textra::Tensor<std::complex<double>,4> &MPO);
    double get_energy();                 /*! Computes the current energy by contracting the current MPS with the Hamiltonian MPO.*/
    double get_entropy();                /*! Computes the current entropy \f$ S = - \sum_n \lambda_n log( \lambda_n) \f$, where \f$\lambda_n \f$ are elements of \f$ \Lambda^A\f$ */
    double get_variance();               /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance1();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance2();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_variance3();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_truncation_error();
    double get_second_cumulant();
    long   get_chi();
    long   get_chain_length();
};


#endif //DMRG_CLASS_OBSERVABLES_H
