//
// Created by david on 2017-11-12.
//

#ifndef DMRG_CLASS_OBSERVABLES_H
#define DMRG_CLASS_OBSERVABLES_H
#include <map>
#include <memory>
#include <IO/class_custom_cout.h>
#include <unsupported/Eigen/CXX11/Tensor>
#include <sim_parameters/nmspc_sim_settings.h>
class class_superblock;

/*!
 * \class class_measurement
 * \brief A class for measuring observables
 * This class extracts observables, expectation values, like energy and entropy from the MPS, as well as other useful numbers like \f$\chi\f$ and truncation errors.
*/

class class_measurement {
public:
    using Scalar = double;
private:
    std::shared_ptr<class_superblock> superblock;
    class_custom_cout ccout;
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic>              lambdaG;
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> eigvecG;
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic>              lambdaID;
    Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> eigvecID;
public:
    SimulationType sim;

    double first_moment();
    explicit class_measurement(std::shared_ptr<class_superblock> superblock_, SimulationType sim_);

    double variance1 = 1; double variance2 = 1; double variance3 = 1;
    double get_expectationvalue(const Eigen::Tensor<double,4> &MPO);
    double get_expectationvalue(const Eigen::Tensor<std::complex<double>,4> &MPO);
    double get_expectationvalue_sq(const Eigen::Tensor<double,4> &MPO);
    double get_expectationvalue_sq(const Eigen::Tensor<std::complex<double>,4> &MPO);
    double get_energy();                 /*! Computes the current energy by contracting the current MPS with the Hamiltonian MPO.*/
    double get_entanglement_entropy();                /*! Computes the current entropy \f$ S = - \sum_n \lambda_n log( \lambda_n) \f$, where \f$\lambda_n \f$ are elements of \f$ \Lambda^A\f$ */
    double get_full_variance();          /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
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
