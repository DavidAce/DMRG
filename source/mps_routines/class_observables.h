//
// Created by david on 2017-11-12.
//

#ifndef DMRG_CLASS_OBSERVABLES_H
#define DMRG_CLASS_OBSERVABLES_H
#include <map>
#include "class_superblock.h"

enum class SimulationType {iDMRG,fDMRG, FES_iDMRG, iTEBD, FES_iTEBD};
class class_observables {
public:
    using Scalar = class_superblock::Scalar;
private:
    class_superblock &superblock;
    SimulationType sim;
    using MapType = std::map<SimulationType, std::string>;
    MapType Sim2String;
public:

    class_observables(class_superblock &superblockRef, SimulationType sim_);

    double get_expectationvalue(const Textra::Tensor<4,Scalar> &MPO);
    double get_energy();                /*! Computes the current energy by contracting the current MPS with the Hamiltonian MPO.*/
    double get_entropy();               /*! Computes the current entropy \f$ S = - \sum_n \lambda_n log( \lambda_n) \f$, where \f$\lambda_n \f$ are elements of \f$ \Lambda^A\f$ */
    double get_variance();              /*! Computes the current variance. A low value tells you that you are close to an eigenstate of the Hamiltonian. */
    double get_truncation_error();
    long   get_chi_max();
    long   get_chi();
    long   get_chain_length();

    void print_status_full(int verbosity);   /*!< Print out status of all observables.*/
    void print_status_update(int step = 0);  /*!< Print out status of all observables.*/


};


#endif //DMRG_CLASS_OBSERVABLES_H
