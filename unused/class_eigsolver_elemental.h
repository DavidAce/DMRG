//
// Created by david on 2018-06-18.
//

#ifndef CLASS_ELEMENTAL_HERMITIANEIG_H
#define CLASS_ELEMENTAL_HERMITIANEIG_H

#include <El.hpp>
enum class LOC {LOW,MID,HIGH};


class class_eigsolver_elemental {
private:
    using Scalar        = El::Complex<double>;
public:
    class_eigsolver_elemental() = default;
    El::HermitianEigCtrl<Scalar> ctrl;
    void eig_sym(int L, std::complex<double> *hermitian_matrix, double *eigvals, std::complex<double> *eigvecs);
    void eig_sym_indexSubset(LOC loc, int numVals, int L, std::complex<double> *hermitian_matrix, double *eigvals, std::complex<double> *eigvecs);
};


#endif //TRAINING_CLASS_ELEMENTAL_HERMITIANEIG_H
