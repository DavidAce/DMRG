//
// Created by david on 2018-08-06.
//

#ifndef CLASS_ARMADILLO_EIGSOLVER_H
#define CLASS_ARMADILLO_EIGSOLVER_H

#include <complex>

class class_eigsolver_armadillo {
public:
    void eig_sym(int L, std::complex<double> *hermitian_matrix, double *eigvals, std::complex<double> *eigvecs, bool remove_phase = false);
    void subtract_phase(std::complex<double>  *eigvecs, int L);
};


#endif //DMRG_CLASS_ARMADILLO_EIGSOLVER_H
