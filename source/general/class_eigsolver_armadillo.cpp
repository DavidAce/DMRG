//
// Created by david on 2018-08-06.
//

#include "class_eigsolver_armadillo.h"
#include <armadillo>

void class_eigsolver_armadillo::eig_sym(int L, std::complex<double> *hermitian_matrix, double *eigvals, std::complex<double> *eigvecs){
    // for matrices with complex elements
    auto uL = static_cast<unsigned int> (L);
    arma::cx_mat arma_matrix(hermitian_matrix, uL, uL ,false, false);
    arma::vec    arma_eigvals(eigvals,uL,false,false);
    arma::cx_mat arma_eigvecs(eigvecs,uL,uL,false,false);
    arma::eig_sym(arma_eigvals, arma_eigvecs, arma_matrix, "dc");
}