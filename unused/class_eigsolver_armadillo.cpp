//
// Created by david on 2018-08-06.
//

#include "class_eigsolver_armadillo.h"
#include <armadillo>


void class_eigsolver_armadillo::subtract_phase(std::complex<double>  *eigvecs, int L) {
        using namespace std::complex_literals;
        auto uL = static_cast<unsigned int> (L);
        arma::cx_mat arma_eigvecs(eigvecs,uL,uL,false,false);
        for (int i = 0; i < L; i++) {
            auto begin = arma_eigvecs.begin() + i * L;
            auto end = begin + L;
            std::complex<double> inv_phase = -1.0i * std::arg(arma_eigvecs[i * L]);
            std::complex<double> exp_inv_phase = std::exp(inv_phase);
            std::transform(begin, end, begin,
                           [exp_inv_phase](std::complex<double> num) -> std::complex<double>
                           {
                               auto result =  num * exp_inv_phase;
                               if (std::abs(result.imag()) < 1e-15){result = result.real();}
                               return result;
                           });
        }
}


void class_eigsolver_armadillo::eig_sym(int L, std::complex<double> *hermitian_matrix, double *eigvals, std::complex<double> *eigvecs, bool remove_phase){
    // for matrices with complex elements
    auto uL = static_cast<unsigned int> (L);
    arma::cx_mat arma_matrix(hermitian_matrix, uL, uL ,false, false);
    arma::vec    arma_eigvals(eigvals,uL,false,false);
    arma::cx_mat arma_eigvecs(eigvecs,uL,uL,false,false);
    arma::eig_sym(arma_eigvals, arma_eigvecs, arma_matrix, "dc");
    if (remove_phase){
        subtract_phase(eigvecs, L);
    }
}

void class_eigsolver_armadillo::eig_sym(int L, double *symmetric_matrix, double *eigvals, double *eigvecs, bool remove_phase){
    // for matrices with complex elements
    auto uL = static_cast<unsigned int> (L);
    arma::mat arma_matrix(symmetric_matrix, uL, uL ,false, false);
    arma::vec    arma_eigvals(eigvals,uL,false,false);
    arma::mat arma_eigvecs(eigvecs,uL,uL,false,false);
    arma::eig_sym(arma_eigvals, arma_eigvecs, arma_matrix, "dc");
}