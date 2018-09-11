//
// Created by david on 2018-06-18.
//

#include "class_eigsolver_elemental.h"


void class_eigsolver_elemental::eig_sym(int L, std::complex<double> *hermitian_matrix, double *eigvals,
                                        std::complex<double> *eigvecs){
    El::Matrix<Scalar> matrix(L,L, static_cast<Scalar*>(hermitian_matrix),L);
    El::Matrix<Scalar> vecs;
    El::Matrix<double> vals;
    El::HermitianEig( El::LOWER, matrix, vals, vecs,ctrl );
    std::move(vecs.Buffer(), vecs.Buffer() + vecs.Height()*vecs.Width(), eigvecs);
    std::move(vals.Buffer(), vals.Buffer() + vals.Height()*vals.Width(), eigvals);
}

  void class_eigsolver_elemental::eig_sym_indexSubset(LOC loc, int numVals, int L,
                                                      std::complex<double> *hermitian_matrix,
                                                      double *eigvals,
                                                      std::complex<double> *eigvecs){

        int lowIdx;
        int highIdx;
        ctrl.tridiagEigCtrl.subset.indexSubset = true;
        switch (loc){
            case LOC::LOW:
                lowIdx = 0;
                highIdx = numVals - 1;
                break;
            case LOC::MID:
                lowIdx = (L - numVals)/2;
                highIdx= lowIdx + numVals - 1;
                break;
            case LOC::HIGH:
                highIdx = L - 1;
                lowIdx = highIdx - numVals + 1;
                break;
        }
        std::cout << "width: " << highIdx - lowIdx << std::endl;
        ctrl.tridiagEigCtrl.subset.lowerIndex = lowIdx;
        ctrl.tridiagEigCtrl.subset.upperIndex = highIdx;

        El::Matrix<Scalar> matrix(L,L, static_cast<Scalar*>(hermitian_matrix),L);
        El::Matrix<Scalar> vecs;
        El::Matrix<double> vals;
        El::HermitianEig         (El::LOWER, matrix, vals, vecs, ctrl);
        std::move(vecs.Buffer(), vecs.Buffer() + vecs.Height()*vecs.Width(), eigvecs);
        std::move(vals.Buffer(), vals.Buffer() + vals.Height()*vals.Width(), eigvals);
    }
