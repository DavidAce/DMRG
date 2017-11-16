//
// Created by david on 2017-11-14.
//

#ifndef DMRG_CLASS_EIG_ARPACK_WRAPPER_H
#define DMRG_CLASS_EIG_ARPACK_WRAPPER_H
#include "n_tensor_extra.h"
#include <Eigen/Core>
#include <unsupported/Eigen/ArpackSupport>


using namespace Textra;
using namespace Eigen;







/**
 * @brief powerIteration Compute the dominant eigenvalue and its relative eigenvector of a square matrix
 * @param A The input matrix
 * @param eigenVector The eigenvector
 * @param tolerance Maximum tolerance
 * @param nIterations Number of iterations
 * @return The dominant eigenvalue
 */
//template < typename Scalar>
double powerIteration(const Textra::Tensor<2,double>& tensor, Textra::Tensor<1,double>& eigvec, double tolerance,  int nIterations)
{
    Textra::Tensor<1,double> approx(tensor.dimension(1));
    approx.setRandom();
    int counter = 0;
    double error=100.0;
    Textra::Tensor<0,double> norm;
    norm(0) = 1.0;
    while (counter < nIterations && error > tolerance  )
    {
        Textra::Tensor<1,double> temp = approx/norm(0);
        approx = tensor.contract(temp, idx<1>({1},{0}));
        norm =  approx.conjugate().contract(approx, idx<1>({0},{0})).sqrt().real();
        error = Tensor1d_to_VectorXd(temp-approx/norm(0)).stableNorm();
        counter++;
    }
    eigvec = approx;
    Textra::Tensor<0,double> dominantExtract = approx.contract(tensor,idx<1>({0},{0})).contract(approx,idx<1>({0},{0}));
    double dominantEigenValue = dominantExtract(0);
    cerr << "Power Iteration:" << endl;
    cerr << "\tTotal iterations= " << counter << endl;
    cerr << "\tError= " << error << endl;
    cerr << "\tDominant Eigenvalue= " << dominantEigenValue << endl;
    cerr << "\tDominant Eigenvector= [" << eigvec<< "]" << endl;
    return dominantEigenValue;
}
#endif //DMRG_CLASS_EIG_ARPACK_WRAPPER_H
