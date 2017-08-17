//
// Created by david on 8/15/17.
//
#include <class_superblock.h>
#include <SymEigsSolver.h>

using namespace std;
using namespace Textra;
using namespace Eigen;


//============================================================================//
// Find smallest eigenvalue using Spectra.
//============================================================================//
void class_superblock::find_ground_state(int eigSteps, double eigThreshold){

    Spectra::SymEigsSolver<double,
            Spectra::SMALLEST_ALGE,
            class_superblock>
            eigs(this, 1, 4);

    Textra::Tensor1 theta = MPS.get_theta().reshape(shape1);
    eigs.init(theta.data()); //Throw in the initial vector v0 here
    eigs.compute(eigSteps,eigThreshold, Spectra::SMALLEST_ALGE);

    if(eigs.info() != Spectra::SUCCESSFUL){
        cout << "Eigenvalue solver failed." << '\n';
        exit(1);
    }
    ground_state =  matrix_to_tensor<2>(eigs.eigenvectors(), shape2);

}

//============================================================================//
// Do SVD decomposition, truncation and normalization of the MPS.
//============================================================================//
void class_superblock::truncate(long chi, double SVDThreshold){

    Eigen::BDCSVD<MatrixType> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(tensor2_to_matrix(ground_state), Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi2       = std::min(SVD.rank(),chi);
//    ArrayXd temp    = SVD.singularValues().head(chi2);
    truncation_error= SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi2).sum();
    double renorm   = SVD.singularValues().head(chi2).norm();

    U               = matrix_to_tensor<3>(SVD.matrixU().leftCols(chi2),{d,chia,chi2});
    MPS.LA          = matrix_to_tensor<1>(SVD.singularValues().head(chi2)/renorm, {chi2});
//    MPS.LA          = matrix_to_tensor<1>(temp/renorm, {chi2});
    V               = matrix_to_tensor<3>(SVD.matrixV().leftCols(chi2),{d,chib,chi2}).shuffle(array3{0,2,1});

}