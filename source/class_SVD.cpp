//
// Created by david on 7/18/17.
//

#include <class_SVD.h>
#include <n_tensor_extra.h>
#include <Eigen/SVD>

using namespace Textra;
using namespace std;
std::tuple<Tensor3, Tensor1, Tensor3> class_SVD::decompose(const Eigen::TensorRef<Tensor2> theta2, const double SVDThreshold, const long d, const long chi, const long chia, const long chib){
    Eigen::BDCSVD<MatrixType> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(tensor2_to_matrix(theta2), Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi2       = std::min(SVD.rank(),chi);
    cout << "chi chi2 chia chib: " << chi  << " " << chi2 << " " << chia << " " << chib << endl;
    renorm  = SVD.singularValues().head(chi2).norm();

    // Return XYZ
    return std::make_tuple<Tensor3,Tensor1,Tensor3>( matrix_to_tensor3(SVD.matrixU().leftCols(chi2),{d,chia,chi2}),
                                                     matrix_to_tensor1(SVD.singularValues().head(chi2)) / renorm,
                                                     matrix_to_tensor3(SVD.matrixV().leftCols(chi2),{d,chib,chi2}).shuffle(array3{0,2,1}));

}

