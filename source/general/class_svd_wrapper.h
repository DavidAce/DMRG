//
// Created by david on 2017-10-04.
//

#ifndef DMRG_CLASS_SVD_H
#define DMRG_CLASS_SVD_H
//#define EIGEN_USE_MKL_ALL
#include "n_tensor_extra.h"
#include <Eigen/SVD>
template<typename Scalar>
class class_SVD{
private:
    double SVDThreshold         = 1e-12;
    double truncation_error     = 0;
    int chi                     = 0;

    using MatrixType  = Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    template <long rankU, long rankS, long rankV>
    using TensorTuple = std::tuple<Tensor<Scalar,rankU>, // Suzuki_Trotter_2nd_order
            Eigen::Tensor<Scalar,rankS>, // S
            Eigen::Tensor<Scalar,rankV>>; // V^T
    Eigen::BDCSVD<MatrixType> SVD;
public:

    double get_truncation_error();
    void setThreshold(double newThreshold);
    TensorTuple<2,1,2> decompose(const Tensor<Scalar,2> &tensor);
    TensorTuple<2,1,2> decompose(const Tensor<Scalar,2> &tensor, const long chi_max);
    TensorTuple<3,1,3> schmidt  (const Tensor<Scalar,2> &tensor, long d, long chiL, long chi_max, long chiR);
};

//
// Definitions
//

template<typename Scalar>
double class_SVD<Scalar>::get_truncation_error(){
    return truncation_error;
}

template<typename Scalar>
void class_SVD<Scalar>::setThreshold(double newThreshold) {
    SVD.setThreshold(newThreshold);
}


template<typename Scalar>
typename class_SVD<Scalar>:: template TensorTuple<2,1,2> class_SVD<Scalar>::decompose(const Tensor<Scalar,2> &tensor) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    return std::make_tuple
            (
                    Matrix_to_Tensor<Scalar,2>(SVD.matrixU(), {SVD.matrixU().rows(), SVD.matrixU().cols() }),
                    Matrix_to_Tensor<Scalar,1>(SVD.singularValues().normalized(), {SVD.singularValues().size()}),
                    Matrix_to_Tensor<Scalar,2>(SVD.matrixV().transpose(),  {SVD.matrixV().cols(), SVD.matrixV().rows() })
            );
}



template<typename Scalar>
typename class_SVD<Scalar>:: template TensorTuple<2,1,2> class_SVD<Scalar>::decompose(const Tensor<Scalar,2> &tensor, const long chi_max) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    long chi = std::min(SVD.rank(), chi_max);
    return std::make_tuple
            (
                    Matrix_to_Tensor<Scalar,2>(SVD.matrixU().leftCols(chi), {SVD.matrixU().rows(), chi}),
                    Matrix_to_Tensor<Scalar,1>(SVD.singularValues().head(chi).normalized() , {chi}),
                    Matrix_to_Tensor<Scalar,2>(SVD.matrixV().leftCols(chi).transpose() ,  {chi, SVD.matrixV().rows() })
            );
}

template<typename Scalar>
typename class_SVD<Scalar>:: template TensorTuple<3,1,3> class_SVD<Scalar>::schmidt(const Tensor<Scalar,2> &tensor, long d, long chiL, long chi_max, long chiR) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi            = std::min(SVD.rank(),chi_max);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    return std::make_tuple
            (
                    Matrix_to_Tensor<Scalar,3>(SVD.matrixU().leftCols(chi), {d, chiL, chi}),
                    Matrix_to_Tensor<Scalar,1>(SVD.singularValues().head(chi).normalized(), { chi }),
                    Matrix_to_Tensor<Scalar,3>(SVD.matrixV().leftCols(chi), { d, chiR, chi }).shuffle(array3{ 0, 2, 1 })
            );
}
#endif //DMRG_CLASS_SVD_H
