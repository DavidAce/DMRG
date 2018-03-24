//
// Created by david on 2017-10-04.
//

#ifndef DMRG_CLASS_SVD_H
#define DMRG_CLASS_SVD_H

#ifdef MKL_AVAILABLE
#define  EIGEN_USE_MKL_ALL
#endif

#include "nmspc_tensor_extra.h"
#include <Eigen/SVD>
template<typename Scalar>
class class_SVD{
private:
    double SVDThreshold         = 1e-12;
    double truncation_error     = 0;
    int chi                     = 0;
    using MatrixType  = Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    Eigen::BDCSVD<MatrixType> SVD;
public:

    double get_truncation_error();
    void setThreshold(double newThreshold);
    auto decompose(const Textra::Tensor<Scalar,2> &tensor);
    auto decompose(const Textra::Tensor<Scalar,2> &tensor, const long chi_max);
    auto schmidt  (const Textra::Tensor<Scalar,2> &tensor, long d, long chiL, long chi_max, long chiR);
    auto schmidt  (const Textra::Tensor<Scalar,2> &tensor, long d, long chiL,               long chiR);
    auto schmidt  (const Textra::Tensor<Scalar,4> &tensor, long d, long chiL, long chi_max,  long chiR);
    auto schmidt  (const Textra::Tensor<Scalar,4> &tensor, long d, long chiL,                long chiR);
    Textra::Tensor<Scalar, 2> pseudo_inverse(const Textra::Tensor<Scalar,2> &tensor);

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
Textra::Tensor<Scalar, 2> class_SVD<Scalar>::pseudo_inverse(const Textra::Tensor<Scalar, 2> &tensor){
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    setThreshold(1e-8);
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    long chi = SVD.rank();
    Textra::Tensor<Scalar, 2> U = Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi);
    Textra::Tensor<Scalar, 1> S = Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).template cast<Scalar>() , chi);
    Textra::Tensor<Scalar, 2> V = Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate() ,  SVD.matrixV().rows(), chi ).shuffle(Textra::array2{1,0});

    return V.conjugate().shuffle(Textra::array2{1,0})
            .contract(Textra::asDiagonalInversed(S),Textra::idx({1},{0}))
            .contract(U.conjugate().shuffle(Textra::array2{1,0}), Textra::idx({1},{0}));
}



template<typename Scalar>
auto class_SVD<Scalar>::decompose(const Textra::Tensor<Scalar,2> &tensor) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    long chi = SVD.rank();
    return std::make_tuple
            <Textra::Tensor<Scalar, 2> ,Textra::Tensor<Scalar, 1>, Textra::Tensor<Scalar, 2> >
            (
                    Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi),
                    Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>() , chi),
                    Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate() ,  chi, SVD.matrixV().rows() ).shuffle(Textra::array2{1,0})
            );
}



template<typename Scalar>
auto class_SVD<Scalar>::decompose(const Textra::Tensor<Scalar,2> &tensor, const long chi_max) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    long chi = std::min(SVD.rank(), chi_max);
    return std::make_tuple
            <Textra::Tensor<Scalar, 2> ,Textra::Tensor<Scalar, 1>, Textra::Tensor<Scalar, 2> >
            (
                Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi),
                Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>() , chi),
                Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate() ,  chi, SVD.matrixV().rows() ).shuffle(Textra::array2{1,0})
            );
}

template<typename Scalar>
auto class_SVD<Scalar>::schmidt(const Textra::Tensor<Scalar,2> &tensor, long d, long chiL, long chi_max, long chiR) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi            = std::min(SVD.rank(),chi_max);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    return std::make_tuple
            <Textra::Tensor<Scalar, 3> ,Textra::Tensor<Scalar, 1>, Textra::Tensor<Scalar, 3> >
            (
                Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), d, chiL, chi),
                Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>(),  chi ),
                Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(),  d, chiR, chi ).shuffle(Textra::array3{ 0, 2, 1 })
            );
}


template<typename Scalar>
auto class_SVD<Scalar>::schmidt(const Textra::Tensor<Scalar,2> &tensor, long d, long chiL, long chiR) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi            = SVD.rank();
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    return std::make_tuple
            <Textra::Tensor<Scalar, 3> ,Textra::Tensor<Scalar, 1>, Textra::Tensor<Scalar, 3> >
            (
                Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), d, chiL, chi),
                Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>(),  chi ),
                Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(),  d, chiR, chi ).shuffle(Textra::array3{ 0, 2, 1 })
            );
}


template<typename Scalar>
auto class_SVD<Scalar>::schmidt(const Textra::Tensor<Scalar,4> &tensor, long d, long chiL, long chiR) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi            = SVD.rank();
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    return std::make_tuple
            <Textra::Tensor<Scalar, 3> ,Textra::Tensor<Scalar, 1>, Textra::Tensor<Scalar, 3> >
            (
                Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), d, chiL, chi),
                Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>(), chi),
                Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(),  d, chiR, chi).shuffle(Textra::array3{0,2,1})
            );
}

template<typename Scalar>
auto class_SVD<Scalar>::schmidt(const Textra::Tensor<Scalar,4> &tensor, long d, long chiL, long chi_max, long chiR) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi            = std::min(SVD.rank(),chi_max);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    return std::make_tuple
            <Textra::Tensor<Scalar, 3> ,Textra::Tensor<Scalar, 1>, Textra::Tensor<Scalar, 3> >
            (
                Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), d, chiL, chi),
                Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>(),  chi ),
                Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(),  d, chiR, chi).shuffle(Textra::array3{0,2,1})
            );
}

#endif //DMRG_CLASS_SVD_H
