//
// Created by david on 2018-05-09.
//

#include "class_svd_wrapper.h"
#include <iomanip>
#include <Eigen/QR>
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
Eigen::Tensor<Scalar, 2>
class_SVD<Scalar>::pseudo_inverse(const Eigen::Tensor<Scalar, 2> &tensor){
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    return Textra::Matrix_to_Tensor2(mat.completeOrthogonalDecomposition().pseudoInverse() );
//    setThreshold(1e-8);
//    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
//    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).squaredNorm();
//    long chi = SVD.rank();
//    Eigen::Tensor<Scalar, 2> U = Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi);
//    Eigen::Tensor<Scalar, 1> S = Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).template cast<Scalar>() , chi);
//    Eigen::Tensor<Scalar, 2> V = Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate() ,  SVD.matrixV().rows(), chi ).shuffle(Textra::array2{1,0});
//
//
//    return V.conjugate().shuffle(Textra::array2{1,0})
//            .contract(Textra::asDiagonalInversed(S),Textra::idx({1},{0}))
//            .contract(U.conjugate().shuffle(Textra::array2{1,0}), Textra::idx({1},{0}));
}



template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
class_SVD<Scalar>::decompose(const Eigen::Tensor<Scalar,2> &tensor) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).squaredNorm();
    long chi = SVD.rank();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>() , chi),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(), SVD.matrixV().rows(), chi).shuffle(Textra::array2{1,0})
                            );
}

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
class_SVD<Scalar>::decompose(const Eigen::Tensor<Scalar,3> &tensor,const long rows,const long cols) {
    Eigen::Map<const MatrixType> mat (tensor.data(), rows, cols);
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi = SVD.rank();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>() , chi),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(), SVD.matrixV().rows(), chi).shuffle(Textra::array2{1,0})
    );
}



template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
class_SVD<Scalar>::decompose(const Eigen::Tensor<Scalar,2> &tensor, const long chi_max) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).squaredNorm();
    long chi = std::min(SVD.rank(), chi_max);
    return std::make_tuple
                           (Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi),
                            Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>() , chi),
                            Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate() ,  chi, SVD.matrixV().rows() ).shuffle(Textra::array2{1,0})
                           );
}

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,2> &tensor, long d, long chiL, long chi_max, long chiR) {
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi            = std::min(SVD.rank(),chi_max);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), d, chiL, chi),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>(),  chi ),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(),  d, chiR, chi ).shuffle(Textra::array3{ 0, 2, 1 })
                           );
}


template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,2> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC            = SVD.rank();
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(),  chiC ),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC ).shuffle(Textra::array3{ 0, 2, 1 })
                            );
}


template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,3> &tensor, long rows,long cols) {
    long d    = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long chiR = tensor.dimension(2);
    Eigen::Map<const MatrixType> mat (tensor.data(), rows,cols);
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC            = SVD.rank();
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), d, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(),  chiC ),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  d, chiR, chiC ).shuffle(Textra::array3{ 0, 2, 1 })
    );
}



template<typename Scalar>
std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,4> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC            = SVD.rank();
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(), chiC),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC).shuffle(Textra::array3{0,2,1})
                           );
}

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,4> &tensor, long chi_max) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC            = std::min(SVD.rank(),chi_max);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(
                           Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(),  chiC ),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC).shuffle(Textra::array3{0,2,1})
                          );
}


template<typename Scalar>
std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>,double  >
class_SVD<Scalar>::schmidt_with_norm(const Eigen::Tensor<Scalar,4> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC            = SVD.rank();


    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(), chiC),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC).shuffle(Textra::array3{0,2,1}),
                           SVD.singularValues().head(chiC).norm()
    );
}

template<typename Scalar>
std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>  >
class_SVD<Scalar>::schmidt_unnormalized(const Eigen::Tensor<Scalar,4> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC            = SVD.rank();


    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).template cast<Scalar>(), chiC),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC).shuffle(Textra::array3{0,2,1})
    );
}
// ============================ //
//   Explicit instantiations
// ============================ //

template class class_SVD<double>;
template class class_SVD<std::complex<double>>;

