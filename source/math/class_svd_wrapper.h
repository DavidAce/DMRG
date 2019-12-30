//
// Created by david on 2017-10-04.
//

#pragma once

#include "general/nmspc_tensor_extra.h"
#include <iomanip>


class class_SVD{
private:
    double SVDThreshold         = 1e-12;
    double truncation_error     = 0;
    template <typename Scalar> using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
    template <typename Scalar> using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;

    template<typename Scalar>
    std::tuple<MatrixType<Scalar>, VectorType<Scalar>,MatrixType<Scalar>, long>
    do_svd_lapacke(const Scalar * mat_ptr, long rows, long cols, long rank_max);


    template<typename Scalar>
    std::tuple<MatrixType<Scalar>, VectorType<Scalar>,MatrixType<Scalar>, long>
    do_svd(const Scalar * mat_ptr, long rows, long cols, long rank_max);

    template<typename Derived>
    std::tuple<MatrixType<typename Derived::Scalar>, VectorType<typename Derived::Scalar>,MatrixType<typename Derived::Scalar>, long>
    do_svd(const Eigen::MatrixBase<Derived> & mat, long rank_max );

    template<typename Derived>
    std::tuple<MatrixType<typename Derived::Scalar>, VectorType<typename Derived::Scalar>,MatrixType<typename Derived::Scalar>, long>
    do_svd(const Eigen::MatrixBase<Derived> & mat );




public:

    class_SVD()=default;
    bool   use_lapacke          = false;
    double get_truncation_error();
    void setThreshold(double newThreshold);

    template<typename Scalar>
    Eigen::Tensor<Scalar, 2> pseudo_inverse(const Eigen::Tensor<Scalar,2> &tensor);

    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,2> &tensor);
    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,2> &tensor, const long chi_max);

    template<typename Scalar>
    std::tuple<MatrixType<Scalar> ,VectorType<Scalar>, MatrixType<Scalar>>
    decompose(const MatrixType<Scalar> &matrix);

    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,3> &tensor, const long rows,const long cols);

    template<typename Scalar,auto rank>
    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,rank> &tensor, long dL, long dR, long chiL, long chi_max, long chiR);

    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,3> &tensor, long rows, long cols);

    template<typename Scalar>
    std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,4> &tensor, long chi_max);
    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,4> &tensor);
    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>, double >
    schmidt_with_norm  (const Eigen::Tensor<Scalar,4> &tensor);
    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>, double >
    schmidt_with_norm  (const Eigen::Tensor<Scalar,4> &tensor, long chi_max);

};


//
// Definitions
//

template<typename Derived>
std::tuple<class_SVD::MatrixType<typename Derived::Scalar>, class_SVD::VectorType<typename Derived::Scalar>,class_SVD::MatrixType<typename Derived::Scalar>, long>
class_SVD::do_svd(const Eigen::MatrixBase<Derived> & mat, long rank_max ){
    return do_svd(mat.derived().data(), mat.rows(),mat.cols(),rank_max);
}

template<typename Derived>
std::tuple<class_SVD::MatrixType<typename Derived::Scalar>, class_SVD::VectorType<typename Derived::Scalar>,class_SVD::MatrixType<typename Derived::Scalar>, long>
class_SVD::do_svd(const Eigen::MatrixBase<Derived> & mat ){
    long rank_max = std::min(mat.rows(),mat.cols());
    return do_svd(mat.derived().data(), mat.rows(),mat.cols(),rank_max);
}





template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
class_SVD::decompose(const Eigen::Tensor<Scalar,2> &tensor) {
    Eigen::Map<const MatrixType<Scalar>> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    auto[U,S,V,rank] = do_svd(mat);
    return std::make_tuple(Textra::MatrixTensorMap(U),
                           Textra::MatrixTensorMap(S.normalized().template cast<Scalar>()),
                           Textra::MatrixTensorMap(V)
    );
}

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
class_SVD::decompose(const Eigen::Tensor<Scalar,3> &tensor,const long rows,const long cols) {
    auto tensormap = Eigen::TensorMap<Eigen::Tensor<Scalar,2>> (tensor.data(), rows,cols);
    return decompose(tensormap);
}



template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
class_SVD::decompose(const Eigen::Tensor<Scalar,2> &tensor, const long chi_max) {
    Eigen::Map<const MatrixType<Scalar>> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    auto[U,S,V,rank] = do_svd(mat,chi_max);
    return std::make_tuple(Textra::MatrixTensorMap(U),
                           Textra::MatrixTensorMap(S.normalized().template cast<Scalar>()),
                           Textra::MatrixTensorMap(V)
    );
}

template<typename Scalar>
std::tuple<class_SVD::MatrixType<Scalar> ,class_SVD::VectorType<Scalar>, class_SVD::MatrixType<Scalar>>
class_SVD::decompose(const class_SVD::MatrixType<Scalar> &matrix){
    auto[U,S,V,rank] = do_svd(matrix);
    return std::make_tuple(U,
                           S.normalized().template cast<Scalar>(),
                           V
    );
}


template<typename Scalar, auto tensor_rank>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD::schmidt(const Eigen::Tensor<Scalar,tensor_rank> &tensor, long dL, long dR, long chiL, long chi_max, long chiR) {
    if (dL*chiL * dR*chiR != tensor.size()){throw std::range_error("schmidt error: tensor size does not match given dimensions.");}
    Eigen::Map<const MatrixType<Scalar>> mat (tensor.data(), dL*chiL, dR*chiR);
    auto [U,S,V,rank] = do_svd(mat,chi_max);
    return std::make_tuple(Textra::MatrixTensorMap(U, dL, chiL, rank),
                           Textra::MatrixTensorMap(S.normalized().template cast<Scalar>(), rank),
                           Textra::MatrixTensorMap(V,  rank, dR, chiR ).shuffle(Textra::array3{ 1, 0, 2 })
    );
}





template<typename Scalar>
std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD::schmidt(const Eigen::Tensor<Scalar,4> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    long chi_max = std::min(dL*chiL, dR*chiR);
    return schmidt(tensor,dL,dR,chiL,chi_max,chiR);

}

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD::schmidt(const Eigen::Tensor<Scalar,4> &tensor, long chi_max) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    return schmidt(tensor,dL,dR,chiL,chi_max,chiR);
}


template<typename Scalar>
std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>,double  >
class_SVD::schmidt_with_norm(const Eigen::Tensor<Scalar,4> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    long chi_max = std::min(dL*chiL,dR*chiR);
    return schmidt_with_norm(tensor,chi_max);
}

template<typename Scalar>
std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>,double  >
class_SVD::schmidt_with_norm(const Eigen::Tensor<Scalar,4> &tensor, long chi_max) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    if (dL*chiL * dR*chiR != tensor.size()){throw std::range_error("schmidt_with_norm error: tensor size does not match given dimensions.");}
    Eigen::Map<const MatrixType<Scalar>> mat (tensor.data(), dL*chiL, dR*chiR);
    auto [U,S,V,rank] = do_svd(mat,chi_max);
//    auto norm = S.norm();
//    auto Snormalized = S/norm;
//    std::cout << std::fixed << std::setprecision(16) << "regular norm: " << norm;
//    std::cout << std::fixed << std::setprecision(16) << "squared norm: " << S.squaredNorm();


    return std::make_tuple(Textra::MatrixTensorMap(U, dL, chiL, rank),
                           Textra::MatrixTensorMap(S.normalized().template cast<Scalar>(), rank),
                           Textra::MatrixTensorMap(V,  rank, dR, chiR ).shuffle(Textra::array3{ 1, 0, 2 }),
                           S.norm()
    );

}



