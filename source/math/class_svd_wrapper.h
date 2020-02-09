//
// Created by david on 2017-10-04.
//

#pragma once

#include "general/nmspc_tensor_extra.h"
#include <iomanip>
#include <optional>

class class_SVD{
private:
    double SVDThreshold         = 1e-12;
    double truncation_error     = 0;
    template <typename Scalar> using MatrixType = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>;
    template <typename Scalar> using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;

    template<typename Scalar>
    std::tuple<MatrixType<Scalar>, VectorType<Scalar>,MatrixType<Scalar>, long>
    do_svd_lapacke(const Scalar * mat_ptr, long rows, long cols, std::optional<long> rank_max = std::nullopt);


    template<typename Scalar>
    std::tuple<MatrixType<Scalar>, VectorType<Scalar>,MatrixType<Scalar>, long>
    do_svd(const Scalar * mat_ptr, long rows, long cols, std::optional<long> rank_max = std::nullopt);

    template<typename Derived>
    std::tuple<MatrixType<typename Derived::Scalar>, VectorType<typename Derived::Scalar>,MatrixType<typename Derived::Scalar>, long>
    do_svd(const Eigen::DenseBase<Derived> & mat, std::optional<long> rank_max = std::nullopt){
        if(not rank_max.has_value()) rank_max = std::min(mat.rows(),mat.cols());
        return do_svd(mat.derived().data(), mat.rows(),mat.cols(),rank_max);
    }

public:

    class_SVD()=default;
    bool   use_lapacke          = false;
    double get_truncation_error();
    void setThreshold(double newThreshold);

    template<typename Scalar>
    Eigen::Tensor<Scalar, 2> pseudo_inverse(const Eigen::Tensor<Scalar,2> &tensor);







    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,2> &tensor, std::optional<long> rank_max = std::nullopt) {
        auto[U,S,V,rank] = do_svd(tensor.data(),tensor.dimension(0), tensor.dimension(1), rank_max );
        return std::make_tuple(Textra::MatrixTensorMap(U),
                               Textra::MatrixTensorMap(S.normalized().template cast<Scalar>()),
                               Textra::MatrixTensorMap(V)
        );
    }



    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,3> &tensor, const long rows,const long cols, std::optional<long> rank_max = std::nullopt){
        auto[U,S,V,rank] = do_svd(tensor.data(),rows,cols, rank_max );
        return std::make_tuple(Textra::MatrixTensorMap(U),
                               Textra::MatrixTensorMap(S.normalized().template cast<Scalar>()),
                               Textra::MatrixTensorMap(V)
        );
    }

    template<typename Derived>
    std::tuple<MatrixType<typename Derived::Scalar>, VectorType<typename Derived::Scalar>,MatrixType<typename Derived::Scalar>>
    decompose(const Eigen::DenseBase<Derived> &matrix, std::optional<long> rank_max = std::nullopt){
        auto[U,S,V,rank] = do_svd(matrix.derived().data(),matrix.rows(),matrix.cols(),rank_max);
        return std::make_tuple(U,
                               S.normalized().template cast<typename Derived::Scalar>(),
                               V
        );
    }


    template<typename Scalar,auto N>
    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,N> &tensor, long dL, long dR, long chiL, long chiR, std::optional<long> rank_max = std::nullopt){
        if (dL*chiL * dR*chiR != tensor.size()){throw std::range_error("schmidt error: tensor size does not match given dimensions.");}
        auto [U,S,V,rank] = do_svd(tensor.data(),dL*chiL, dR*chiR,rank_max);
        return std::make_tuple(Textra::MatrixTensorMap(U, dL, chiL, rank),
                               Textra::MatrixTensorMap(S.normalized().template cast<Scalar>(), rank),
                               Textra::MatrixTensorMap(V,  rank, dR, chiR ).shuffle(Textra::array3{ 1, 0, 2 })
        );
    }


    template<typename Scalar>
    std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,4> &tensor, std::optional<long> rank_max = std::nullopt){
        long dL   = tensor.dimension(0);
        long chiL = tensor.dimension(1);
        long dR   = tensor.dimension(2);
        long chiR = tensor.dimension(3);
        if (dL*chiL * dR*chiR != tensor.size()){throw std::range_error("schmidt error: tensor size does not match given dimensions.");}
        auto [U,S,V,rank] = do_svd(tensor.data(),dL*chiL, dR*chiR,rank_max);
        return std::make_tuple(Textra::MatrixTensorMap(U, dL, chiL, rank),
                               Textra::MatrixTensorMap(S.normalized().template cast<Scalar>(), rank),
                               Textra::MatrixTensorMap(V,  rank, dR, chiR ).shuffle(Textra::array3{ 1, 0, 2 })
        );
    }


    template<typename Scalar>
    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>, double >
    schmidt_with_norm  (const Eigen::Tensor<Scalar,4> &tensor,std::optional<long> rank_max = std::nullopt) {
        long dL   = tensor.dimension(0);
        long chiL = tensor.dimension(1);
        long dR   = tensor.dimension(2);
        long chiR = tensor.dimension(3);
        if (dL*chiL * dR*chiR != tensor.size()){throw std::range_error("schmidt_with_norm error: tensor size does not match given dimensions.");}
        auto [U,S,V,rank] = do_svd(tensor.data(),dL*chiL, dR*chiR,rank_max);
        return std::make_tuple(Textra::MatrixTensorMap(U, dL, chiL, rank),
                               Textra::MatrixTensorMap(S.normalized().template cast<Scalar>(), rank),
                               Textra::MatrixTensorMap(V,  rank, dR, chiR ).shuffle(Textra::array3{ 1, 0, 2 }),
                               S.norm()
        );

    }


};





