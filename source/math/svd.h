//
// Created by david on 2017-10-04.
//

#pragma once

#include "general/nmspc_tensor_extra.h"
#include <iomanip>
#include <optional>
#include <io/nmspc_logger.h>

namespace svd {
    inline std::shared_ptr<spdlog::logger> log;
    class solver {
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
        solver(size_t logLevel = 2);
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
        schmidt(const Eigen::Tensor<Scalar,N> &tensor, long dL, long chiL,long dR, long chiR, std::optional<long> rank_max = std::nullopt){
            /*
             * When using this function, it is important that the index order is correct
             * The order is
             * 1) left physical index
             * 2) left bond index
             * 3) right physical index
             * 4) right bond index
             *
             * It is assumed that this order is built into tensor already. Hard to debug errors will occurr
             * if this is not the case!
             * The function call argument order dL, chiL, dR,chiR is meant as a hint for how to use this function.
             */

            if (dL*chiL * dR*chiR != tensor.size())
                throw std::range_error("schmidt error: tensor size does not match given dimensions.");
            auto [U,S,V,rank] = do_svd(tensor.data(),dL*chiL, dR*chiR,rank_max);

//            Eigen::Tensor<Scalar, 3> U_map = Textra::MatrixTensorMap(U, dL, chiL, rank);
//            Eigen::Tensor<Scalar, 1> S_map = Textra::MatrixTensorMap(S.normalized().template cast<Scalar>(), rank);
//            Eigen::Tensor<Scalar, 3> V_map = Textra::MatrixTensorMap(V,  rank, dR, chiR ).shuffle(Textra::array3{ 1, 0, 2 });
//            return std::make_tuple(U_map,
//            S_map,
//            V_map );

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
            if (dL*chiL * dR*chiR != tensor.size())
                throw std::range_error("schmidt error: tensor size does not match given dimensions.");
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
            if (dL*chiL * dR*chiR != tensor.size())
                throw std::range_error("schmidt_with_norm error: tensor size does not match given dimensions.");
            auto [U,S,V,rank] = do_svd(tensor.data(),dL*chiL, dR*chiR,rank_max);
            return std::make_tuple(Textra::MatrixTensorMap(U, dL, chiL, rank),
                                   Textra::MatrixTensorMap(S.normalized().template cast<Scalar>(), rank),
                                   Textra::MatrixTensorMap(V,  rank, dR, chiR ).shuffle(Textra::array3{ 1, 0, 2 }),
                                   S.norm());
        }



        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
        schmidt_multisite(const Eigen::Tensor<Scalar,3> &tensor, long dL, long dR, long chiL, long chiR, std::optional<long> rank_max = std::nullopt){
            /* This function assumes that the tensor is given with indices in multisite standard form, i.e., with order
             *
             * 1) physical indices (any number of them, contracted from left to right)
             * 2) left bond index
             * 4) right bond index
             *
             * (1)chiL ---[tensor]--- (2)chiR
             *               |
             *          (0)d*d*d...
             *
             * we start by transposing the tensor into "left-right" order suitable for schmidt decomposition
             *
             * (1)chiL---[  tensor  ]---(3)chiR
             *            |        |
             *        (0)dL      (2)dR
             *
             * The function call argument order dL, dR, chiL,chiR is meant as a hint for how to use this function.
             *
             */

            Eigen::Tensor<Scalar,4> tensor_for_schmidt = tensor.reshape(Textra::array4{dL,dR,chiL,chiR}).shuffle(Textra::array4{0,2,1,3});
            return schmidt(tensor_for_schmidt,dL,chiL,dR,chiR, rank_max);
        }


        template<typename Scalar>
        std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
        schmidt_from_right  (const Eigen::Tensor<Scalar,3> &tensor,long dR , std::optional<long> rank_max = std::nullopt){
            /* This schmidt decomposition is used to pull a site out of rank-3 tensor from the right
             * We obtain USV matrices from a tensor containing N sites with spin dim = d in two steps:
             *
             * (1)chiL ---[mps]--- (2)chiR
             *             |
             *        (0) dL*dR
             *
             * where the index order is in parentheses. This is is reinterpreted as
             *
             *
             * (2)chiL---[    mps    ]---(3)chiR
             *            |         |
             *  (0)dL=d^(N-1)   (1)dR=d
             *
             * and reshuffled into
             *
             * (1)chiL---[    mps    ]---(3)chiR
             *            |         |
             *  (0)dL=d^(N-1)   (2)dR=d
             *
             * which then is split up into
             *
             * chiL ---[U]----chi chi---[S]---chi chi---[V]---chiR
             *          |                                |
             *     d^(N-1)=dL                           d=dR
             *
             * Here V is an "M" matrix of type B = Gamma * Lambda
             *
             */
            if(tensor.dimension(0) % dR != 0)
                throw std::runtime_error("Tensor dim 0 is not divisible by the given spin dimension " + std::to_string(dR));

            long dLdR = tensor.dimension(0);
            long dL   = dLdR/dR;
            long chiL = tensor.dimension(1);
            long chiR = tensor.dimension(2);
            return schmidt_multisite(tensor,dL,dR,chiL,chiR,rank_max);
        }



        template<typename Scalar>
        std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
        schmidt_from_left  (const Eigen::Tensor<Scalar,3> &tensor,long dL , std::optional<long> rank_max = std::nullopt){
            /* This schmidt decomposition is used to pull a site out of rank-3 tensor from the left
             * We obtain USV matrices from a tensor containing N sites with spin dim = d in two steps:
             *
             * (1)chiL ---[mps]--- (2)chiR
             *             |
             *        (0) dL*dR
             *
             * where the index order is in parentheses. This is is reinterpreted as
             *
             *
             * (2)chiL---[    mps    ]---(3)chiR
             *            |         |
             *        (0)dL=d  (1)dR=d^(N-1)
             *
             * and reshuffled into
             *
             * (1)chiL---[    mps    ]---(3)chiR
             *            |         |
             *        (0)dL=d   (2)dR=d^(N-1)
             *
             * which then is split up into
             *
             * chiL ---[U]----chi chi---[S]---chi chi---[V]---chiR
             *          |                                |
             *     d^(N-1)=dL                           d=dR
             *
             * Here V is an "M" matrix of type B = Gamma * Lambda
             *
             */
            if(tensor.dimension(0) % dL != 0)
                throw std::runtime_error("Tensor dim 0 is not divisible by the given spin dimension " + std::to_string(dL));


            long dLdR = tensor.dimension(0);
            long dR   = dLdR/dL;
            long chiL = tensor.dimension(1);
            long chiR = tensor.dimension(2);
            return schmidt_multisite(tensor,dL,dR,chiL,chiR,rank_max);
        }
    };
}
