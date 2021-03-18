//
// Created by david on 2017-10-04.
//

#pragma once

#include "svd/settings.h"
#include <general/nmspc_tensor_extra.h>
#include <math/num.h>
#include <optional>
#include <tools/common/log.h>

class class_tic_toc;

namespace svd {
    inline std::shared_ptr<spdlog::logger> log;
    class solver {
        private:

        void copy_settings(const svd::settings &svd_settings);

        template<typename Scalar>
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        template<typename Scalar>
        using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

        template<typename Scalar>
        std::tuple<MatrixType<Scalar>, VectorType<Scalar>, MatrixType<Scalar>, long> do_svd_lapacke(const Scalar *mat_ptr, long rows, long cols,
                                                                                                    std::optional<long> rank_max = std::nullopt);

        template<typename Scalar>
        std::tuple<MatrixType<Scalar>, VectorType<Scalar>, MatrixType<Scalar>, long> do_svd(const Scalar *mat_ptr, long rows, long cols,
                                                                                            std::optional<long> rank_max = std::nullopt);

        template<typename Derived>
        std::tuple<MatrixType<typename Derived::Scalar>, VectorType<typename Derived::Scalar>, MatrixType<typename Derived::Scalar>, long>
            do_svd(const Eigen::DenseBase<Derived> &mat) {
            return do_svd(mat.derived().data(), mat.rows(), mat.cols());
        }

        public:
        solver();
        solver(const svd::settings & svd_settings);
        solver(std::optional<svd::settings> svd_settings);
        double threshold    = 1e-12;
        size_t switchsize   = 16;
        bool   use_lapacke = true;
        bool   use_bdc = true;

        static std::optional<long long> count;
        double truncation_error = 0;

        std::shared_ptr<class_tic_toc> t_wrk;
        std::shared_ptr<class_tic_toc> t_adj;
        std::shared_ptr<class_tic_toc> t_jac;
        std::shared_ptr<class_tic_toc> t_svd;

        void   setLogLevel(size_t logLevel);
        void   enableProfiling();
        void   disableProfiling();

        template<typename Scalar>
        Eigen::Tensor<Scalar, 2> pseudo_inverse(const Eigen::Tensor<Scalar, 2> &tensor);

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 2>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2>> decompose(const Eigen::Tensor<Scalar, 2> &tensor,
                                                                                                           std::optional<long> rank_max = std::nullopt) {
            auto [U, S, V, rank] = do_svd(tensor.data(), tensor.dimension(0), tensor.dimension(1), rank_max);
            return std::make_tuple(Textra::TensorMap(U), Textra::TensorMap(S.normalized().template cast<Scalar>()), Textra::TensorMap(V));
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 2>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2>>
            decompose(const Eigen::Tensor<Scalar, 3> &tensor, const long rows, const long cols, std::optional<long> rank_max = std::nullopt) {
            auto [U, S, V, rank] = do_svd(tensor.data(), rows, cols, rank_max);
            return std::make_tuple(Textra::TensorMap(U), Textra::TensorMap(S.normalized().template cast<Scalar>()), Textra::TensorMap(V));
        }

        template<typename Derived>
        std::tuple<MatrixType<typename Derived::Scalar>, VectorType<typename Derived::Scalar>, MatrixType<typename Derived::Scalar>>
            decompose(const Eigen::DenseBase<Derived> &matrix, std::optional<long> rank_max = std::nullopt) {
            auto [U, S, V, rank] = do_svd(matrix.derived().data(), matrix.rows(), matrix.cols(), rank_max);
            return std::make_tuple(U, S.normalized().template cast<typename Derived::Scalar>(), V);
        }

        template<typename Scalar, auto N>
        std::tuple<Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>>
            schmidt(const Eigen::Tensor<Scalar, N> &tensor, long dL, long chiL, long dR, long chiR, std::optional<long> rank_max = std::nullopt) {
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
            if(dL * chiL * dR * chiR != tensor.size()) throw std::range_error("schmidt error: tensor size does not match given dimensions.");
            auto [U, S, V, rank] = do_svd(tensor.data(), dL * chiL, dR * chiR, rank_max);
            return std::make_tuple(Textra::TensorMap(U, dL, chiL, rank), Textra::TensorMap(S.normalized().template cast<Scalar>(), rank),
                                   Textra::TensorMap(V, rank, dR, chiR).shuffle(Textra::array3{1, 0, 2}));
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>> schmidt(const Eigen::Tensor<Scalar, 4> &tensor,
                                                                                                         std::optional<long> rank_max = std::nullopt) {
            long dL   = tensor.dimension(0);
            long chiL = tensor.dimension(1);
            long dR   = tensor.dimension(2);
            long chiR = tensor.dimension(3);
            if(dL * chiL * dR * chiR != tensor.size()) throw std::range_error("schmidt error: tensor size does not match given dimensions.");
            auto [U, S, V, rank] = do_svd(tensor.data(), dL * chiL, dR * chiR, rank_max);
            return std::make_tuple(Textra::TensorMap(U, dL, chiL, rank), Textra::TensorMap(S.normalized().template cast<Scalar>(), rank),
                                   Textra::TensorMap(V, rank, dR, chiR).shuffle(Textra::array3{1, 0, 2}));
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>, double>
            schmidt_with_norm(const Eigen::Tensor<Scalar, 4> &tensor, std::optional<long> rank_max = std::nullopt) {
            long dL   = tensor.dimension(0);
            long chiL = tensor.dimension(1);
            long dR   = tensor.dimension(2);
            long chiR = tensor.dimension(3);
            if(dL * chiL * dR * chiR != tensor.size()) throw std::range_error("schmidt_with_norm error: tensor size does not match given dimensions.");
            auto [U, S, V, rank] = do_svd(tensor.data(), dL * chiL, dR * chiR, rank_max);
            return std::make_tuple(Textra::TensorMap(U, dL, chiL, rank), Textra::TensorMap(S.normalized().template cast<Scalar>(), rank),
                                   Textra::TensorMap(V, rank, dR, chiR).shuffle(Textra::array3{1, 0, 2}), S.norm());
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>>
            schmidt_multisite(const Eigen::Tensor<Scalar, 3> &tensor, long dL, long dR, long chiL, long chiR, std::optional<long> rank_max = std::nullopt) {
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

            Eigen::Tensor<Scalar, 4> tensor_for_schmidt = tensor.reshape(Textra::array4{dL, dR, chiL, chiR}).shuffle(Textra::array4{0, 2, 1, 3});
            return schmidt(tensor_for_schmidt, dL, chiL, dR, chiR, rank_max);
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>>
            schmidt_into_right_normalized(const Eigen::Tensor<Scalar, 3> &tensor, long dR, std::optional<long> rank_max = std::nullopt) {
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
             * Here V is an "M" matrix of type B
             *
             */
            if(tensor.dimension(0) % dR != 0) throw std::runtime_error("Tensor dim 0 is not divisible by the given spin dimension " + std::to_string(dR));

            long dLdR = tensor.dimension(0);
            long dL   = dLdR / dR;
            long chiL = tensor.dimension(1);
            long chiR = tensor.dimension(2);
            return schmidt_multisite(tensor, dL, dR, chiL, chiR, rank_max);
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>>
            schmidt_into_left_normalized(const Eigen::Tensor<Scalar, 3> &tensor, long dL, std::optional<long> rank_max = std::nullopt) {
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
             *        dL=d                            dR=d^(N-1)
             *
             * Here U is an "M" matrix of type A
             * Here V is an "M" matrix of type B
             *
             */
            if(tensor.dimension(0) % dL != 0) throw std::runtime_error("Tensor dim 0 is not divisible by the given spin dimension " + std::to_string(dL));

            long dLdR = tensor.dimension(0);
            long dR   = dLdR / dL;
            long chiL = tensor.dimension(1);
            long chiR = tensor.dimension(2);
            return schmidt_multisite(tensor, dL, dR, chiL, chiR, rank_max);
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 4>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2>> split_mpo_l2r(const Eigen::Tensor<Scalar, 4> &mpo) {
            /*
             * Compress an MPO left to right using SVD as described in https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.035129
             *
             *          (2) d
             *             |
             *  (0) m ---[mpo]--- (1) m
             *            |
             *        (3) d
             *
             * is shuffled into
             *
             *          (0) d
             *             |
             *  (2) m ---[mpo]--- (3) m
             *            |
             *        (1) d
             *
             * and reshaped like
             *
             * ddm (012) ---[mpo]--- (3) m
             *
             * and subjected to the typical SVD so mpo = USV. This is then reshaped back into
             *
             *             d
             *             |
             *     m ---[mpo]---  m'   m'---[S]---'m m'---[V]---m
             *            |
             *            d
             *
             * where hopefully m' < m and the transfer matrix T = SV is multiplied onto the mpo on the right later.
             *
             * To stablize the compression, it is useful to insert avgS *  1/avgS, where one factor is put into U and the other into S.
             *
             *
             */

            auto dim0    = mpo.dimension(2);
            auto dim1    = mpo.dimension(3);
            auto dim2    = mpo.dimension(0);
            auto dim3    = mpo.dimension(1);
            auto dim_ddm = dim0 * dim1 * dim2;
            Eigen::Tensor<double, 2> mpo_rank2 = mpo.shuffle(Textra::array4{2, 3, 0, 1}).reshape(Textra::array2{dim_ddm, dim3}).real();
            auto [U, S, V, rank]               = do_svd(mpo_rank2.data(), mpo_rank2.dimension(0), mpo_rank2.dimension(1));
            auto avgS                          = num::next_power_of_two<double>(S.mean()); // Nearest power of two larger than S.mean();
            U *= avgS;
            S /= avgS;
            /* clang-format off */
            return std::make_tuple(
                Textra::TensorMap(U).reshape(Textra::array4{dim0, dim1, dim2, rank}).shuffle(Textra::array4{2, 3, 0, 1}).template cast<Scalar>(),
                Textra::TensorMap(S).template cast<Scalar>(),
                Textra::TensorMap(V).template cast<Scalar>());
            /* clang-format off */
        }
        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 2>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 4>> split_mpo_r2l(const Eigen::Tensor<Scalar, 4> &mpo) {
            /*
             * Splits an MPO right to left using SVD as described in https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.035129
             *
             *          (2) d
             *             |
             *  (0) m ---[mpo]--- (1) m
             *            |
             *        (3) d
             *
             * is shuffled into
             *
             *          (1) d
             *             |
             *  (0) m ---[mpo]--- (3) m
             *            |
             *        (2) d
             *
             * and reshaped like
             *
             * d (0) ---[mpo]--- (123) ddm
             *
             * and subjected to the typical SVD so mpo = USV. This is then reshaped back into
             *
             *                                          d
             *                                          |
             *   m---[U]---m'   m'---[S]---'m   m' ---[mpo]---  m
             *                                          |
             *                                          d
             *
             * where hopefully m' < m and the transfer matrix T = US is multiplied onto the mpo on the left later.
             *
             * To stablize the compression, it is useful to insert avgS *  1/avgS, where one factor is put into V and the other into S.
             *
             *
             */

            auto dim0    = mpo.dimension(0);
            auto dim1    = mpo.dimension(2);
            auto dim2    = mpo.dimension(3);
            auto dim3    = mpo.dimension(1);
            auto dim_ddm = dim1 * dim2 * dim3;

            Eigen::Tensor<double, 2> mpo_rank2 = mpo.shuffle(Textra::array4{0, 2, 3, 1}).reshape(Textra::array2{dim0, dim_ddm}).real();
            auto [U, S, V, rank]               = do_svd(mpo_rank2.data(), mpo_rank2.dimension(0), mpo_rank2.dimension(1));
            auto avgS = num::next_power_of_two<double>(S.mean()); // Nearest power of two larger than S.mean();
            V *= avgS;                                            // Rescaled singular values
            S /= avgS;                                            // Rescaled singular values
            /* clang-format off */
            return std::make_tuple(
                Textra::TensorMap(U).template cast<Scalar>(),
                Textra::TensorMap(S).template cast<Scalar>(),
                Textra::TensorMap(V).reshape(Textra::array4{rank, dim1, dim2, dim3}).shuffle(Textra::array4{0, 3, 1, 2}).template cast<Scalar>());
            /* clang-format off */
        }
    };
}