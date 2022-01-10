#pragma once

#include "svd/settings.h"
#include <math/num.h>
#include <math/tenx.h>
#include <optional>
#include <tools/common/log.h>
namespace tid {
    class ur;
}

namespace svd {
    inline std::shared_ptr<spdlog::logger> log;
    namespace internal {
        // LAPACK uses internal workspace arrays which can be reused for the duration of the program.
        // Call clear() to recover this memory space
        void clear_lapack();
    }

    class solver {
        private:
        void copy_settings(const svd::settings &svd_settings);
        template<typename Scalar>
        using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
        template<typename Scalar>
        using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

        template<typename Scalar>
        void save_svd(const MatrixType<Scalar> &A, const MatrixType<Scalar> &U, const VectorType<Scalar> &S, const MatrixType<Scalar> &VT, long rank_max,
                      const std::string &lib, const std::vector<std::pair<std::string, std::string>> &details) const;

        template<typename Scalar>
        std::tuple<MatrixType<Scalar>, VectorType<Scalar>, MatrixType<Scalar>, long> do_svd_lapacke(const Scalar *mat_ptr, long rows, long cols,
                                                                                                    std::optional<long> rank_max = std::nullopt);

        template<typename Scalar>
        std::tuple<MatrixType<Scalar>, VectorType<Scalar>, MatrixType<Scalar>, long> do_svd_eigen(const Scalar *mat_ptr, long rows, long cols,
                                                                                                  std::optional<long> rank_max = std::nullopt);

        template<typename Scalar>
        std::tuple<MatrixType<Scalar>, VectorType<Scalar>, MatrixType<Scalar>, long> do_svd_rsvd(const Scalar *mat_ptr, long rows, long cols,
                                                                                                 std::optional<long> rank_max = std::nullopt);

        template<typename Scalar>
        std::tuple<MatrixType<Scalar>, VectorType<Scalar>, MatrixType<Scalar>, long> do_svd_ptr(const Scalar *mat_ptr, long rows, long cols,
                                                                                                std::optional<long> rank_max = std::nullopt);

        template<typename Scalar>
        void print_matrix(const Scalar *mat_ptr, long rows, long cols, long dec = 8);
        template<typename Scalar>
        void print_vector(const Scalar *vec_ptr, long size, long dec = 8);

        public:
        solver();
        solver(const svd::settings &svd_settings);
        solver(std::optional<svd::settings> svd_settings);
        double                          threshold      = 1e-12;
        size_t                          switchsize_bdc = 16;   // Use Jacobi algorithm when rows < switchsize_bdc and BDC otherwise
        size_t                          switchsize_rnd = 1024; // Use Randomized SVD algorithm when rows < switchsize_bdc and BDC otherwise
        SVDLib                          svd_lib        = SVDLib::lapacke;
        bool                            use_bdc        = true; // Use fast bi-diagonal divide and conquer algorithm if rows >= switchsize_bdc
        bool                            save_fail      = false;
        bool                            save_result    = false;
        bool                            benchmark      = false;
        static std::optional<long long> count;
        double                          truncation_error = 0;

        void setLogLevel(size_t logLevel);

        template<typename Scalar>
        Eigen::Tensor<Scalar, 2> pseudo_inverse(const Eigen::Tensor<Scalar, 2> &tensor);

        template<typename Derived>
        auto do_svd(const Eigen::DenseBase<Derived> &mat, std::optional<long> rank_max = std::nullopt) {
            return do_svd_ptr(mat.derived().data(), mat.rows(), mat.cols(), rank_max);
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 2>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2>> decompose(const Eigen::Tensor<Scalar, 2> &tensor,
                                                                                                           std::optional<long> rank_max = std::nullopt) {
            auto [U, S, V, rank] = do_svd_ptr(tensor.data(), tensor.dimension(0), tensor.dimension(1), rank_max);
            return std::make_tuple(tenx::TensorMap(U), tenx::TensorMap(S.normalized().template cast<Scalar>()), tenx::TensorMap(V));
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 2>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2>>
            decompose(const Eigen::Tensor<Scalar, 3> &tensor, const long rows, const long cols, std::optional<long> rank_max = std::nullopt) {
            auto [U, S, V, rank] = do_svd_ptr(tensor.data(), rows, cols, rank_max);
            return std::make_tuple(tenx::TensorMap(U), tenx::TensorMap(S.normalized().template cast<Scalar>()), tenx::TensorMap(V));
        }

        template<typename Derived>
        std::tuple<MatrixType<typename Derived::Scalar>, VectorType<typename Derived::Scalar>, MatrixType<typename Derived::Scalar>>
            decompose(const Eigen::DenseBase<Derived> &matrix, std::optional<long> rank_max = std::nullopt) {
            auto [U, S, V, rank] = do_svd_ptr(matrix.derived().data(), matrix.rows(), matrix.cols(), rank_max);
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
            auto [U, S, V, rank] = do_svd_ptr(tensor.data(), dL * chiL, dR * chiR, rank_max);
            return std::make_tuple(tenx::TensorMap(U, dL, chiL, rank), tenx::TensorMap(S.normalized().template cast<Scalar>(), rank),
                                   tenx::TensorMap(V, rank, dR, chiR).shuffle(tenx::array3{1, 0, 2}));
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>> schmidt(const Eigen::Tensor<Scalar, 4> &tensor,
                                                                                                         std::optional<long> rank_max = std::nullopt) {
            long dL   = tensor.dimension(0);
            long chiL = tensor.dimension(1);
            long dR   = tensor.dimension(2);
            long chiR = tensor.dimension(3);
            if(dL * chiL * dR * chiR != tensor.size()) throw std::range_error("schmidt error: tensor size does not match given dimensions.");
            auto [U, S, V, rank] = do_svd_ptr(tensor.data(), dL * chiL, dR * chiR, rank_max);
            return std::make_tuple(tenx::TensorMap(U, dL, chiL, rank), tenx::TensorMap(S.normalized().template cast<Scalar>(), rank),
                                   tenx::TensorMap(V, rank, dR, chiR).shuffle(tenx::array3{1, 0, 2}));
        }

        template<typename Scalar>
        std::tuple<Eigen::Tensor<Scalar, 3>, Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>, double>
            schmidt_with_norm(const Eigen::Tensor<Scalar, 4> &tensor, std::optional<long> rank_max = std::nullopt) {
            long dL   = tensor.dimension(0);
            long chiL = tensor.dimension(1);
            long dR   = tensor.dimension(2);
            long chiR = tensor.dimension(3);
            if(dL * chiL * dR * chiR != tensor.size()) throw std::range_error("schmidt_with_norm error: tensor size does not match given dimensions.");
            auto [U, S, V, rank] = do_svd_ptr(tensor.data(), dL * chiL, dR * chiR, rank_max);
            return std::make_tuple(tenx::TensorMap(U, dL, chiL, rank), tenx::TensorMap(S.normalized().template cast<Scalar>(), rank),
                                   tenx::TensorMap(V, rank, dR, chiR).shuffle(tenx::array3{1, 0, 2}), S.norm());
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

            Eigen::Tensor<Scalar, 4> tensor_for_schmidt = tensor.reshape(tenx::array4{dL, dR, chiL, chiR}).shuffle(tenx::array4{0, 2, 1, 3});
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
             * where the index order is in parentheses. This is reinterpreted as
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

            auto                     dim0      = mpo.dimension(2);
            auto                     dim1      = mpo.dimension(3);
            auto                     dim2      = mpo.dimension(0);
            auto                     dim3      = mpo.dimension(1);
            auto                     dim_ddm   = dim0 * dim1 * dim2;
            Eigen::Tensor<Scalar, 2> mpo_rank2 = mpo.shuffle(tenx::array4{2, 3, 0, 1}).reshape(tenx::array2{dim_ddm, dim3});
            auto [U, S, V, rank]               = do_svd_ptr(mpo_rank2.data(), mpo_rank2.dimension(0), mpo_rank2.dimension(1));
            auto avgS                          = num::next_power_of_two<double>(S.mean()); // Nearest power of two larger than S.mean();
            U *= avgS;
            S /= avgS;
            /* clang-format off */
            return std::make_tuple(
                tenx::TensorMap(U).reshape(tenx::array4{dim0, dim1, dim2, rank}).shuffle(tenx::array4{2, 3, 0, 1}).template cast<Scalar>(),
                tenx::TensorMap(S).template cast<Scalar>(),
                tenx::TensorMap(V).template cast<Scalar>());
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

            Eigen::Tensor<Scalar, 2> mpo_rank2 = mpo.shuffle(tenx::array4{0, 2, 3, 1}).reshape(tenx::array2{dim0, dim_ddm});
            auto [U, S, V, rank]               = do_svd_ptr(mpo_rank2.data(), mpo_rank2.dimension(0), mpo_rank2.dimension(1));
            auto avgS = num::next_power_of_two<double>(S.mean()); // Nearest power of two larger than S.mean();
            V *= avgS;                                            // Rescaled singular values
            S /= avgS;                                            // Rescaled singular values
            /* clang-format off */
            return std::make_tuple(
                tenx::TensorMap(U).template cast<Scalar>(),
                tenx::TensorMap(S).template cast<Scalar>(),
                tenx::TensorMap(V).reshape(tenx::array4{rank, dim1, dim2, dim3}).shuffle(tenx::array4{0, 3, 1, 2}).template cast<Scalar>());
            /* clang-format off */
        }
    };
}