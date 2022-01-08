#pragma once
#include "common.h"
#include <array>
#include <fmt/core.h>
#include <unsupported/Eigen/CXX11/Tensor>

namespace linalg::tensor {

    // Shorthand for the list of index pairs.
    template<auto N>
    using idxlistpair = std::array<Eigen::IndexPair<Eigen::Index>, N>;

    inline constexpr idxlistpair<0> idx() { return {}; }

    template<std::size_t N, typename idxType>
    constexpr idxlistpair<N> idx(const idxType (&list1)[N], const idxType (&list2)[N]) {
        // Use numpy-style indexing for contraction. Each list contains a list of indices to be contracted for the respective
        // tensors. This function zips them together into pairs as used in Eigen::Tensor module. This does not sort the indices in decreasing order.
        idxlistpair<N> pairlistOut;
        for(size_t i = 0; i < N; i++) {
            pairlistOut[i] = Eigen::IndexPair<Eigen::Index>{static_cast<Eigen::Index>(list1[i]), static_cast<Eigen::Index>(list2[i])};
        }
        return pairlistOut;
    }

    template<typename T>
    std::string to_string(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors> &expr, int prec = 1, int width = 2, std::string_view sep = ", ") {
        using Evaluator = Eigen::TensorEvaluator<const Eigen::TensorForcedEvalOp<const T>, Eigen::DefaultDevice>;
        using Scalar    = typename Eigen::internal::remove_const<typename Evaluator::Scalar>::type;

        // Evaluate the expression if needed
        Eigen::TensorForcedEvalOp<const T> eval = expr.eval();
        Evaluator                          tensor(eval, Eigen::DefaultDevice());
        tensor.evalSubExprsIfNeeded(NULL);
        Eigen::Index total_size = Eigen::internal::array_prod(tensor.dimensions());

        if(total_size > 0 and tensor.dimensions().size() > 0) {
            Eigen::Index first_dim = tensor.dimensions()[0];
            if constexpr(T::NumDimensions == 4) first_dim = tensor.dimensions()[0] * tensor.dimensions()[1];
            if constexpr(T::NumDimensions == 6) first_dim = tensor.dimensions()[0] * tensor.dimensions()[1] * tensor.dimensions()[2];
            if constexpr(T::NumDimensions == 8) first_dim = tensor.dimensions()[0] * tensor.dimensions()[1] * tensor.dimensions()[2] * tensor.dimensions()[3];
            Eigen::Index other_dim = total_size / first_dim;
            auto matrix = Eigen::Map<Eigen::Array<typename T::Scalar, Eigen::Dynamic, Eigen::Dynamic, Evaluator::Layout>>(tensor.data(), first_dim, other_dim);
            std::string str;

            int comma = 1;
            if constexpr(std::is_integral_v<Scalar>) {
                comma = 0;
                prec  = 0;
            }
            if constexpr(linalg::is_std_complex_v<Scalar>)
                if constexpr(std::is_integral_v<typename Scalar::value_type>) {
                    comma = 0;
                    prec  = 0;
                }

            auto max_val   = static_cast<double>(matrix.cwiseAbs().maxCoeff());
            int  min_width = std::max(width, static_cast<int>(1 + std::max(0.0, std::log10(max_val))) + comma + prec);

            int min_width_real = min_width;
            int min_width_imag = min_width;
            if constexpr(linalg::is_std_complex_v<Scalar>) {
                auto max_val_real = static_cast<double>(matrix.real().cwiseAbs().maxCoeff());
                auto max_val_imag = static_cast<double>(matrix.imag().cwiseAbs().maxCoeff());
                min_width_real    = std::max(width, static_cast<int>(1 + std::max(0.0, std::log10(max_val_real))) + comma + prec);
                min_width_imag    = std::max(width, static_cast<int>(1 + std::max(0.0, std::log10(max_val_imag))) + comma + prec);
                if(matrix.real().minCoeff() < 0) min_width_real += 1;
                if(matrix.imag().minCoeff() < 0) min_width_imag += 1;
            } else {
                if(matrix.minCoeff() < 0) min_width += 1;
            }

            for(long i = 0; i < first_dim; i++) {
                str += fmt::format("[");
                for(long j = 0; j < other_dim; j++) {
                    if constexpr(linalg::is_std_complex_v<Scalar>) {
                        if constexpr(std::is_floating_point_v<typename Scalar::value_type>) {
                            std::string real = fmt::format("({0:.{1}f}", matrix(i, j).real(), prec);
                            std::string imag = fmt::format("{0:.{1}f})", matrix(i, j).imag(), prec);
                            std::string cplx = fmt::format("{:>},{:<}", real, imag);
                            str += fmt::format("{0:>{1}}", cplx, min_width_real + min_width_imag + 3); // Two doubles, comma, and parentheses

                        } else if constexpr(std::is_integral_v<typename Scalar::value_type>) {
                            std::string real = fmt::format("({}", matrix(i, j).real());
                            std::string imag = fmt::format("{})", matrix(i, j).imag());
                            std::string cplx = fmt::format("{:>},{:<}", real, imag);
                            str += fmt::format("{0:>{1}}", cplx, min_width_real + min_width_imag + 3); // Two doubles, comma, and parentheses
                        }
                    } else if constexpr(std::is_floating_point_v<Scalar>)
                        str += fmt::format("{0:>{1}.{2}f}", matrix(i, j), min_width, prec);
                    else if constexpr(std::is_integral_v<Scalar>)
                        str += fmt::format("{0:>{1}}", matrix(i, j), min_width);
                    if(j < other_dim - 1) str += sep;
                }
                if(i < first_dim - 1)
                    str += fmt::format("]\n");
                else
                    str += fmt::format("]");
            }
            return str;
        } else
            return "[]";
    }

    template<typename Scalar = double>
    Eigen::Tensor<Scalar, 2> identity(const Eigen::Index &dim) {
        Eigen::Tensor<double, 1> tensor(dim);
        tensor.setConstant(1);
        return tensor.inflate(std::array<long, 1>{tensor.size() + 1}).reshape(std::array<long, 2>{tensor.size(), tensor.size()}).template cast<Scalar>();
    }

    template<typename Scalar, int rank>
    Eigen::Tensor<Scalar, rank> mirror(const Eigen::Tensor<Scalar, rank> &tensor) {
        /*
         Returns a mirrored tensor

         Example: Starting with A
                0 1 2 3
                | | | |
               [  A   ]
               | | | |
               4 5 6 7

         returns

                3 2 1 0
                | | | |
               [  A   ]
               | | | |
               7 6 5 4

         This is useful for compaitibility with kronecker products which gives results indexed right to left.

        */
        if constexpr(rank <= 2)
            return tensor;
        else {
            std::array<Eigen::Index, rank> shf_idx{};
            for(size_t i = 0; i < static_cast<size_t>(rank); i++) { shf_idx[i] = static_cast<Eigen::Index>(i); }
            std::reverse(shf_idx.begin(), shf_idx.begin() + rank / 2);
            std::reverse(shf_idx.begin() + rank / 2, shf_idx.end());
            return tensor.shuffle(shf_idx);
        }
    }

    template<typename Scalar, int rank1, int rank2>
    Eigen::Tensor<Scalar, rank1 + rank2> outer(const Eigen::Tensor<Scalar, rank1> &tensor1, const Eigen::Tensor<Scalar, rank2> &tensor2) {
        std::array<Eigen::IndexPair<Eigen::Index>, 0> idx{};
        return tensor1.contract(tensor2, idx);
    }
    template<typename Scalar>
    Eigen::Tensor<Scalar, 4> outer(const EigenMatrix<Scalar> &matrix1, const EigenMatrix<Scalar> &matrix2) {
        std::array<Eigen::IndexPair<Eigen::Index>, 0> idx{};
        auto tensor1 = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(matrix1.data(), matrix1.rows(), matrix1.cols());
        auto tensor2 = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(matrix2.data(), matrix2.rows(), matrix2.cols());
        return tensor1.contract(tensor2, idx);
    }

    template<typename Scalar, int rankA, int rankB>
    Eigen::Tensor<Scalar, rankA + rankB> kronecker(const Eigen::Tensor<Scalar, rankA> &tensorA, const Eigen::Tensor<Scalar, rankB> &tensorB) {
        /*
         Returns the equivalent kronecker product for a tensor following left-to-right index order

                0  1       0 1 2         0 1 4 5 6           0 1 2 3 4
                |  |       | | |         | | | | |           | | | | |
               [ A ]  ⊗  [  B  ]   =   [   AB    ]  ===>  [   AB    ] (shuffle 0,1,4,5,6,2,3,7,8,9)
               |  |       | | |          | | | | |          | | | | |
               2  3       3 4 5          2 3 7 8 9          5 6 7 8 9

         */
        constexpr Eigen::Index                  topA  = static_cast<Eigen::Index>(rankA) / 2;
        constexpr Eigen::Index                  topB  = static_cast<Eigen::Index>(rankB) / 2;
        constexpr Eigen::Index                  topAB = topA + topB;
        std::array<Eigen::Index, rankA + rankB> shf{};
        for(size_t i = 0; i < shf.size(); i++) {
            if(i < topAB) {
                if(i < topA)
                    shf[i] = static_cast<Eigen::Index>(i);
                else
                    shf[i] = static_cast<Eigen::Index>(i) + (rankA - topA);
            } else {
                if(i - topAB < topA)
                    shf[i] = static_cast<Eigen::Index>(i) - topB;
                else
                    shf[i] = shf[i] = static_cast<Eigen::Index>(i);
            }
        }
        return linalg::tensor::outer(tensorA, tensorB).shuffle(shf);
    }

    template<typename Scalar>
    Eigen::Tensor<Scalar, 4> kronecker(const EigenMatrix<Scalar> &A, const EigenMatrix<Scalar> &B) {
        Eigen::Tensor<Scalar, 2> tensorA = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(A.data(), A.rows(), A.cols());
        Eigen::Tensor<Scalar, 2> tensorB = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(B.data(), B.rows(), B.cols());
        return linalg::tensor::kronecker(tensorA, tensorB);
    }
    template<typename Scalar, int rankA>
    Eigen::Tensor<Scalar, rankA + 2> kronecker(const Eigen::Tensor<Scalar, rankA> &tensorA, const EigenMatrix<Scalar> &B) {
        Eigen::Tensor<Scalar, 2> tensorB = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(B.data(), B.rows(), B.cols());
        return linalg::tensor::kronecker(tensorA, tensorB);
    }
    template<typename Scalar, int rankB>
    Eigen::Tensor<Scalar, rankB + 2> kronecker(const EigenMatrix<Scalar> &A, const Eigen::Tensor<Scalar, rankB> &tensorB) {
        Eigen::Tensor<Scalar, 2> tensorA = Eigen::TensorMap<const Eigen::Tensor<const Scalar, 2>>(A.data(), A.rows(), A.cols());
        return linalg::tensor::kronecker(tensorA, tensorB);
    }

    template<int rank>
    Eigen::IndexPair<Eigen::Index> mirror_idx_pair(const Eigen::IndexPair<Eigen::Index> &idx_pair) {
        if constexpr(rank <= 2)
            return idx_pair;
        else {
            std::array<Eigen::Index, rank> shf_idx{};
            for(size_t i = 0; i < static_cast<size_t>(rank); i++) { shf_idx[i] = static_cast<Eigen::Index>(i); }
            std::reverse(shf_idx.begin(), shf_idx.begin() + rank / 2);
            std::reverse(shf_idx.begin() + rank / 2, shf_idx.end());
            Eigen::IndexPair<Eigen::Index> idx_pair_mirrored{};
            idx_pair_mirrored.first  = std::distance(shf_idx.begin(), std::find(shf_idx.begin(), shf_idx.end(), idx_pair.first));
            idx_pair_mirrored.second = std::distance(shf_idx.begin(), std::find(shf_idx.begin(), shf_idx.end(), idx_pair.second));
            return idx_pair_mirrored;
        }
    }

    template<typename Scalar, int rank, size_t npair>
    Eigen::Tensor<Scalar, rank - 2 * npair> trace(const Eigen::Tensor<Scalar, rank> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, npair> &idx_pair,
                                                  bool mirror = false) {
        /*
         * Returns the partial trace of a tensor
         * Note that the tensor given here may be mirrored!
         */
        static_assert(rank >= 2 * npair, "Rank must be large enough");
        if constexpr(npair == 1) {
            auto idx_pair_r = idx_pair[0];
            if(mirror) idx_pair_r = mirror_idx_pair<rank>(idx_pair[0]); // Find the corresponding indices if this tensor was mirrored

            // Collect indices and dimensions traced
            auto                        idx1 = static_cast<size_t>(idx_pair_r.first);
            auto                        idx2 = static_cast<size_t>(idx_pair_r.second);
            std::array<Eigen::Index, 2> dim_tr{tensor.dimension(idx1), tensor.dimension(idx2)};
            //            Eigen::Tensor<Scalar, rank - 2> result;
            Eigen::Tensor<Scalar, rank - 2> result2;
            if(dim_tr[0] != dim_tr[1]) throw std::runtime_error("Traced dimensions must be equal size");
            //
            //            auto t_trace1 = tid::tic_scope(fmt::format("1<rank-{},npair-{}>", rank, npair));
            //            result = linalg::tensor::identity<Scalar>(dim_tr[0]).contract(tensor, idx({1ul, 0ul}, {idx1, idx2}));
            //
            //
            //            auto t_trace2 = tid::tic_scope(fmt::format("2<rank-{},npair-{}>", rank, npair));
            //            result2 = tensor.trace(std::array<Eigen::Index,2>{idx_pair_r.first, idx_pair_r.second});

            return tensor.trace(std::array<Eigen::Index, 2>{idx_pair_r.first, idx_pair_r.second});
        } else if constexpr(npair == 2) {
            std::array<long, 2> pair1{idx_pair[1].first, idx_pair[1].second};
            std::array<long, 2> pair0{idx_pair[0].first, idx_pair[0].second};
            pair0[0] -= std::count_if(pair1.begin(), pair1.end(), [&pair0](auto i) { return i < pair0[0]; });
            pair0[1] -= std::count_if(pair1.begin(), pair1.end(), [&pair0](auto i) { return i < pair0[1]; });
            //            auto res1 = linalg::tensor::trace(tensor, idx({pair1[0]}, {pair1[1]}));
            //            return linalg::tensor::trace(res1, idx({pair0[0]}, {pair0[1]}));
            return tensor.trace(pair1).trace(pair0);
        } else if constexpr(npair == 3) {
            std::array<long, 2> pair2{idx_pair[2].first, idx_pair[2].second};
            std::array<long, 2> pair1{idx_pair[1].first, idx_pair[1].second};
            std::array<long, 2> pair0{idx_pair[0].first, idx_pair[0].second};
            pair1[0] -= std::count_if(pair2.begin(), pair2.end(), [&pair1](auto i) { return i < pair1[0]; });
            pair1[1] -= std::count_if(pair2.begin(), pair2.end(), [&pair1](auto i) { return i < pair1[1]; });
            pair0[0] -= std::count_if(pair2.begin(), pair2.end(), [&pair0](auto i) { return i < pair0[0]; });
            pair0[1] -= std::count_if(pair2.begin(), pair2.end(), [&pair0](auto i) { return i < pair0[1]; });
            auto res = linalg::tensor::trace(tensor, idx({pair2[0]}, {pair2[1]}));
            return linalg::tensor::trace(tensor, idx({pair0[0], pair1[0]}, {pair0[1], pair1[1]}));
        } else
            throw std::runtime_error("Trace not implemented");
    }


}