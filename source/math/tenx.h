#pragma once
#ifndef lapack_complex_float
    #define lapack_complex_float std::complex<float>
#endif
#ifndef lapack_complex_double
    #define lapack_complex_double std::complex<double>
#endif

#include "tenx/eval.h"
#include "tenx/omp.h"
#include "tenx/sfinae.h"
#include "tenx/span.h"
#include <array>
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

/*! \brief **Textra** stands for "Tensor Extra". Provides extra functionality to Eigen::Tensor.*/

/*!
 *  \namespace Textra
 *  This namespace makes shorthand typedef's to Eigen's unsupported Tensor module, and provides handy functions
 *  to interface between `Eigen::Tensor` and `Eigen::Matrix` objects.
 *  The contents of this namespace is co clear it is self-documenting ;)
 */

/*clang-format off */
namespace tenx {
    template<typename Scalar>
    using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    template<typename Scalar>
    using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    template<Eigen::Index rank>
    using array = std::array<Eigen::Index, rank>;

    using array8 = array<8>;
    using array7 = array<7>;
    using array6 = array<6>;
    using array5 = array<5>;
    using array4 = array<4>;
    using array3 = array<3>;
    using array2 = array<2>;
    using array1 = array<1>;
    using array0 = array<0>;

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

    template<std::size_t N, typename idxType>
    constexpr idxlistpair<N> idx(const std::array<idxType, N> &list1, const std::array<idxType, N> &list2) {
        // Use numpy-style indexing for contraction. Each list contains a list of indices to be contracted for the respective
        // tensors. This function zips them together into pairs as used in Eigen::Tensor module. This does not sort the indices in decreasing order.
        idxlistpair<N> pairlistOut;
        for(size_t i = 0; i < N; i++) {
            pairlistOut[i] = Eigen::IndexPair<Eigen::Index>{static_cast<Eigen::Index>(list1[i]), static_cast<Eigen::Index>(list2[i])};
        }
        return pairlistOut;
    }

    struct idx_dim_pair {
        Eigen::Index idxA;
        Eigen::Index idxB;
        Eigen::Index dimB;
    };

    template<std::size_t NB, std::size_t N>
    constexpr idxlistpair<N> sortIdx(const std::array<Eigen::Index, NB> &dimensions, const Eigen::Index (&idx_ctrct_A)[N],
                                     const Eigen::Index (&idx_ctrct_B)[N]) {
        // When doing contractions, some indices may be larger than others. For performance, you want to
        // contract the largest indices first. This will return a sorted index list in decreasing order.
        std::array<idx_dim_pair, N> idx_dim_pair_list;
        for(size_t i = 0; i < N; i++) { idx_dim_pair_list[i] = {idx_ctrct_A[i], idx_ctrct_B[i], dimensions[idx_ctrct_B[i]]}; }
        std::sort(idx_dim_pair_list.begin(), idx_dim_pair_list.end(), [](const auto &i, const auto &j) { return i.dimB > j.dimB; });
        idxlistpair<N> pairlistOut;
        for(size_t i = 0; i < N; i++) { pairlistOut[i] = Eigen::IndexPair<long>{idx_dim_pair_list[i].idxA, idx_dim_pair_list[i].idxB}; }
        return pairlistOut;
    }

    //
    //    //***************************************//
    //    //Different views for rank 1 and 2 tensors//
    //    //***************************************//
    //

    template<typename T, typename Device = Eigen::DefaultDevice>
    auto extractDiagonal(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors> &expr, const Device &device = Device()) {
        auto tensor    = tenx::asEval(expr, device);
        auto tensorMap = tensor.map();
        auto size      = tensorMap.size();
        auto rows      = static_cast<Eigen::Index>(std::floor(std::sqrt(size)));
        return tensorMap.reshape(array1{size}).stride(array1{rows + 1});
    }

    template<typename T, typename Device = Eigen::DefaultDevice>
    auto asDiagonal(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors> &expr, const Device &device = Device()) {
        auto tensor    = tenx::asEval(expr, device);
        auto tensorMap = tensor.map();
        auto size      = tensorMap.size();
        using Scalar   = typename decltype(tensorMap)::Scalar;
        return static_cast<Eigen::Tensor<Scalar, 2>>(tensorMap.inflate(array1{size + 1}).reshape(array2{size, size}));
    }

    template<typename Scalar>
    Eigen::Tensor<Scalar, 2> asDiagonalSquared(const Eigen::Tensor<Scalar, 1> &tensor) {
        return tensor.square().inflate(array1{tensor.size() + 1}).reshape(array2{tensor.size(), tensor.size()});
    }
    template<typename Scalar>
    Eigen::Tensor<Scalar, 2> asDiagonalSqrt(const Eigen::Tensor<Scalar, 1> &tensor) {
        return tensor.sqrt().inflate(array1{tensor.size() + 1}).reshape(array2{tensor.size(), tensor.size()});
    }

    template<typename Scalar>
    Eigen::Tensor<Scalar, 2> asDiagonalInversed(const Eigen::Tensor<Scalar, 1> &tensor) {
        return tensor.inverse().inflate(array1{tensor.size() + 1}).reshape(array2{tensor.size(), tensor.size()});
    }

    template<typename Scalar>
    Eigen::Tensor<Scalar, 2> asDiagonalInversed(const Eigen::Tensor<Scalar, 2> &tensor) {
        if(not tensor.dimension(0) == tensor.dimension(1)) throw std::runtime_error("tenx::asDiagonalInversed expects a square tensor");
        return asDiagonalInversed(extractDiagonal(tensor));
    }

    template<typename T, typename Device = Eigen::DefaultDevice>
    double norm(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors> &expr, const Device &device = Device()) {
        return asEval(expr.square().sum().sqrt().abs(), device)->coeff(0);
    }

    inline Eigen::Tensor<std::complex<double>, 1> broadcast(Eigen::Tensor<std::complex<double>, 1> &tensor, const std::array<long, 1> &bcast) {
        // Use this function to avoid a bug in Eigen when broadcasting complex tensors of rank 1, with compiler option -mfma
        // See more here https://gitlab.com/libeigen/eigen/-/issues/2351
        std::array<long, 2> bcast2 = {bcast[0], 1};
        std::array<long, 2> shape2 = {tensor.size(), 1};
        std::array<long, 1> shape1 = {tensor.size() * bcast[0]};
        return tensor.reshape(shape2).broadcast(bcast2).reshape(shape1);
    }

    template<typename T, typename Device = Eigen::DefaultDevice>
    auto asNormalized(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors> &expr, const Device &device = Device()) {
        auto           tensor    = tenx::asEval(expr, device);
        auto           tensorMap = tensor.map();
        constexpr auto rank      = tensor.rank();
        using Scalar             = typename decltype(tensorMap)::Scalar;
        using TensorType         = typename Eigen::Tensor<Scalar, rank>;

        Eigen::Tensor<double, 0> normInverse = tensorMap.square().sum().sqrt().abs().inverse();
        return static_cast<TensorType>(tensorMap * tensorMap.constant(normInverse.coeff(0)));
    }

    template<typename Scalar, auto rank>
    void normalize(Eigen::Tensor<Scalar, rank> &tensor) {
        Eigen::Map<VectorType<Scalar>> map(tensor.data(), tensor.size());
        map.normalize();
    }

    template<typename Scalar, auto rank>
    auto TensorConstant(Scalar constant, const std::array<Eigen::Index, rank> &dims) {
        Eigen::Tensor<Scalar, rank> tensor(dims);
        tensor.setConstant(constant);
        return tensor;
    }

    template<typename Scalar, typename... Dims>
    auto TensorConstant(Scalar constant, const Dims... dims) {
        return TensorConstant(constant, std::array<Eigen::Index, sizeof...(Dims)>{dims...});
    }

    template<typename Scalar, auto rank>
    auto TensorRandom(const std::array<Eigen::Index, rank> &dims) {
        Eigen::Tensor<Scalar, rank> tensor(dims);
        tensor.setRandom();
        return tensor;
    }

    template<typename Scalar, typename... Dims>
    auto TensorRandom(const Dims... dims) {
        return TensorRandom<Scalar>(std::array<Eigen::Index, sizeof...(Dims)>{dims...});
    }

    template<typename Scalar>
    auto TensorIdentity(long diagSize) {
        Eigen::Tensor<Scalar, 1> tensor(diagSize);
        tensor.setConstant(1.0);
        return asDiagonal(tensor);
    }

    template<typename Scalar, auto rank>
    Eigen::Tensor<Scalar, rank - 2> trace(const Eigen::Tensor<Scalar, rank> &tensor, const idxlistpair<1> &idx_pair) {
        static_assert(rank >= 2, "Rank must be >= 2 for trace of an index pair");
        long idx0 = idx_pair[0].first;
        long idx1 = idx_pair[0].second;
        printf("Tracing pair %ld %ld of rank %ld tensor", idx0, idx1, rank);
        if(tensor.dimension(idx0) != tensor.dimension(idx1)) throw std::logic_error("Can't trace index pair of different dimensions");
        long                     dim0 = tensor.dimension(idx0);
        Eigen::Tensor<Scalar, 1> id(dim0);
        id.setConstant(1.0);
        return asDiagonal(id).contract(tensor, tenx::idx({0l, 1l}, {idx0, idx1}));
    }

    template<typename Scalar, auto rank>
    Eigen::Tensor<Scalar, rank - 4> trace(const Eigen::Tensor<Scalar, rank> &tensor, const idxlistpair<2> &idx_pair) {
        static_assert(rank >= 4, "Rank must be >= 4 for trace of 2 index pairs");
        long idx00 = idx_pair[0].first;
        long idx01 = idx_pair[0].second;
        long idx10 = idx_pair[1].first;
        long idx11 = idx_pair[1].second;
        // idx_pair[0] connect to each other, and then idx_pair[1] connect to each other
        if(tensor.dimension(idx00) != tensor.dimension(idx01)) throw std::logic_error("Can't trace index pair of different dimensions");
        if(tensor.dimension(idx10) != tensor.dimension(idx11)) throw std::logic_error("Can't trace index pair of different dimensions");
        long                     dim00 = tensor.dimension(idx00);
        long                     dim11 = tensor.dimension(idx11);
        Eigen::Tensor<Scalar, 1> id00(dim00), id11(dim11);
        id00.setConstant(1.0);
        id11.setConstant(1.0);
        return asDiagonal(id00)
            .contract(asDiagonal(id11), tenx::idx()) // Outer product
            .contract(tensor, tenx::idx({0l, 1l, 2l, 3l}, {idx00, idx01, idx10, idx11}));
    }

    //    //****************************//
    //    //Matrix to tensor conversions//
    //    //****************************//

    // Detects if Derived is a plain object, like "MatrixXd" or similar.
    // std::decay removes pointer or ref qualifiers if present
    template<typename Derived>
    using is_plainObject = std::is_base_of<Eigen::PlainObjectBase<std::decay_t<Derived>>, std::decay_t<Derived>>;

    template<typename Derived, typename T, auto rank>
    Eigen::Tensor<typename Derived::Scalar, rank> TensorCast(const Eigen::EigenBase<Derived> &matrix, const std::array<T, rank> &dims) {
        return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>(matrix.derived().eval().data(), dims);
    }

    template<typename Derived, typename T, auto rank>
    Eigen::Tensor<typename Derived::Scalar, rank> TensorCast(const Eigen::EigenBase<Derived> &matrix, const Eigen::DSizes<T, rank> &dims) {
        return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>(matrix.derived().eval().data(), dims);
    }

    // Helpful overload
    template<typename Derived, typename... Dims>
    auto TensorCast(const Eigen::EigenBase<Derived> &matrix, const Dims... dims) {
        static_assert(sizeof...(Dims) > 0, "TensorCast: sizeof... (Dims) must be larger than 0");
        return TensorCast(matrix, std::array<Eigen::Index, sizeof...(Dims)>{dims...});
    }

    template<typename Derived>
    auto TensorCast(const Eigen::EigenBase<Derived> &matrix) {
        if constexpr(Derived::ColsAtCompileTime == 1 or Derived::RowsAtCompileTime == 1) {
            return TensorCast(matrix, matrix.size());
        } else {
            return TensorCast(matrix, matrix.rows(), matrix.cols());
        }
    }

    template<typename Derived, auto rank>
    auto TensorMap(const Eigen::PlainObjectBase<Derived> &matrix, const std::array<long, rank> &dims) {
        return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>(matrix.derived().data(), dims);
    }

    template<typename Derived, typename... Dims>
    auto TensorMap(const Eigen::PlainObjectBase<Derived> &matrix, const Dims... dims) {
        return TensorMap(matrix, std::array<long, static_cast<int>(sizeof...(Dims))>{dims...});
    }

    template<typename Derived>
    auto TensorMap(const Eigen::PlainObjectBase<Derived> &matrix) {
        if constexpr(Derived::ColsAtCompileTime == 1 or Derived::RowsAtCompileTime == 1) {
            return TensorMap(matrix, matrix.size());
        } else {
            return TensorMap(matrix, matrix.rows(), matrix.cols());
        }
    }

    //
    //    //****************************//
    //    //Tensor to matrix conversions//
    //    //****************************//
    //

    template<typename T, typename sizeType, typename Device = Eigen::DefaultDevice>
    auto MatrixCast(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors> &expr, const sizeType rows, const sizeType cols, const Device &device = Device()) {
        auto tensor    = asEval(expr, device);
        auto tensorMap = tensor.map();
        using Scalar   = typename decltype(tensorMap)::Scalar;
        return static_cast<MatrixType<Scalar>>(Eigen::Map<const MatrixType<Scalar>>(tensorMap.data(), rows, cols));
    }

    template<typename T, typename Device = Eigen::DefaultDevice>
    auto VectorCast(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors> &expr, const Device &device = Device()) {
        auto tensor    = asEval(expr, device);
        auto tensorMap = tensor.map();
        auto size      = Eigen::internal::array_prod(tensorMap.dimensions());
        using Scalar   = typename decltype(tensorMap)::Scalar;
        return static_cast<VectorType<Scalar>>(Eigen::Map<const VectorType<Scalar>>(tensorMap.data(), size));
    }

    template<typename Scalar, auto rank, typename sizeType>
    auto MatrixMap(const Eigen::Tensor<Scalar, rank> &tensor, const sizeType rows, const sizeType cols) {
        return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), rows, cols);
    }

    template<typename Scalar, auto rank, typename sizeType>
    auto MatrixMap(const Eigen::Tensor<Scalar, rank> &&tensor, const sizeType rows, const sizeType cols) = delete; // Prevent map from temporary

    template<typename Scalar>
    auto MatrixMap(const Eigen::Tensor<Scalar, 2> &tensor) {
        return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    }
    template<typename Scalar>
    auto MatrixMap(const Eigen::Tensor<Scalar, 2> &&tensor) = delete; // Prevent map from temporary

    template<typename Scalar, auto rank>
    auto VectorMap(const Eigen::Tensor<Scalar, rank> &tensor) {
        return Eigen::Map<const VectorType<Scalar>>(tensor.data(), tensor.size());
    }
    template<typename Scalar, auto rank>
    auto VectorMap(const Eigen::Tensor<Scalar, rank> &&tensor) = delete; // Prevent map from temporary

    //************************//
    // change storage layout //
    //************************//
    template<typename Scalar, auto rank>
    Eigen::Tensor<Scalar, rank, Eigen::RowMajor> to_RowMajor(const Eigen::Tensor<Scalar, rank, Eigen::ColMajor> tensor) {
        std::array<long, rank> neworder;
        std::iota(std::begin(neworder), std::end(neworder), 0);
        std::reverse(neworder.data(), neworder.data() + neworder.size());
        return tensor.swap_layout().shuffle(neworder);
    }

    template<typename Derived>
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> to_RowMajor(const Eigen::MatrixBase<Derived> &matrix) {
        if(matrix.IsRowMajor) { return matrix; }
        Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrowmajor = matrix;
        return matrowmajor;
    }

    template<typename Scalar, auto rank>
    Eigen::Tensor<Scalar, rank, Eigen::ColMajor> to_ColMajor(const Eigen::Tensor<Scalar, rank, Eigen::RowMajor> tensor) {
        std::array<long, rank> neworder;
        std::iota(std::begin(neworder), std::end(neworder), 0);
        std::reverse(neworder.data(), neworder.data() + neworder.size());
        return tensor.swap_layout().shuffle(neworder);
    }

    template<typename Derived>
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> to_ColMajor(const Eigen::MatrixBase<Derived> &matrix) {
        if(not matrix.IsRowMajor) { return matrix; }
        Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> matrowmajor = matrix;
        return matrowmajor;
    }

    //******************************************************//
    // Tests for real/complex //
    //******************************************************//

    template<typename Derived>
    bool isPositive(const Eigen::EigenBase<Derived> &obj) {
        return (obj.derived().array().real() >= 0).all();
    }

    template<typename Scalar, auto rank>
    bool isPositive(const Eigen::Tensor<Scalar, rank> &tensor) {
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> vector(tensor.data(), tensor.size());
        return isPositive(vector);
    }

    template<typename Derived>
    bool isReal(const Eigen::EigenBase<Derived> &obj, double threshold = std::numeric_limits<double>::epsilon()) {
        using Scalar = typename Derived::Scalar;
        if constexpr(sfinae::is_std_complex_v<Scalar>) {
            auto imag_sum = obj.derived().imag().cwiseAbs().sum();
            threshold *= std::max<double>(1.0, static_cast<double>(obj.derived().size()));
            //            if(imag_sum >= threshold) {
            //                std::printf("thr*size : %.20f imag_sum : %.20f | isreal %d \n", threshold, imag_sum, imag_sum < threshold);
            //            }
            return imag_sum < threshold;

        } else {
            return true;
        }
    }

    template<typename Scalar, auto rank>
    bool isReal(const Eigen::Tensor<Scalar, rank> &tensor, double threshold = std::numeric_limits<double>::epsilon()) {
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> vector(tensor.data(), tensor.size());
        return isReal(vector, threshold);
    }

    template<typename Scalar, auto rank>
    bool isZero(const Eigen::Tensor<Scalar, rank> &tensor, double threshold = std::numeric_limits<double>::epsilon()) {
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> vector(tensor.data(), tensor.size());
        return vector.isZero(threshold);
    }

    template<typename Scalar>
    bool isIdentity(const Eigen::Tensor<Scalar, 2> &tensor, double threshold = std::numeric_limits<double>::epsilon()) {
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> matrix(tensor.data(), tensor.dimension(0), tensor.dimension(1));
        return matrix.isIdentity(threshold);
    }

    template<typename T, typename Device = Eigen::DefaultDevice>
    auto isIdentity(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors> &expr, double threshold = std::numeric_limits<double>::epsilon(),
                    const Device &device = Device()) {
        auto tensor = tenx::asEval(expr, device);
        static_assert(tensor.rank() == 2 and "isIdentity(): expression must be a tensor of rank 2");
        auto tensorMap = tensor.map();
        using Scalar   = typename decltype(tensorMap)::Scalar;
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> matrix(tensorMap.data(), tensorMap.dimension(0), tensorMap.dimension(1));
        return matrix.isIdentity(threshold);
    }

    template<typename Derived>
    bool hasNaN(const Eigen::EigenBase<Derived> &obj) {
        return obj.derived().hasNaN();
    }

    template<typename Scalar, auto rank>
    bool hasNaN(const Eigen::Tensor<Scalar, rank> &tensor) {
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> vector(tensor.data(), tensor.size());
        return hasNaN(vector);
    }

    template<typename Scalar>
    bool hasNaN(const std::vector<Scalar> &vec) {
        Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> vector(vec.data(), static_cast<long>(vec.size()));
        return hasNaN(vector);
    }

    template<typename Derived>
    auto subtract_phase(Eigen::MatrixBase<Derived> &v) {
        using Scalar = typename Derived::Scalar;
        std::vector<double> angles;
        if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
            for(int i = 0; i < v.cols(); i++) {
                if(v.col(i)(0).imag() == 0.0) {
                    angles.emplace_back(0.0);
                    continue;
                }
                angles.emplace_back(std::arg(v.col(i)(0)));
                Scalar inv_phase     = Scalar(0.0, -1.0) * angles.back();
                Scalar exp_inv_phase = std::exp(inv_phase);
                v.col(i) *= exp_inv_phase;
                v.col(i) = (v.col(i).array().imag().cwiseAbs() > 1e-15).select(v.col(i), v.col(i).real());
            }
        }
        return angles;
    }

    template<typename Scalar, auto rank>
    auto subtract_phase(Eigen::Tensor<Scalar, rank> &tensor) {
        auto map = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(tensor.data(), tensor.size());
        return subtract_phase(map);
    }

    template<typename Derived>
    void add_phase(Eigen::MatrixBase<Derived> &v, std::vector<double> &angles) {
        using Scalar = typename Derived::Scalar;
        if constexpr(std::is_same<Scalar, std::complex<double>>::value) {
            if(v.cols() != angles.size()) { throw std::runtime_error("Mismatch in columns and angles supplied"); }
            for(int i = 0; i < v.cols(); i++) {
                Scalar exp_phase = std::exp(Scalar(0.0, 1.0) * angles[i]);
                v.col(i) *= exp_phase;
                v.col(i) = (v.col(i).array().imag().cwiseAbs() > 1e-15).select(v.col(i), v.col(i).real());
            }
        }
    }

    template<typename Scalar, auto rank>
    void add_phase(Eigen::Tensor<Scalar, rank> &tensor, std::vector<double> &angles) {
        auto map = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(tensor.data(), tensor.size());
        add_phase(map, angles);
    }

    // Compute sparcity
    template<typename Scalar, auto rank>
    double sparcity(const Eigen::Tensor<Scalar, rank> &tensor) {
        auto map = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(tensor.data(), tensor.size());
        return static_cast<double>((map.array() != 0.0).count()) / static_cast<double>(map.size());
    }
    template<typename Derived>
    double sparcity(const Eigen::EigenBase<Derived> &matrix) {
        return static_cast<double>((matrix.derived().array() != 0.0).count()) / static_cast<double>(matrix.derived().size());
    }

}
/*clang-format on */