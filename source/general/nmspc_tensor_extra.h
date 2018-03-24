//
// Created by david on 6/7/17.
//

#ifndef TENSOR_EXTRA_H
#define TENSOR_EXTRA_H

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iterator>
#include <iostream>


/*! \brief **Textra** stands for "Tensor Extra". Provides extra functionality to Eigen::Tensor.*/

/*!
 *  \namespace Textra
 *  This namespace makes shorthand typedef's to Eigen's unsupported Tensor module, and provides handy functions
 *  to interface between `Eigen::Tensor` and `Eigen::Matrix` objects.
 *  The contents of this namespace are quite self-explanatory.
 */


namespace Textra {
    using namespace Eigen;
    using cdouble       = std::complex<double>;

    template<typename Scalar> using MatrixType    = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    template<typename Scalar> using VectorType    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    template<long rank>       using array               = Eigen::array<long, rank>;
    template <typename Scalar, long length>
    using idxlistpair = Eigen::array<Eigen::IndexPair<Scalar>,length>;


    template<typename T, std::size_t N>
    constexpr idxlistpair<long,N> idx (const T (&list1)[N], const T (&list2)[N]){
        //Use numpy-style indexing for contraction. Each list contains a list of indices to be contracted for the respective
        //tensors. This function zips them together into pairs as used in Eigen::Tensor module.
        static_assert(std::is_integral_v<T>);
        Eigen::array<Eigen::IndexPair<long>,N> pairlistOut;
        for(unsigned long i = 0; i < N; i++){
            pairlistOut[i] = Eigen::IndexPair<long>{list1[i], list2[i]};
        }
        return pairlistOut;
    }
    template<typename T>
    struct idx_dim_pair{
        T idxA;
        T idxB;
        T dimB;
    };

    template<std::size_t NB, std::size_t N>
    constexpr idxlistpair<long,N> sortIdx (const Eigen::array<long,NB> &dimensions, const long (&idx_ctrct_A)[N],const long (&idx_ctrct_B)[N]){
        static_assert(std::is_integral_v<long>);
        static_assert(std::is_integral_v<long>);
        Eigen::array<idx_dim_pair<long>,N> idx_dim_pair_list;
        for (unsigned long i = 0; i < N; i++ ){
            idx_dim_pair_list[i] = {idx_ctrct_A[i], idx_ctrct_B[i], dimensions[idx_ctrct_B[i]]};
        }
        std::sort(idx_dim_pair_list.begin(), idx_dim_pair_list.end(), [](const auto& i, const auto& j) { return i.dimB > j.dimB; } );
        idxlistpair<long,N> pairlistOut;
        for(unsigned long i = 0; i< N; i++){
            pairlistOut[i] = Eigen::IndexPair<long>{idx_dim_pair_list[i].idxA, idx_dim_pair_list[i].idxB};
        }
        return pairlistOut;
    }



    using array8        = array<8>;
    using array7        = array<7>;
    using array6        = array<6>;
    using array5        = array<5>;
    using array4        = array<4>;
    using array3        = array<3>;
    using array2        = array<2>;
    using array1        = array<1>;



//
//    //***************************************//
//    //Different views fo rank 1 and 2 tensors//
//    //***************************************//
//
    template <typename Scalar>
    constexpr auto asDiagonal(const Tensor<Scalar,1> &tensor) {
        return tensor.inflate(array1{tensor.size()+1}).reshape(array2{tensor.size(),tensor.size()});
    }

    template <typename Scalar>
    constexpr auto asDiagonalSquared(const Tensor<Scalar,1> &tensor) {
        return tensor.square().inflate(array1{tensor.size()+1}).reshape(array2{tensor.size(), tensor.size()});
    }

    template <typename Scalar>
    constexpr auto asDiagonalInversed(const Tensor<Scalar,1> &tensor) {
        return tensor.inverse().inflate(array1{tensor.size()+1}).reshape(array2{tensor.size(),tensor.size()});
    }

    template<typename Scalar>
    constexpr auto asNormalized(const Tensor<Scalar,1> &tensor) {
        Eigen::Map<const VectorType<Scalar>> map (tensor.data(),tensor.size());
        return Eigen::TensorMap<Tensor<const Scalar,1>>(map.normalized().eval().data(), array1{map.size()});
    }


//
//
//
//    //****************************//
//    //Matrix to tensor conversions//
//    //****************************//
//
//    template<typename Scalar,auto rank>
//    constexpr auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix, const array<rank> &dims){
////        using Scalar =  typename Derived::Scalar;
//        return Eigen::TensorMap<Tensor<const Scalar, rank>>(matrix.derived().eval().data(), dims);
//    }
//
//    template<typename Scalar, typename... Dims>
//    constexpr auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix, Dims... dims){
////        using Scalar =  typename Derived::Scalar;
//        constexpr int rank = sizeof... (Dims);
//        return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.derived().eval().data(), {dims...});
//    }

//    template<typename Derived,auto rank>
//    auto Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, const array<rank> &dims){
//        using Scalar =  typename Derived::Scalar;
//        return Eigen::TensorMap<Tensor<const Scalar, rank>>(matrix.derived().eval().data(), dims);
//    }
//
//    template<typename Derived, typename... Dims>
//    auto Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, Dims... dims){
//        using Scalar =  typename Derived::Scalar;
//        constexpr int rank = sizeof... (Dims);
//        return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.derived().eval().data(), {dims...});
//    }
//
//    template<typename Derived,auto rank>
//    constexpr auto Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, const array<rank> &dims){
//        using Scalar =  typename Derived::Scalar;
//        return Eigen::TensorMap<Tensor<const Scalar, rank>>(matrix.derived().eval().data(), dims);
//    }
//
//    template<typename Derived, typename... Dims>
//    auto Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, Dims... dims){
//        using Scalar =  typename Derived::Scalar;
//        constexpr int rank = sizeof... (Dims);
//        return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.derived().eval().data(), dims...);
//    }
//



    //Detects if Derived is a plain object, like "MatrixXd" or similar.
    //std::decay removes pointer or ref qualifiers if present
    template<typename Derived>
    using is_plainObject_v  = std::is_base_of<Eigen::PlainObjectBase<std::decay_t<Derived> >, std::decay_t<Derived> > ;


    template<typename Derived,auto rank>
    constexpr Tensor<typename Derived::Scalar, rank> Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, const array<rank> &dims){
        if constexpr (is_plainObject_v<Derived>::value) {
            //Return map from raw input.
            return TensorMap<const Tensor<const typename Derived::Scalar, rank>>(matrix.derived().eval().data(), dims);
        }
        else{
            //Create a temporary
            MatrixType<typename Derived::Scalar> matref = matrix;
            return TensorMap<Tensor<typename Derived::Scalar, rank>>(matref.data(), dims);
        }
    }

    template<typename Derived, typename... Dims>
    constexpr Tensor<typename Derived::Scalar, sizeof... (Dims)> Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, Dims... dims) {
        return Matrix_to_Tensor(matrix, array<sizeof...(Dims)>{dims...});
    }



//
//    template<typename Derived,auto rank>
//    constexpr Tensor<typename Derived::Scalar, rank> Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, const array<rank> &dims){
//        //        return TensorMap<Tensor<const Scalar, rank>>(matrix.derived().eval().data(), dims...);
//        return TensorMap<Tensor<typename Derived::Scalar, rank>>(const_cast<typename Derived::Scalar*>(matrix.derived().eval().data()), dims);
//    }
//
//    template<typename Derived, typename... Dims>
//    constexpr Tensor<typename Derived::Scalar, sizeof... (Dims)> Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, Dims... dims){
//       return TensorMap<Tensor<typename Derived::Scalar, sizeof... (Dims)>>(const_cast<typename Derived::Scalar*>(matrix.derived().eval().data()), dims...);
//    }


    template <typename Scalar>
    constexpr auto Matrix_to_Tensor1(const MatrixType<Scalar> &matrix) {
        return Eigen::TensorMap<Tensor<const Scalar,1>>(matrix.data(), array1{matrix.size()});
    }

    template <typename Scalar>
    constexpr auto Matrix_to_Tensor2(const MatrixType<Scalar> &matrix) {
        return Eigen::TensorMap<Tensor<const Scalar,2 >>(matrix.data(), array2{matrix.rows(),matrix.cols()});
    }



//
//    //****************************//
//    //Tensor to matrix conversions//
//    //****************************//
//
//
    template <typename Scalar>
    constexpr MatrixType<Scalar> Tensor2_to_Matrix(const Tensor<Scalar,2> &tensor) {
        return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    }

    template <typename Scalar>
    constexpr MatrixType<Scalar> Tensor1_to_Vector(const Tensor<Scalar,1> &tensor) {
        return Eigen::Map<const VectorType<Scalar>>(tensor.data(), tensor.size());
    }

    template<typename Scalar,auto rank, typename sizeType>
    constexpr MatrixType<Scalar> Tensor_to_Matrix(const Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols){
        return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
    }


    //******************************************************//
    //std::cout overloads for dimension() and array objects //
    //******************************************************//



    template <typename T, int L>
    std::ostream& operator<< (std::ostream& out, const Eigen::DSizes<T,L>& v) {
        if ( !v.empty() ) {
            out << "[ ";
            std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
            out << "]";
        }
        return out;
    }


    template <typename T, int L>
    std::ostream& operator<< (std::ostream& out, const Eigen::array<T,L>& v) {
        if ( !v.empty() ) {
            out << "[ ";
            std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
            out << "]";
        }
        return out;
    }



};




#endif //TENSOR_EXTRA_H
