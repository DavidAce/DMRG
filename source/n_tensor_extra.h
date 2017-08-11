//
// Created by david on 6/7/17.
//

#ifndef TEBD_EIGEN_N_TENSOR_EXTRA_H
#define TEBD_EIGEN_N_TENSOR_EXTRA_H
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iterator>
#include <vector>
#include <iostream>
/*! \brief **Textra** stands for "Tensor Extra". Provides extra functionality to Eigen::Tensor.*/

/*!
 *  \namespace Textra
 *  This namespace makes shorthand typedef's to Eigen's unsupported Tensor module, and provides handy functions
 *  to interface between `Eigen::Tensor` and `Eigen::Matrix` objects.
 *  The contents of this namespace are quite self-explanatory.
 */

namespace Textra {
    using Scalar        = double;

    using idx2          = Eigen::IndexPair<long>;
    using MatrixType    = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>;
    using VectorType    = Eigen::Matrix<Scalar,Eigen::Dynamic, 1, Eigen::ColMajor>;

    template<long rank> using idxlist       = Eigen::array<idx2, rank>;
    template<long rank> using Tensor        = Eigen::Tensor<Scalar,rank,Eigen::ColMajor>;
    template<long rank> using const_Tensor  = Eigen::Tensor<const Scalar,rank,Eigen::ColMajor>;
    template<long rank> using array         = Eigen::array<long, rank>;



    using Tensor4       = Eigen::Tensor<Scalar,4,Eigen::ColMajor>;
    using Tensor3       = Eigen::Tensor<Scalar,3,Eigen::ColMajor>;
    using Tensor2       = Eigen::Tensor<Scalar,2,Eigen::ColMajor>;
    using Tensor1       = Eigen::Tensor<Scalar,1,Eigen::ColMajor>;
    using Tensor0       = Eigen::Tensor<Scalar,0,Eigen::ColMajor>;
    using const_Tensor4 = Eigen::Tensor<const Scalar,4,Eigen::ColMajor>;
    using const_Tensor3 = Eigen::Tensor<const Scalar,3,Eigen::ColMajor>;
    using const_Tensor2 = Eigen::Tensor<const Scalar,2,Eigen::ColMajor>;
    using const_Tensor1 = Eigen::Tensor<const Scalar,1,Eigen::ColMajor>;

    using array4        = Eigen::array<long, 4>;
    using array3        = Eigen::array<long, 3>;
    using array2        = Eigen::array<long, 2>;
    using array1        = Eigen::array<long, 1>;

    using idxlist4      = Eigen::array<idx2, 4>;
    using idxlist3      = Eigen::array<idx2, 3>;
    using idxlist2      = Eigen::array<idx2, 2>;
    using idxlist1      = Eigen::array<idx2, 1>;


//    extern Tensor<2> asDiagonal        (const Eigen::TensorRef<const Tensor<1>> tensor);
//    extern Tensor<2> asInverseDiagonal (const Eigen::TensorRef<const Tensor<1>> tensor);
//    extern Tensor<1> matrix_to_tensor1(const Eigen::Ref<const MatrixType> matrix);
//    extern Tensor<2> matrix_to_tensor2(const Eigen::Ref<const MatrixType> matrix);
//    extern MatrixType tensor2_to_matrix(const Eigen::TensorRef<const Tensor<2>> tensor);

    inline Tensor<2> asDiagonal(const Eigen::TensorRef<const Tensor<1>> tensor) {
        MatrixType mat = Eigen::Map<const VectorType>(tensor.data(), tensor.size()).asDiagonal();
        return Eigen::TensorMap<Tensor<2>>(mat.data(), mat.rows(), mat.cols());
    }

    inline Tensor<2> asInverseDiagonal(const Eigen::TensorRef<const Tensor<1>> tensor) {
        return asDiagonal((Tensor<1>)tensor.inverse());
    }

    inline Tensor<1> matrix_to_tensor1(const Eigen::Ref<const MatrixType> matrix) {
        return Eigen::TensorMap<const_Tensor<1>>(matrix.data(), array<1>{matrix.size()});
    }


    inline Tensor<2> matrix_to_tensor2(const Eigen::Ref<const MatrixType> matrix) {
        return Eigen::TensorMap<const_Tensor<2>>(matrix.data(), array<2>{matrix.rows(),matrix.cols()});
    }

    inline MatrixType tensor2_to_matrix(const Eigen::TensorRef<const Tensor<2>> tensor) {
        return Eigen::Map<const MatrixType>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    }

    template<long rank>
    Tensor<rank> matrix_to_tensor(const Eigen::Ref<const MatrixType> matrix, array<rank> dims){
        return Eigen::TensorMap<const_Tensor<rank>>(matrix.data(), dims);
    }




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


#endif //TEBD_EIGEN_N_TENSOR_EXTRA_H
