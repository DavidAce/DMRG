//
// Created by david on 6/7/17.
//

#include <n_tensor_extra.h>

namespace Textra {

    Tensor<2> asDiagonal(const Eigen::TensorRef<const Tensor<1>> tensor) {
        MatrixType mat = Eigen::Map<const VectorType>(tensor.data(), tensor.size()).asDiagonal();
        return Eigen::TensorMap<Tensor<2>>(mat.data(), mat.rows(), mat.cols());
    }

    Tensor<2> asInverseDiagonal(const Eigen::TensorRef<const Tensor<1>> tensor) {
        return asDiagonal((Tensor<1>)tensor.inverse());
    }

    Tensor<1> matrix_to_tensor1(const Eigen::Ref<const MatrixType> matrix) {
        return Eigen::TensorMap<const_Tensor<1>>(matrix.data(), array<1>{matrix.size()});
    }
    Tensor<2> matrix_to_tensor2(const Eigen::Ref<const MatrixType> matrix) {
        return Eigen::TensorMap<const_Tensor<2>>(matrix.data(), array<2>{matrix.rows(),matrix.cols()});
    }

//    Tensor<2> matrix_to_tensor2(const Eigen::Ref<const MatrixType> matrix, array<2> dims) {
//        return Eigen::TensorMap<const_Tensor<2>>(matrix.data(), dims);
//    }
//
//
//
//
//    Tensor3 matrix_to_tensor3(const Eigen::Ref<const MatrixType> matrix, array3 dims) {
//        return Eigen::TensorMap<const_Tensor3>(matrix.data(), dims);
//    }
//
//
//    Tensor4 matrix_to_tensor4(const Eigen::Ref<const MatrixType> matrix, array4 dims) {
//        return Eigen::TensorMap<const_Tensor4>(matrix.data(), dims);
//    }

//    MatrixType tensor2_to_matrix(const Tensor2 &tensor) {
//        return Eigen::Map<const MatrixType>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
//    }
    MatrixType tensor2_to_matrix(const Eigen::TensorRef<const Tensor<2>> tensor) {
        return Eigen::Map<const MatrixType>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    }
}