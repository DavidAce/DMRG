//
// Created by david on 6/7/17.
//

#include <n_tensor_extra.h>

namespace Textra {
    Tensor3 reshape(MatrixType matrix, array3 dims) {
        return Eigen::TensorMap<Tensor3>(matrix.data(), dims);
    }


    Tensor2 asDiagonal(const Eigen::TensorRef<const Tensor1> tensor) {
        MatrixType mat = Eigen::Map<const VectorType>(tensor.data(), tensor.size()).asDiagonal();
        return Eigen::TensorMap<Tensor2>(mat.data(), mat.rows(), mat.cols());
    }

    Tensor2 asInverseDiagonal(const Eigen::TensorRef<const Tensor1> tensor) {
        return asDiagonal((Tensor1)tensor.inverse());
    }

    Tensor1 matrix_to_tensor1(const Eigen::Ref<const MatrixType> matrix) {
        return Eigen::TensorMap<const_Tensor1>(matrix.data(), array1{matrix.size()});
    }

    Tensor2 matrix_to_tensor2(const Eigen::Ref<const MatrixType> matrix, array2 dims) {
        return Eigen::TensorMap<const_Tensor2>(matrix.data(), dims);
    }

    Tensor2 matrix_to_tensor2(const Eigen::Ref<const MatrixType> matrix) {
        return Eigen::TensorMap<const_Tensor2>(matrix.data(), array2{matrix.rows(),matrix.cols()});
    }


    Tensor3 matrix_to_tensor3(const Eigen::Ref<const MatrixType> matrix, array3 dims) {
        return Eigen::TensorMap<const_Tensor3>(matrix.data(), dims);
    }


    Tensor4 matrix_to_tensor4(const Eigen::Ref<const MatrixType> matrix, array4 dims) {
        return Eigen::TensorMap<const_Tensor4>(matrix.data(), dims);
    }

    MatrixType tensor2_to_matrix(const Tensor2 &tensor) {
        return Eigen::Map<const MatrixType>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    }

}