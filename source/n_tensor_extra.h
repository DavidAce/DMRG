//
// Created by david on 6/7/17.
//

#ifndef TENSOR_EXTRA_H
#define TENSOR_EXTRA_H
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
//    using Scalar        = double;
    using cdouble       = std::complex<double>;
    using idx2          = Eigen::IndexPair<long>;

    template<typename Scalar = double> using MatrixType    = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    template<typename Scalar = double> using VectorType    = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    template<long rank> using idxlist                                 = Eigen::array<idx2, rank>;
    template<long rank, typename Scalar = double> using Tensor        = Eigen::Tensor<Scalar,rank,Eigen::ColMajor>;
    template<long rank, typename Scalar = double> using const_Tensor  = Eigen::Tensor<const Scalar,rank,Eigen::ColMajor>;
    template<long rank> using array                                   = Eigen::array<long, rank>;



    using Tensor4d       = Tensor<4,double>;
    using Tensor3d       = Tensor<3,double>;
    using Tensor2d       = Tensor<2,double>;
    using Tensor1d       = Tensor<1,double>;
    using Tensor0d       = Tensor<0,double>;
    using Tensor4cd      = Tensor<4,cdouble>;
    using Tensor3cd      = Tensor<3,cdouble>;
    using Tensor2cd      = Tensor<2,cdouble>;
    using Tensor1cd      = Tensor<1,cdouble>;
    using Tensor0cd      = Tensor<0,cdouble>;
    using const_Tensor4d = Tensor<4,const double>;
    using const_Tensor3d = Tensor<3,const double>;
    using const_Tensor2d = Tensor<2,const double>;
    using const_Tensor1d = Tensor<1,const double>;

    using array4        = array<4>;
    using array3        = array<3>;
    using array2        = array<2>;
    using array1        = array<1>;

    using idxlist4      = idxlist<4>;
    using idxlist3      = idxlist<3>;
    using idxlist2      = idxlist<2>;
    using idxlist1      = idxlist<1>;

    template <typename Scalar = double>
    inline Tensor<2,Scalar> asDiagonal(const Tensor<1, Scalar> &tensor) {
        MatrixType<Scalar> mat = Eigen::Map<const VectorType<Scalar>>(tensor.data(), tensor.size()).asDiagonal();
        return Eigen::TensorMap<Tensor<2,Scalar>>(mat.data(), mat.rows(), mat.cols());
    }

    template <typename Scalar = double>
    inline Tensor<2, Scalar> asInverseDiagonal(const Tensor<1, Scalar> &tensor) {
        return asDiagonal((Tensor<1, Scalar>)tensor.inverse());
    }

    inline Tensor2d Tensor2d_inverse(const Tensor2d &tensor) {
        Eigen::Map<const Eigen::MatrixXd> map =  Eigen::Map<const Eigen::MatrixXd>(tensor.data(),tensor.dimension(0), tensor.dimension(1));
        return Eigen::TensorMap<const_Tensor2d>(map.inverse().eval().data(), array2{map.rows(),map.cols()});
    }

    inline Tensor1d Tensor1d_normalize(const Tensor1d &tensor) {
        Eigen::Map<const Eigen::VectorXd> map =  Eigen::Map<const Eigen::VectorXd>(tensor.data(),tensor.dimension(0));
        return Eigen::TensorMap<const_Tensor1d>(map.normalized().eval().data(), array1{map.size()});
    }


    //Matrix to tensor conversions
    template<long rank, typename Scalar>
    Tensor<rank, Scalar> Matrix_to_Tensor(const MatrixType<Scalar> &matrix, array<rank> dims){
        return Eigen::TensorMap<const_Tensor<rank, Scalar>>(matrix.data(), dims);
    }



    template <typename Scalar>
    inline Tensor<1, Scalar> Matrix_to_Tensor1(const MatrixType<Scalar> &matrix) {
        return Eigen::TensorMap<const_Tensor<1,Scalar>>(matrix.data(), array1{matrix.size()});
    }

    inline Tensor1d MatrixXd_to_Tensor1d(const Eigen::MatrixXd &matrix) {
        return Eigen::TensorMap<const_Tensor1d>(matrix.data(), array1{matrix.size()});
    }

    template <typename Scalar>
    inline Tensor<2,Scalar> Matrix_to_Tensor2(const MatrixType<Scalar> &matrix) {
        return Eigen::TensorMap<const_Tensor<2,Scalar>>(matrix.data(), array2{matrix.rows(),matrix.cols()});
    }

    inline Tensor2d MatrixXd_to_Tensor2d(const Eigen::MatrixXd &matrix) {
        return Eigen::TensorMap<const_Tensor2d>(matrix.data(), array2{matrix.rows(),matrix.cols()});
    }



    //Tensor to matrix conversions


    template <typename Scalar>
    inline MatrixType<Scalar> tensor2_to_matrix(const Tensor<2,Scalar> &tensor) {
        return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    }

    inline Eigen::MatrixXd Tensor2d_to_MatrixXd(const Tensor2d &tensor){
        return Eigen::Map<const Eigen::MatrixXd>(tensor.data(),tensor.dimension(0), tensor.dimension(1));
    }


    template <typename Scalar>
    inline VectorType<Scalar> Tensor1_to_Vector(const Tensor<1, Scalar> &tensor) {
        return Eigen::Map<const VectorType<Scalar>>(tensor.data(), tensor.size());
    }

    inline Eigen::VectorXd Tensor1d_to_VectorXd(const Tensor1d &tensor){
        return Eigen::Map<const Eigen::VectorXd>(tensor.data(), tensor.size());
    }

    template<long rank, typename sizeType>
    inline Eigen::MatrixXd Tensor_to_MatrixXd(const Tensor<rank,double> &tensor, sizeType rows, sizeType cols){
        Tensor<2,double> tensor2d  = tensor.reshape(array2{rows,cols});
        return Eigen::Map<const Eigen::MatrixXd> (tensor2d.data(), rows ,cols );
    }




//    template<long rank, typename ...Args>
//    Tensor<rank, double> matrixXd_to_tensor(Args && ...args){
//        matrix_to_tensor<rank,double>(std::forward<Args>(args)...);
//    }



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
