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
    template<long rank> using array               = Eigen::array<long, rank>;
    template <typename Scalar, long length>
    using idxlistpair = Eigen::array<Eigen::IndexPair<Scalar>,length>;
    template <long length>
    inline idxlistpair<long,length> idx (std::initializer_list<long> list1,std::initializer_list<long> list2){
        //Use numpy-style indexing for contraction. Each list contains a list of indices to be contracted for the respective
        //tensors. This function zips them together into pairs as used in Eigen::Tensor module.
        Eigen::array<Eigen::IndexPair<long>,length> pairlistOut;
        auto it1 = list1.begin();
        auto it2 = list2.begin();
        for(int i = 0; i< length; i++){
            pairlistOut[i] = Eigen::IndexPair<long>{*it1++, *it2++};
        }
        return pairlistOut;
    }
    struct idx_dim_pair{
        long idxA;
        long idxB;
        long dimB;
    };
    template< long length, long rankB>
    inline idxlistpair<long,length> __attribute__((hot)) sortIdx (const array<rankB>  dimsB, const std::initializer_list<long> listA, const std::initializer_list<long> listB){
        Eigen::array<idx_dim_pair,length> idx_dim_pair_list;
        auto itlistA = listA.begin();
        auto itlistB = listB.begin();
        for (int i = 0; i < length; i++ ){
            idx_dim_pair_list[i] = {*itlistA++, *itlistB, dimsB[*itlistB++]};
        }
        std::sort(idx_dim_pair_list.begin(), idx_dim_pair_list.end(), [](const auto& i, const auto& j) { return i.dimB > j.dimB; } );
        idxlistpair<long,length> pairlistOut;
        for(int i = 0; i< length; i++){
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
    inline Tensor<Scalar,2> __attribute__((always_inline)) asDiagonal(const Tensor<Scalar,1> &tensor) {
        return tensor.inflate(array1{tensor.size()+1}).reshape(array2{tensor.size(),tensor.size()});
    }

    template <typename Scalar>
    inline Tensor<Scalar,2> __attribute__((always_inline)) asDiagonalSquared(const Tensor<Scalar,1> &tensor) {
        return tensor.square().inflate(array1{tensor.size()+1}).reshape(array2{tensor.size(), tensor.size()});
    }

    template <typename Scalar>
    inline Tensor<Scalar,2> __attribute__((always_inline)) asDiagonalInversed(const Tensor<Scalar,1> &tensor) {
        return tensor.inverse().inflate(array1{tensor.size()+1}).reshape(array2{tensor.size(),tensor.size()});
    }

    template<typename Scalar>
    inline Tensor<Scalar,1>  __attribute__((always_inline)) asNormalized(const Tensor<Scalar,1> &tensor) {
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

    template<typename Scalar,int rank>
    inline Eigen::Tensor<Scalar,rank> __attribute__((always_inline)) Matrix_to_Tensor(const MatrixType<Scalar> &matrix, const array<rank> &dims){
        return Eigen::TensorMap<Tensor<const Scalar, rank>>(matrix.data(), dims);
    }


    template <typename Scalar>
    inline Eigen::Tensor<Scalar,1> __attribute__((always_inline)) Matrix_to_Tensor1(const MatrixType<Scalar> &matrix) {
        return Eigen::TensorMap<Tensor<const Scalar,1>>(matrix.data(), array1{matrix.size()});
    }

    template <typename Scalar>
    inline Eigen::Tensor<Scalar,2> __attribute__((always_inline)) Matrix_to_Tensor2(const MatrixType<Scalar> &matrix) {
        return Eigen::TensorMap<Tensor<const Scalar,2 >>(matrix.data(), array2{matrix.rows(),matrix.cols()});
    }



//
//    //****************************//
//    //Tensor to matrix conversions//
//    //****************************//
//
//
    template <typename Scalar>
    inline MatrixType<Scalar> __attribute__((always_inline)) Tensor2_to_Matrix(const Tensor<Scalar,2> &tensor) {
        return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    }

    template <typename Scalar>
    inline VectorType<Scalar> __attribute__((always_inline)) Tensor1_to_Vector(const Tensor<Scalar,1> &tensor) {
        return Eigen::Map<const VectorType<Scalar>>(tensor.data(), tensor.size());
    }

    template<typename Scalar,long rank, typename sizeType>
    inline MatrixType<Scalar> __attribute__((always_inline)) Tensor_to_Matrix(const Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols){
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
