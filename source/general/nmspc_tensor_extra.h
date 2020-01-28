//
// Created by david on 6/7/17.
//

#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <general/nmspc_type_check.h>


/*! \brief **Textra** stands for "Tensor Extra". Provides extra functionality to Eigen::Tensor.*/

/*!
 *  \namespace Textra
 *  This namespace makes shorthand typedef's to Eigen's unsupported Tensor module, and provides handy functions
 *  to interface between `Eigen::Tensor` and `Eigen::Matrix` objects.
 *  The contents of this namespace is co clear it is self-documenting ;)
 */


namespace Textra {
    using cdouble       = std::complex<double>;


    template<typename Scalar> using MatrixType          = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    template<typename Scalar> using VectorType          = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    template<typename Scalar> using SparseMatrixType    = Eigen::SparseMatrix<Scalar>;
    template<long rank>       using array               = Eigen::array<long, rank>;

    //Shorthand for the list of index pairs.
    template <typename Scalar, long length>
    using idxlistpair = Eigen::array<Eigen::IndexPair<Scalar>,length>;

    constexpr idxlistpair<long,0> idx(){
        Eigen::array<Eigen::IndexPair<long>,0> empty_index_list = {};
        return empty_index_list;
    }


    template<typename T, std::size_t N>
    constexpr idxlistpair<long,N> idx (const T (&list1)[N], const T (&list2)[N]){
        //Use numpy-style indexing for contraction. Each list contains a list of indices to be contracted for the respective
        //tensors. This function zips them together into pairs as used in Eigen::Tensor module. This does not sort the indices in decreasing order.
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
        //When doing contractions, some indices may be larger than others. For performance, you want to
        // contract the largest indices first. This will return a sorted index list in decreasing order.
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



    using array8        = Eigen::array<long,8>;
    using array7        = Eigen::array<long,7>;
    using array6        = Eigen::array<long,6>;
    using array5        = Eigen::array<long,5>;
    using array4        = Eigen::array<long,4>;
    using array3        = Eigen::array<long,3>;
    using array2        = Eigen::array<long,2>;
    using array1        = Eigen::array<long,1>;



//
//    //***************************************//
//    //Different views for rank 1 and 2 tensors//
//    //***************************************//
//

    template <typename Scalar>
    constexpr Eigen::Tensor<Scalar,1> extractDiagonal(const Eigen::Tensor<Scalar,2> &tensor) {
        auto rows = tensor.dimension(0);
        auto cols = tensor.dimension(1);
        assert(tensor.dimension(0) == tensor.dimension(1) and "extractDiagonal expects a square tensor");

        Eigen::Tensor<Scalar,1> diagonals(rows);
        for (auto i = 0; i < rows; i++){
            diagonals(i) = tensor(i,i);
        }
        std::cout << "diagonals: \n" << diagonals;
        return diagonals;
//        return tensor.reshape(array1{rows*cols}).stride(array1{cols+1});
    }


    template <typename Scalar>
    constexpr auto asDiagonal(const Eigen::Tensor<Scalar,1> &tensor) {
        return tensor.inflate(array1{tensor.size()+1}).reshape(array2{tensor.size(),tensor.size()});
    }

    template <typename Scalar>
    constexpr auto asDiagonalSquared(const Eigen::Tensor<Scalar,1> &tensor) {
        return tensor.square().inflate(array1{tensor.size()+1}).reshape(array2{tensor.size(), tensor.size()});
    }

    template <typename Scalar>
    constexpr auto asDiagonalInversed(const Eigen::Tensor<Scalar,1> &tensor) {
        return tensor.inverse().inflate(array1{tensor.size()+1}).reshape(array2{tensor.size(),tensor.size()});
    }

    template <typename Scalar>
    constexpr auto asDiagonalInversed(const Eigen::Tensor<Scalar,2> &tensor) {
        assert(tensor.dimension(0) == tensor.dimension(1) and "Textra::asDiagonalInversed expects a square tensor");
        Eigen::Tensor<Scalar,2> inversed = asDiagonalInversed(extractDiagonal(tensor));
        std::cout << "inversed:\n" << inversed << std::endl;
        return inversed;
    }

    template<typename Scalar>
    constexpr auto asNormalized(const Eigen::Tensor<Scalar,1> &tensor) {
        Eigen::Map<const VectorType<Scalar>> map (tensor.data(),tensor.size());
        return Eigen::TensorMap<Eigen::Tensor<const Scalar,1>>(map.normalized().eval().data(), array1{map.size()});
    }







//    //****************************//
//    //Matrix to tensor conversions//
//    //****************************//


    //Detects if Derived is a plain object, like "MatrixXd" or similar.
    //std::decay removes pointer or ref qualifiers if present
    template<typename Derived>
    using is_plainObject  = std::is_base_of<Eigen::PlainObjectBase<std::decay_t<Derived> >, std::decay_t<Derived> > ;

    template<typename Derived,auto rank>
    constexpr Eigen::Tensor<typename Derived::Scalar, rank> Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, const Eigen::array<long,rank> &dims){
        if constexpr (is_plainObject<Derived>::value) {
            //Return map from raw input.
            return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>(matrix.derived().eval().data(), dims);
        }
        else{
            //Create a temporary
            MatrixType<typename Derived::Scalar> matref = matrix;
            return Eigen::TensorMap<Eigen::Tensor<typename Derived::Scalar, rank>>(matref.data(), dims);
        }
    }

    //Helpful overload
    template<typename Derived, typename... Dims>
    constexpr Eigen::Tensor<typename Derived::Scalar, sizeof... (Dims)> Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, const Dims... dims) {
        return Matrix_to_Tensor(matrix, array<sizeof...(Dims)>{dims...});
    }
    //Helpful overload
    template<typename Derived, auto rank>
    constexpr Eigen::Tensor<typename Derived::Scalar, rank> Matrix_to_Tensor(const Eigen::EigenBase<Derived> &matrix, const Eigen::DSizes<long,rank> &dims) {
        Eigen::array<long,rank> dim_array = dims;
        std::copy(std::begin(dims), std::end(dims), std::begin(dim_array));
        return Matrix_to_Tensor(matrix, dim_array);
    }




    template <typename Derived>
    constexpr auto Matrix_to_Tensor1(const Eigen::EigenBase<Derived> &matrix) {
        return Matrix_to_Tensor(matrix, matrix.size());
    }

    template <typename Derived>
    constexpr auto Matrix_to_Tensor2(const Eigen::EigenBase<Derived> &matrix) {
        return Matrix_to_Tensor(matrix, matrix.rows(),matrix.cols());
    }

    template<typename Derived, auto rank>
    constexpr auto MatrixTensorMap(const Eigen::EigenBase<Derived> &matrix,const Eigen::DSizes<long,rank> &dims){
        if constexpr (is_plainObject<Derived>::value) {
            return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>(matrix.derived().data(), dims);
        }
        else{
            return Matrix_to_Tensor(matrix,dims);
        }
    }

    template<typename Derived, typename... Dims>
    constexpr auto MatrixTensorMap(const Eigen::EigenBase<Derived> &matrix,const Dims... dims){
        return MatrixTensorMap(matrix, Eigen::DSizes<long,sizeof...(Dims)>{dims...});
    }

    template<typename Derived>
    constexpr auto MatrixTensorMap(const Eigen::EigenBase<Derived> &matrix){
        if constexpr(Derived::ColsAtCompileTime ==  1){
            return MatrixTensorMap(matrix, matrix.size());
        }else{
            return MatrixTensorMap(matrix, matrix.rows(),matrix.cols());
        }
    }


//
//    //****************************//
//    //Tensor to matrix conversions//
//    //****************************//
//


    template <typename Scalar>
    constexpr MatrixType<Scalar> Tensor2_to_Matrix(const Eigen::Tensor<Scalar,2> &tensor) {
        return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    }

    template <typename Scalar>
    constexpr MatrixType<Scalar> Tensor1_to_Vector(const Eigen::Tensor<Scalar,1> &tensor) {
        return Eigen::Map<const VectorType<Scalar>>(tensor.data(), tensor.size());
    }

    template<typename Scalar,auto rank, typename sizeType>
    constexpr MatrixType<Scalar> Tensor_to_Matrix(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols){
        return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
    }

    template <typename Scalar>
    constexpr SparseMatrixType<Scalar> Tensor2_to_SparseMatrix(const Eigen::Tensor<Scalar,2> &tensor, double prune_threshold = 1e-15) {
        return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1)).sparseView().pruned(prune_threshold);
    }


    template<typename Scalar,auto rank, typename sizeType>
    constexpr auto TensorMatrixMap(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols){
        return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
    }

    template<typename Scalar>
    constexpr auto TensorMatrixMap(const Eigen::Tensor<Scalar,2> &tensor){
        return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), tensor.dimension(0),tensor.dimension(1));
    }

    template<typename Scalar>
    constexpr auto TensorVectorMap(const Eigen::Tensor<Scalar,2> &tensor){
        return Eigen::Map<const VectorType<Scalar>> (tensor.data(), tensor.size());
    }






    //************************//
    // change storage layout //
    //************************//
    template<typename Scalar,auto rank>
    Eigen::Tensor<Scalar,rank, Eigen::RowMajor> to_RowMajor(const Eigen::Tensor<Scalar,rank, Eigen::ColMajor> tensor){
        std::array<long,rank> neworder;
        std::iota(std::begin(neworder), std::end(neworder), 0);
        std::reverse(neworder.data(), neworder.data()+neworder.size());
        return tensor.swap_layout().shuffle(neworder);
    }


    template<typename Derived>
    Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> to_RowMajor(const Eigen::MatrixBase<Derived> &matrix){
        if(matrix.IsRowMajor) {return matrix;}
        Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor> matrowmajor = matrix;
        return matrowmajor;
    }

    template<typename Scalar,auto rank>
    Eigen::Tensor<Scalar,rank, Eigen::ColMajor> to_ColMajor(const Eigen::Tensor<Scalar,rank, Eigen::RowMajor> tensor){
        std::array<long,rank> neworder;
        std::iota(std::begin(neworder), std::end(neworder), 0);
        std::reverse(neworder.data(), neworder.data()+neworder.size());
        return tensor.swap_layout().shuffle(neworder);
    }

    template<typename Derived>
    Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::ColMajor> to_ColMajor(const Eigen::MatrixBase<Derived> &matrix){
        if(not matrix.IsRowMajor) {return matrix;}
        Eigen::Matrix<typename Derived::Scalar,Eigen::Dynamic,Eigen::Dynamic, Eigen::ColMajor> matrowmajor = matrix;
        return matrowmajor;
    }


//******************************************************//
// Tests for real/complex //
//******************************************************//



    template<typename Derived>
    bool isReal(const Eigen::EigenBase<Derived> &obj,[[maybe_unused]]const std::string &name = "", double threshold = 1e-14) {
        using Scalar = typename Derived::Scalar;
        if constexpr (TypeCheck::is_specialization<Scalar, std::complex>::value){
            auto imag_sum = obj.derived().imag().cwiseAbs().sum();
//            std::cout <<"imag sum " << name << " :" << std::fixed << std::setprecision(16)<< imag_sum << std::endl;
            return imag_sum < threshold;
        }else{
            return true;
        }
    }


    template<typename Scalar, auto rank>
    bool isReal(const Eigen::Tensor<Scalar,rank> &tensor, const std::string &name = "", double threshold = 1e-14) {
        Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> vector (tensor.data(),tensor.size());
        return isReal(vector, name, threshold);
    }

    template<typename Derived>
    bool hasNaN(const Eigen::EigenBase<Derived> &obj,[[maybe_unused]]const std::string &name = "") {
        return obj.derived().hasNaN();
    }


    template<typename Scalar, auto rank>
    bool hasNaN(const Eigen::Tensor<Scalar,rank> &tensor, const std::string &name = "") {
        Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> vector (tensor.data(),tensor.size());
        return hasNaN(vector, name);
    }


    template <typename Derived>
    auto subtract_phase(Eigen::MatrixBase<Derived> &v){
        using Scalar = typename Derived::Scalar;
        std::vector<double> angles;
        if constexpr (std::is_same<Scalar, std::complex<double>>::value) {
            for (int i = 0; i < v.cols(); i++){
                if (v.col(i)(0).imag() == 0.0){
                    angles.emplace_back(0.0);
                    continue;
                }
                angles.emplace_back(std::arg(v.col(i)(0)));
                Scalar inv_phase = Scalar(0.0,-1.0) * angles.back();
                Scalar exp_inv_phase = std::exp(inv_phase);
//                    std::cout << std::setprecision(20) << std::fixed << "inv_phase    : " << inv_phase << std::endl;
//                    std::cout << std::setprecision(20) << std::fixed << "exp_inv_phase: " << exp_inv_phase << std::endl;
                v.col(i) *= exp_inv_phase;
                v.col(i) = (v.col(i).array().imag().cwiseAbs() > 1e-15  ).select(v.col(i), v.col(i).real());
            }
        }
        return angles;
    }



    template<typename Scalar, auto rank>
    auto subtract_phase(Eigen::Tensor<Scalar,rank> &tensor){
        auto map = Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(tensor.data(),tensor.size());
        return subtract_phase(map);
    }


    template <typename Derived>
    void add_phase(Eigen::MatrixBase<Derived> &v, std::vector<double> & angles){
        using Scalar = typename Derived::Scalar;
        if constexpr (std::is_same<Scalar, std::complex<double>>::value) {
            if (v.cols() != angles.size() ){throw std::runtime_error("Mismatch in columns and angles supplied");}
            for (int i = 0; i < v.cols(); i++){
                Scalar exp_phase = std::exp(Scalar(0.0,1.0) * angles[i]);
                v.col(i) *= exp_phase;
                v.col(i) = (v.col(i).array().imag().cwiseAbs() > 1e-15  ).select(v.col(i), v.col(i).real());
            }
        }
    }


    template<typename Scalar, auto rank>
    void add_phase(Eigen::Tensor<Scalar,rank> &tensor, std::vector<double> &angles){
        auto map = Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(tensor.data(),tensor.size());
        add_phase(map,angles);
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


/*! \brief Prints the content of a list nicely */
    template<typename T>
    std::ostream &operator<<(std::ostream &out, const std::list<T> &v) {
        if (!v.empty()) {
            out << "[ ";
            std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
            out << "]";
        }
        return out;
    }

/*! \brief Prints the content of a list nicely */
    template<typename T>
    std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
        if (!v.empty()) {
            out << "[ ";
            std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
            out << "]";
        }
        return out;
    }

}


