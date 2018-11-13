//
// Created by david on 2018-10-29.
//

#include "matrix_recast.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iterator>

/*! \brief Prints the content of a vector nicely */
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    if (!v.empty()) {
//        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, "  "));
        out << '\n';
    }
    return out;
}

template<typename Scalar>
matrix_recast<Scalar>::matrix_recast(const Scalar *matrix_ptr_, int L_):matrix_ptr(matrix_ptr_), L(L_){
    recheck_all();
}

template<typename Scalar>
void matrix_recast<Scalar>::prune(double threshold){
    matrix_pruned.resize(L*L);
    Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> matmap   (matrix_ptr,L,L);
    Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> matpruned(matrix_pruned.data(),L,L);
    matpruned = (matmap.array().cwiseAbs() < threshold).select(Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>::Zero(L,L), matmap);
    recheck_all();
    pruned = true;
}


template<typename Scalar>
void matrix_recast<Scalar>:: recheck_all(){
    check_if_sparse();
    check_if_hermitian();
    check_if_real();
}




template<typename Scalar>
void matrix_recast<Scalar>::check_if_real() {
    Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> matrix (matrix_ptr,L,L);
    if constexpr (std::is_same<Scalar, double>::value){isReal = true;}
    else {isReal = matrix.imag().isZero(1e-14);}
}


template<typename Scalar>
void matrix_recast<Scalar>::check_if_hermitian() {
    Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> matrix (matrix_ptr,L,L);
    isHermitian = matrix.isApprox(matrix.adjoint(), 1e-14);
}


template<typename Scalar>
void matrix_recast<Scalar>::check_if_sparse() {
    Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> matrix (matrix_ptr,L,L);
    sparcity = (matrix.array().cwiseAbs() > 1e-14 )
                       .select(Eigen::MatrixXd::Ones(L,L),0).sum() / matrix.size();
    isSparse =  sparcity < 0.1;
}

template<typename Scalar>
DenseMatrixProduct<double> matrix_recast<Scalar>::get_as_real_dense() {
    if constexpr(std::is_same<Scalar, double>::value){
        if (pruned){return DenseMatrixProduct<double>(matrix_pruned.data(),L);}
        else       {return DenseMatrixProduct<double>(matrix_ptr,L);}
    }else{
//        assert(isReal and "ERROR: The given matrix has a nonzero imaginary part. Can't convert to real.");
        if (not isReal){
            double sum = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_ptr,L,L).imag().array().cwiseAbs().sum();
            std::cerr << "WARNING: The given matrix has a nonzero imaginary part, yet converting to real. Imag sum: " << sum << std::endl;
        }
        Eigen::MatrixXd matrix_recast;
        if (pruned){matrix_recast = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_pruned.data(),L,L).real();}
        else       {matrix_recast = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_ptr,L,L).real();}
        return DenseMatrixProduct<double>(matrix_recast.data(),L);


    }
}

template<typename Scalar>
DenseMatrixProduct<std::complex<double>> matrix_recast<Scalar>::get_as_cplx_dense() {
    if constexpr(std::is_same<Scalar, double>::value){
        Eigen::MatrixXcd matrix_recast(L,L);
        matrix_recast.setZero();
        if (pruned){matrix_recast.real() = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_pruned.data(),L,L);}
        else       {matrix_recast.real() = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_ptr,L,L);}
        return DenseMatrixProduct<std::complex<double>>(matrix_recast.data(),L);
    }else{
        return DenseMatrixProduct<std::complex<double>>(matrix_ptr,L);
    }
}

template<typename Scalar>
SparseMatrixProduct<double> matrix_recast<Scalar>::get_as_real_sparse() {
    if constexpr(std::is_same<Scalar, double>::value){
        if(pruned){return SparseMatrixProduct<double>(matrix_pruned.data(),L);}
        else      {return SparseMatrixProduct<double>(matrix_ptr,L);}
    }else{
        if (not isReal){
            double sum = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_ptr,L,L).imag().array().cwiseAbs().sum();
            std::cerr << "WARNING: The given matrix has a nonzero imaginary part, yet converting to real. Imag sum: " << sum << std::endl;
        }
        Eigen::MatrixXd matrix_recast;
        if(pruned){matrix_recast = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_pruned.data(),L,L).real();}
        else      {matrix_recast = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_ptr,L,L).real();}
        return SparseMatrixProduct<double>(matrix_recast);
    }
}

template<typename Scalar>
SparseMatrixProduct<std::complex<double>> matrix_recast<Scalar>::get_as_cplx_sparse() {
    if constexpr(std::is_same<Scalar, double>::value){
        Eigen::MatrixXcd matrix_recast(L,L);
        matrix_recast.setZero();
        if(pruned){matrix_recast.real() = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_pruned.data(),L,L);}
        else      {matrix_recast.real() = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_ptr,L,L);}
        return SparseMatrixProduct<std::complex<double>>(matrix_recast);
    }else{
        if(pruned){return SparseMatrixProduct<std::complex<double>>(matrix_pruned.data(),L);}
        else      {return SparseMatrixProduct<std::complex<double>>(matrix_ptr,L);}
    }
}





//
//template<typename Scalar>
//void matrix_recast<Scalar>::convert_to_real_dense(){
//    // READ THIS TO LEARN MORE http://atantet.github.io/ATSuite_cpp/atspectrum_8hpp_source.html
//
//    matrix_real_dense.clear();
//    matrix_cplx_dense.clear();
////    matrix_real_sparse.clear();
////    matrix_cplx_sparse.clear();
//
//
////    matrix_real_dense.L    = L;
////    matrix_real_dense.N    = L*L;
//    if constexpr(std::is_same<Scalar, double>::value){
//        matrix_real_dense = DenseMatrixProduct<double>(matrix_ptr,L);
//    }else{
//        Eigen::MatrixXd matrix_recast = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_ptr,L,L).real();
//        matrix_real_dense = DenseMatrixProduct<double>(matrix_recast.data(),L);
//    }
//
//}
//
//
//template<typename Scalar>
//void matrix_recast<Scalar>::convert_to_cplx_dense(){
//    // READ THIS TO LEARN MORE http://atantet.github.io/ATSuite_cpp/atspectrum_8hpp_source.html
//
//    matrix_real_dense.clear();
//    matrix_cplx_dense.clear();
////    matrix_real_sparse.clear();
////    matrix_cplx_sparse.clear();
//
//
////    matrix_cplx_dense.L    = L;
////    matrix_cplx_dense.N    = L*L;
//    if constexpr(std::is_same<Scalar, double>::value){
//        Eigen::MatrixXcd matrix_recast = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> (matrix_ptr,L,L);
//        matrix_cplx_dense =  DenseMatrixProduct<std::complex<double>>(matrix_recast.data(),L);
//    }else{
//        matrix_cplx_dense = DenseMatrixProduct<std::complex<double>>(matrix_ptr,L);
//    }
//}


//template<typename Scalar>
//void matrix_recast<Scalar>::convert_to_real_sparse(){
    // READ THIS TO LEARN MORE http://atantet.github.io/ATSuite_cpp/atspectrum_8hpp_source.html

//    matrix_real_dense.clear();
//    matrix_cplx_dense.clear();
//    matrix_real_sparse.clear();
//    matrix_cplx_sparse.clear();
//
//    assert(isReal and "Matrix is not real!");
//    Eigen::SparseMatrix<double> matrix_recast = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(matrix_ptr,L,L).real().sparseView();
//    matrix_recast.makeCompressed();
//    matrix_real_sparse.nnz = matrix_recast.nonZeros();
//    matrix_real_sparse.L   = L;
//    matrix_real_sparse.N   = L*L;
//    matrix_real_sparse.irow = std::vector<int>   (matrix_recast.innerIndexPtr(), matrix_recast.innerIndexPtr() + matrix_recast.nonZeros() );
//    matrix_real_sparse.pcol = std::vector<int>   (matrix_recast.outerIndexPtr(), matrix_recast.outerIndexPtr() + matrix_recast.outerSize() +1);
//    matrix_real_sparse.vals = std::vector<double>(matrix_recast.valuePtr()     , matrix_recast.valuePtr()      + matrix_recast.nonZeros());
//

//}


//template<typename Scalar>
//void matrix_recast<Scalar>::convert_to_cplx_sparse(){
    // READ THIS TO LEARN MORE http://atantet.github.io/ATSuite_cpp/atspectrum_8hpp_source.html



//    matrix_real_dense.clear();
//    matrix_cplx_dense.clear();
//    matrix_real_sparse.clear();
//    matrix_cplx_sparse.clear();
//    assert(not isReal and "Matrix is not cplx!");
//    Eigen::SparseMatrix<std::complex<double>> matrix_recast = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>>(matrix_ptr,L,L).template cast<std::complex<double>>().sparseView();
//    matrix_recast.makeCompressed();
//    matrix_recast.finalize();
//    std::cout << "Matrix compressed: \n" << matrix_recast << std::endl;
//
//    matrix_cplx_sparse.nnz = matrix_recast.nonZeros();
//    matrix_cplx_sparse.L   = L;
//    matrix_cplx_sparse.N   = L*L;
//    matrix_cplx_sparse.irow = std::vector<int>   (matrix_recast.innerIndexPtr(), matrix_recast.innerIndexPtr() + matrix_recast.nonZeros() );
//    matrix_cplx_sparse.pcol = std::vector<int>   (matrix_recast.outerIndexPtr(), matrix_recast.outerIndexPtr() + matrix_recast.outerSize() +1);
//    matrix_cplx_sparse.vals = std::vector<std::complex<double>>(matrix_recast.valuePtr()       , matrix_recast.valuePtr()       + matrix_recast.nonZeros());
//
//    std::cout << "nonzeros : " << matrix_recast.nonZeros() << std::endl;
//    std::cout << "innerSize: " << matrix_recast.innerSize() << std::endl;
//    std::cout << "outerSize: " << matrix_recast.outerSize() << std::endl;
//    std::cout << "irow \n" << matrix_cplx_sparse.irow << std::endl;
//    std::cout << "pcol \n" << matrix_cplx_sparse.pcol << std::endl;
//    std::cout << "vals \n" << matrix_cplx_sparse.vals << std::endl;
//    auto innernnz = std::vector<int>   (matrix_recast.innerNonZeroPtr(), matrix_recast.innerNonZeroPtr() + matrix_recast.innerSize() );
//    auto innernnz = std::vector<int>   (matrix_recast.innerNonZeroPtr(), matrix_recast.innerIndexPtr() + matrix_recast.innerSize() );
//    std::cout << "innernnz: \n" << innernnz << std::endl;

//}



template class matrix_recast<double>;
template class matrix_recast<std::complex<double>>;