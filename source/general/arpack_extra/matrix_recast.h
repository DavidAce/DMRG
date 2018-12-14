//
// Created by david on 2018-10-29.
//

#ifndef MATRIX_RECAST_H
#define MATRIX_RECAST_H

#include <complex>
#include <vector>
#include "general/nmspc_eigutils.h"
#include "matrix_product_dense.h"
#include "matrix_product_sparse.h"


template<typename Scalar>
class matrix_recast {
private:


    bool isSparse;
    bool isReal;
    bool isHermitian;
    double sparcity;


    const Scalar *matrix_ptr;
    std::vector<Scalar> matrix_pruned;
    bool pruned = false;
    int     L;

    void check_if_real();
    void check_if_sparse();
    void check_if_hermitian();
    void recheck_all();

//    eigutils::DenseType<double>                matrix_real_dense;
//    eigutils::DenseType<std::complex<double>>  matrix_cplx_dense;
//    eigutils::SparseType<double>               matrix_real_sparse;
//    eigutils::SparseType<std::complex<double>> matrix_cplx_sparse;
//    DenseMatrixProduct<double>                matrix_real_dense;
//    DenseMatrixProduct<std::complex<double>>  matrix_cplx_dense;
//    SparseMatrixProduct<double>               matrix_real_sparse;
//    SparseMatrixProduct<std::complex<double>> matrix_cplx_sparse;


public:
    matrix_recast(const Scalar *matrix_ptr_, int L_);
    void prune(double threshold = 1e-14);
    DenseMatrixProduct<double>                get_as_real_dense();
    DenseMatrixProduct<std::complex<double>>  get_as_cplx_dense();
    SparseMatrixProduct<double>               get_as_real_sparse();
    SparseMatrixProduct<std::complex<double>> get_as_cplx_sparse();

//    void convert_to_real_dense();
//    void convert_to_cplx_dense();
//    void convert_to_real_sparse();
//    void convert_to_cplx_sparse();

    bool is_sparse      ()   {return isSparse;}
    bool is_real        ()   {return isReal;}
    bool is_symmetric()   {return isHermitian;}
};


#endif //MATRIX_RECAST_H
