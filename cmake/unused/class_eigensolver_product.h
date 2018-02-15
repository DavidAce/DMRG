//
// Created by david on 2017-11-15.
//

#ifndef DMRG_EIGENSOLVER_PRODUCT_H
#define DMRG_EIGENSOLVER_PRODUCT_H


#ifdef MKL_AVAILABLE
#define  EIGEN_USE_MKL_ALL
#endif

#include "general/nmspc_tensor_extra.h"
/*!
 * \class class_eigensolver_product
 * \brief Class that the tensor-vector product needed in the eigensolver Spectra.
 *  Defines the matrix-vector product in the left side of \f$Av = \lambda \f$
 */
class class_eigensolver_product {
private:
    using Scalar = double;
    const Textra::Tensor<Scalar,1> &theta;
    const Textra::Tensor<Scalar,3> &Lblock;
    const Textra::Tensor<Scalar,3> &Rblock;
    const Textra::Tensor<Scalar,6> &MM;
    const Textra::array4           &shape4;
          Textra::array1           shape1;
    Textra::Tensor<Scalar,2> test2;
public:
    Eigen::MatrixXd test;
//    Eigen::SparseMatrix<double> test2;

    Textra::idxlistpair<long, 3> sortedIndices1;
    Textra::idxlistpair<long, 2> sortedIndices2;

    typedef Eigen::Array<Scalar, Eigen::Dynamic,1> Vector;
    typedef Eigen::Map<const Vector> MapConstVec;
    typedef Eigen::Map<Vector> MapVec;
//    typedef Eigen::SparseMatrix<Scalar, Flags, StorageIndex> SparseMatrix;
public:
    explicit class_eigensolver_product(const Textra::Tensor<Scalar,1> &theta_,
                                       const Textra::Tensor<Scalar,3> &Lblock_,
                                       const Textra::Tensor<Scalar,3> &Rblock_,
                                       const Textra::Tensor<Scalar,6> &MM_,
                                       const Textra::array4 &shape4_);


    int rows()const;
    int cols()const;
    void perform_op(const double *x_in, double *y_out);
};


#endif //DMRG_CLASS_GROUNDSTATE_SEARCH_H
