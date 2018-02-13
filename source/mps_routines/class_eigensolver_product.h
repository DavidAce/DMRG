//
// Created by david on 2017-11-15.
//

#ifndef DMRG_CLASS_GROUNDSTATE_SEARCH_H
#define DMRG_CLASS_GROUNDSTATE_SEARCH_H



#include "general/n_tensor_extra.h"

/*!
 * \class class_eigensolver_product
 * \brief Class that the tensor-vector product needed in the eigensolver Spectra.
 *  Defines the matrix-vector product in the left side of \f$Av = \lambda \f$
 */
class class_eigensolver_product {
private:
    using Scalar = double;
    const Textra::Tensor<Scalar,3> &Lblock;
    const Textra::Tensor<Scalar,3> &Rblock;
    const Textra::Tensor<Scalar,6> &MM;
    const Textra::array4         &shape4;
          Textra::array1          shape1;

    Textra::idxlistpair<long, 3> sortedIndices1;
    Textra::idxlistpair<long, 2> sortedIndices2;
public:
    explicit class_eigensolver_product(const Textra::Tensor<Scalar,3> &Lblock_,
                                       const Textra::Tensor<Scalar,3> &Rblock_,
                                       const Textra::Tensor<Scalar,6> &MM_,
                                       const Textra::array4 &shape4_);


    int rows()const;
    int cols()const;
    void perform_op(const double *x_in, double *y_out) const;
};


#endif //DMRG_CLASS_GROUNDSTATE_SEARCH_H
