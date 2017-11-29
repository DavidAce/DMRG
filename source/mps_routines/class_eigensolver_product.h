//
// Created by david on 2017-11-15.
//

#ifndef DMRG_CLASS_GROUNDSTATE_SEARCH_H
#define DMRG_CLASS_GROUNDSTATE_SEARCH_H



#include "general/n_tensor_extra.h"
using namespace Textra;
class class_eigensolver_product {
private:
    using Scalar = double;
    const Tensor<Scalar,3> &Lblock;
    const Tensor<Scalar,3> &Rblock;
    const Tensor<Scalar,6> &MM;
    const array4         &shape4;
          array1          shape1;

    idxlistpair<long, 3> sortedIndices1;
    idxlistpair<long, 2> sortedIndices2;
public:
    explicit class_eigensolver_product(const Tensor<Scalar,3> &Lblock_,
                                       const Tensor<Scalar,3> &Rblock_,
                                       const Tensor<Scalar,6> &MM_,
                                       const array4 &shape4_);


    /*! Function for eigenvalue solver Spectra
     *  Defines the matrix-vector product in the left side of \f$Av = \lambda \f$
     */
    int rows()const;
    int cols()const;
    void perform_op(const double *x_in, double *y_out) const;
};


#endif //DMRG_CLASS_GROUNDSTATE_SEARCH_H
