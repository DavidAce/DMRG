//
// Created by david on 2017-11-15.
//

#include "class_eigensolver_product.h"
class_eigensolver_product::class_eigensolver_product(
        const Textra::Tensor<Scalar,3> &Lblock_,
        const Textra::Tensor<Scalar,3> &Rblock_,
        const Textra::Tensor<Scalar,6> &MM_,
        const Textra::array4           &shape4_):
            Lblock(Lblock_),
            Rblock(Rblock_),
            MM(MM_),
            shape4(shape4_)
{
    shape1[0] = shape4[0] * shape4[1] * shape4[2] * shape4[3];
    sortedIndices1 = Textra::sortIdx<3, 6>(MM.dimensions(), {1, 2, 3}, {0, 4, 5});
    sortedIndices2 = Textra::sortIdx<2, 3>(Rblock.dimensions(), {1, 2}, {1, 2});
}

//Functions for eigenvalue solver Spectra
int  class_eigensolver_product::rows()const {return (int)shape1[0];}
int  class_eigensolver_product::cols()const {return (int)shape1[0];}
void class_eigensolver_product::perform_op(const double *x_in, double *y_out) const {
    Eigen::TensorMap<Textra::Tensor<const double, 4>> x_input (x_in, shape4);
    Eigen::TensorMap<Textra::Tensor<double, 1>> y_output (y_out, shape1);
    //Best yet!
    y_output = Lblock
            .contract(x_input, Textra::idx<1>({1},{1}))
            .contract(MM ,     sortedIndices1)//  idx<3>({1,2,3},{0,4,5}))
            .contract(Rblock,  sortedIndices2)
            .shuffle(Textra::array4{1,0,2,3})
            .reshape(shape1);
}