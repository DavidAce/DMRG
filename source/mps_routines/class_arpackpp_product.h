//
// Created by david on 2018-02-14.
//

#ifndef DMRG_CLASS_EIGENSOLVER_PRODUCT_ARPACK_H
#define DMRG_CLASS_EIGENSOLVER_PRODUCT_ARPACK_H

#ifdef MKL_AVAILABLE
#define  EIGEN_USE_MKL_ALL
#endif


#include "general/nmspc_tensor_extra.h"

template<class T>
class class_arpackpp_product {

private:
    using Scalar = double;
    const Textra::Tensor<Scalar,1> &theta;
    const Textra::Tensor<Scalar,3> &Lblock;
    const Textra::Tensor<Scalar,3> &Rblock;
    const Textra::Tensor<Scalar,6> &MM;
    const Textra::array4           &shape4;
    Textra::array1                 shape1;

    Textra::idxlistpair<long, 3> sortedIndices1;
    Textra::idxlistpair<long, 2> sortedIndices2;
public:
    int rows()const {return (int)shape1[0];};
    int cols()const {return (int)shape1[0];};
    void MultMv(T* v, T* w);

//    explicit class_arpackpp_product(int nxval): MatrixWithProduct<T>(nxval*nxval) { nx = nxval; }
    class_arpackpp_product(
            const Textra::Tensor<Scalar,1> &theta_,
            const Textra::Tensor<Scalar,3> &Lblock_,
            const Textra::Tensor<Scalar,3> &Rblock_,
            const Textra::Tensor<Scalar,6> &MM_,
            const Textra::array4           &shape4_):
//            MatrixWithProduct<T>((int)theta_.size()),
            theta(theta_),
            Lblock(Lblock_),
            Rblock(Rblock_),
            MM(MM_),
            shape4(shape4_)
    {
        shape1[0] = shape4[0] * shape4[1] * shape4[2] * shape4[3];
        sortedIndices1 = Textra::sortIdx<3, 6>(MM.dimensions(), {1, 2, 3}, {0, 4, 5});
        sortedIndices2 = Textra::sortIdx<2, 3>(Rblock.dimensions(), {1, 2}, {1, 2});
    }

};


template<class ART>
void class_arpackpp_product<ART>::MultMv(ART* x_in, ART* y_out)
/*
  Computes w <- A*v.
*/

{
    Eigen::TensorMap<Textra::Tensor<const double, 4>> x_input (x_in, shape4);
    Eigen::TensorMap<Textra::Tensor<double, 1>>       y_output(y_out, shape1);
    //Best yet! The sparcity of the effective hamiltonian (Lblock MM Rblock) is about 58% nonzeros.
    y_output = Lblock
            .contract(x_input, Textra::idx<1>({1},{1}))
            .contract(MM ,     sortedIndices1)//  idx<3>({1,2,3},{0,4,5}))
            .contract(Rblock,  sortedIndices2)
            .shuffle(Textra::array4{1,0,2,3})
            .reshape(shape1);

} //  MultMv.





#endif //DMRG_CLASS_EIGENSOLVER_PRODUCT_ARPACK_H
