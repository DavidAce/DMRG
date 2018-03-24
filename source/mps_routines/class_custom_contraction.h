//
// Created by david on 2018-02-14.
//

#ifndef DMRG_CLASS_CUSTOM_CONTRACTION_H
#define DMRG_CLASS_CUSTOM_CONTRACTION_H

#ifdef MKL_AVAILABLE
#define  EIGEN_USE_MKL_ALL
#endif

#include "general/nmspc_tensor_extra.h"
#include <iomanip>

/*!
 \class class_custom_contraction
 \brief A class that defines a custom product for the iterative eigenvalue solver in arpackpp.

 This class contracts the tensors involved in the eigenvalue equation of the form \f$ Ax = kx \f$ .
 In this case, the operator \f$ A \f$ is the effective Hamiltonian at the current sites,
 \f$ x \f$ is the sought eigenvector <code> theta </code> and the eigenvalue \f$ k \f$ is the energy \f$ E \f$ for the whole chain.
 Thus the eigenvalue equation can be written as \f$ H\theta = E\theta \f$

 The definition of <code> theta </code> is shown below. The index order of theta is arbitrary as long as \f$ H \f$ follows the same
 convention, but this particular choice facilitates the SVD decomposition that will follow next.

    @verbatim
             0   2                           0                   0
             |   |      =                    |                   |
        1--[ theta ]--3   0--[LB]--1  1--[  GA   ]--2     1--[  GB   ]--2 0--[LB]--1
    @endverbatim

 The definition of \f$ H \f$ is shown below following the same index order as theta.

    @verbatim
                                  [      ]--0                                            0--[      ]
                                  [      ]                                                  [      ]
                                  [      ]                                                  [      ]
        |-1         3-|           [      ]                                                  [      ]
        |    0   2    |           [      ]              2                   2               [      ]
        |    |   |    |           [      ]              |                   |               [      ]
        |--[   H   ]--|    =      [ left ]--2    0--[   M    ]--1    0--[   M    ]--1    2--[ right]
        |    |   |    |           [      ]              |                   |               [      ]
        |    4   6    |           [      ]              3                   3               [      ]
        |-5         7-|           [      ]                                                  [      ]
                                  [      ]                                                  [      ]
                                  [      ]                                                  [      ]
                                  [      ]--1                                            1--[      ]
    @endverbatim

 */

template<class T>
class class_custom_contraction {

private:
    using Scalar = std::complex<double>;
    const Textra::Tensor<Scalar,3> &Lblock;
    const Textra::Tensor<Scalar,3> &Rblock;
    const Textra::Tensor<Scalar,4> &M;
    const Textra::array4           &shape4;
    Textra::array1                 shape1;

//    Textra::idxlistpair<long, 3> sortedIndices1;
//    Textra::idxlistpair<long, 2> sortedIndices2;
public:
    int rows()const {return (int)shape1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
    int cols()const {return (int)shape1[0];};               /*!< The "matrix" \f$ H \f$ a has rows = columns = \f$d^2 \times \chi_L \times \chi_R \f$  */
    void MultMv(T* theta_in_, T* theta_out_);               /*!< The function that contracts.  */

    class_custom_contraction(
            const Textra::Tensor<Scalar,3> &Lblock_,        /*!< The left block tensor.  */
            const Textra::Tensor<Scalar,3> &Rblock_,        /*!< The right block tensor.  */
            const Textra::Tensor<Scalar,4> &MM_,            /*!< The two Hamiltonian MPO's already contracted into one  */
            const Textra::array4           &shape4_         /*!< An array containing the shapes of theta  */
    ):                                                      /*!< Initializes the custom contraction. */
            Lblock(Lblock_),
            Rblock(Rblock_),
            M(MM_),
            shape4(shape4_),
            shape1({shape4[0] * shape4[1] * shape4[2] * shape4[3]})
    {
//        sortedIndices1 = Textra::sortIdx(M.dimensions(), {1, 2, 3}, {0, 2, 3});
//        sortedIndices2 = Textra::sortIdx(Rblock.dimensions(), {1, 2}, {0, 2});

    }

};


template<class T>
void class_custom_contraction<T>::MultMv(T* theta_in_, T* theta_out_) {
    Eigen::TensorMap<Textra::Tensor<const T, 4>> theta_in (theta_in_, shape4);
    Eigen::TensorMap<Textra::Tensor<T, 1>>       theta_out(theta_out_, shape1);
    //Best yet! The sparcity of the effective hamiltonian (Lblock M Rblock) is about 58% nonzeros.
//    theta_out = Lblock
//            .contract(theta_in, Textra::idx({0},{1}))
//            .contract(M ,     sortedIndices1)//  idx({1,2,3},{0,4,5}))
//            .contract(Rblock,  sortedIndices2)
//            .shuffle(Textra::array4{1,0,2,3})
//            .reshape(shape1);
    theta_out = Lblock
            .contract(theta_in, Textra::idx({0},{1}))
            .contract(M ,     Textra::idx({1,2},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(M ,     Textra::idx({3,1},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(Rblock,  Textra::idx({1,3},{0,2}))
            .shuffle(Textra::array4{1,0,2,3})
            .reshape(shape1);
}





#endif
