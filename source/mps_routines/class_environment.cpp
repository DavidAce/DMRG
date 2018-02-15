//
// Created by david on 2017-12-02.
//


#ifdef MKL_AVAILABLE
#define  EIGEN_USE_MKL_ALL
#endif

#include "class_environment.h"
#include <mps_routines/class_mps.h>

using namespace std;
using namespace Textra;

template<typename Scalar>
void class_environment<Scalar>::enlarge(const std::shared_ptr<class_mps<Scalar>> &MPS, const Tensor<Scalar,4> &M){
    /*!< Contracts a site into the block. */
    Tensor<Scalar,3> block_enlarged;
    if (side == "L"){

            /*! # Left environment contraction
             * [      ]--0 0--[LB]--1 1--[ GA conj ]--2
             * [      ]                      |
             * [      ]                      0
             * [      ]
             * [      ]                      2
             * [      ]                      |
             * [ left ]--2            0--[   M    ]--1
             * [      ]                      |
             * [      ]                      3
             * [      ]
             * [      ]                      0
             * [      ]                      |
             * [      ]--1 0--[LB]--1  1--[  GA   ]--2
             */
            size++;
            block_enlarged =
                    block.contract(asDiagonal(MPS->L_tail), idx<1>({0},{0}))
                            .contract(MPS->GA.conjugate(),     idx<1>({2},{1}))
                            .contract(M,                      idx<2>({1,2},{0,2}))
                            .contract(asDiagonal(MPS->L_tail), idx<1>({0},{0}))
                            .contract(MPS->GA,                 idx<2>({3,2},{1,0}))
                            .shuffle(array3{0,2,1});
            block = block_enlarged;
    }else if (side== "R"){
            /*! # Right environment contraction
             *  1--[ GB conj ]--2 0--[LB]--1  0--[      ]
             *          |                        [      ]
             *          0                        [      ]
             *                                   [      ]
             *          2                        [      ]
             *          |                        [      ]
             *   0--[   M    ]--1             2--[ right]
             *          |                        [      ]
             *          3                        [      ]
             *                                   [      ]
             *          0                        [      ]
             *          |                        [      ]
             *    1--[  GB   ]--2 0--[LB]--1  1--[      ]
            */

            size++;
            block_enlarged =
                    block.contract(asDiagonal(MPS->LB),    idx<1>({0},{1}))
                            .contract(MPS->GB.conjugate(), idx<1>({2},{2}))
                            .contract(M,                  idx<2>({1,2},{1,2}))
                            .contract(asDiagonal(MPS->LB), idx<1>({0},{1}))
                            .contract(MPS->GB,             idx<2>({3,2},{2,0}))
                            .shuffle(array3{0,2,1});
            block = block_enlarged;
    }
}


template<typename Scalar>
void class_environment_var<Scalar>::enlarge(const std::shared_ptr<class_mps<Scalar>> &MPS, const Tensor<Scalar,4> &M){
    /*!< Contracts a site into the block. */
    Tensor<Scalar,4> block_enlarged;
    if (side == "L"){

        /*! # Left environment contraction
         * [      ]--0 0--[LB]--1 1--[ GA conj ]--2
         * [      ]                      |
         * [      ]                      0
         * [      ]
         * [      ]                      2
         * [      ]                      |
         * [      ]--2            0--[   M    ]--1
         * [      ]                      |
         * [      ]                      3
         * [ left ]
         * [      ]                      2
         * [      ]                      |
         * [      ]--3            0--[   M    ]--1
         * [      ]                      |
         * [      ]                      3
         * [      ]
         * [      ]                      0
         * [      ]                      |
         * [      ]--1 0--[LB]--1  1--[  GA   ]--2
         */
        size++;
        block_enlarged =
                block.contract(asDiagonal(MPS->L_tail),    idx<1>({0},{0}))
                        .contract(MPS->GA.conjugate(),     idx<1>({3},{1}))
                        .contract(M,                      idx<2>({1,3},{0,2}))
                        .contract(M,                      idx<2>({1,4},{0,2}))
                        .contract(asDiagonal(MPS->L_tail), idx<1>({0},{0}))
                        .contract(MPS->GA,                 idx<2>({4,3},{1,0}))
                        .shuffle(array4{0,3,1,2});
        block = block_enlarged;

    }
    if (side == "R"){
        /*! # Right environment contraction
         *  1--[ GB conj ]--2 0--[LB]--1  0--[      ]
         *          |                        [      ]
         *          0                        [      ]
         *                                   [      ]
         *          2                        [      ]
         *          |                        [      ]
         *   0--[   M    ]--1             2--[      ]
         *          |                        [      ]
         *          3                        [      ]
         *                                   [ right]
         *          2                        [      ]
         *          |                        [      ]
         *   0--[   M    ]--1             2--[      ]
         *          |                        [      ]
         *          3                        [      ]
         *                                   [      ]
         *          0                        [      ]
         *          |                        [      ]
         *    1--[  GB   ]--2 0--[LB]--1  1--[      ]
        */

        size++;
        block_enlarged =
                block.contract(asDiagonal(MPS->LB), idx<1>({0},{1}))
                        .contract(MPS->GB.conjugate(), idx<1>({3},{2}))
                        .contract(M,                  idx<2>({1,3},{1,2}))
                        .contract(M,                  idx<2>({1,4},{1,2}))
                        .contract(asDiagonal(MPS->LB), idx<1>({0},{1}))
                        .contract(MPS->GB,             idx<2>({4,3},{2,0}))
                        .shuffle(array4{0, 3, 1, 2});
        block = block_enlarged;
    }
}



template class class_environment<double>;
template class class_environment_var<double>;

//template class class_environment<std::complex<double>, Side::L>;
//template class class_environment<std::complex<double>, Side::R>;
//template class class_environment_var<std::complex<double>, Side::L>;
//template class class_environment_var<std::complex<double>, Side::R>;