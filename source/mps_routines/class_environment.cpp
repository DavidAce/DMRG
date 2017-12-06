//
// Created by david on 2017-12-02.
//


#include "class_environment.h"
template<Side side>
void class_environment<side>::enlarge(const class_mps &MPS, const Tensor<Scalar,4> &M){
    /*!< Contracts a site into the block. */
    Tensor<Scalar,3> block_enlarged;
    if constexpr(side == Side::L){

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
                    block.contract(asDiagonal(MPS.L_tail), idx<1>({0},{0}))
                            .contract(MPS.GA.conjugate(),     idx<1>({2},{1}))
                            .contract(M,                      idx<2>({1,2},{0,2}))
                            .contract(asDiagonal(MPS.L_tail), idx<1>({0},{0}))
                            .contract(MPS.GA,                 idx<2>({2,3},{0,1}))
                            .shuffle(array3{0,2,1});
            block = block_enlarged;
    }else if constexpr(side== Side::R){
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
                    block.contract(asDiagonal(MPS.LB), idx<1>({0},{1}))
                            .contract(MPS.GB.conjugate(), idx<1>({2},{2}))
                            .contract(M,                  idx<2>({1,2},{1,2}))
                            .contract(asDiagonal(MPS.LB), idx<1>({0},{1}))
                            .contract(MPS.GB,             idx<2>({2,3},{0,2}))
                            .shuffle(array3{0,2,1});
            block = block_enlarged;
    }
}

template class class_environment<Side::L>;
template class class_environment<Side::R>;
