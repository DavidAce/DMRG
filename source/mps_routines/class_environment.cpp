//
// Created by david on 2017-12-02.
//


#include "class_environment.h"
#include <mps_routines/class_mps.h>

using namespace std;
using namespace Textra;
using Scalar = class_environment::Scalar;
void class_environment::enlarge(const std::shared_ptr<class_mps> &MPS, const Tensor<Scalar,4> &M){
    /*!< Contracts a site into the block. */
    Tensor<Scalar,3> block_enlarged;
    if (side == "L"){

            /*! # Left environment contraction
             * [      ]--0 0--[LB]--1 1--[  GA    ]--2
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
             * [      ]--1 0--[LB]--1  1--[GA conj ]--2
             */
            size++;
            block_enlarged =
                    block.contract(asDiagonal(MPS->LB_left),   idx({0},{0}))
                            .contract(MPS->GA,                idx({2},{1}))
                            .contract(M,                      idx({1,2},{0,2}))
                            .contract(asDiagonal(MPS->LB_left),idx({0},{0}))
                            .contract(MPS->GA.conjugate(),    idx({3,2},{1,0}))
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
                    block.contract(asDiagonal(MPS->LB),    idx({0},{1}))
                            .contract(MPS->GB,             idx({2},{2}))
                            .contract(M,                   idx({1,2},{1,2}))
                            .contract(asDiagonal(MPS->LB), idx({0},{1}))
                            .contract(MPS->GB.conjugate(), idx({3,2},{2,0}))
                            .shuffle(array3{0,2,1});
            block = block_enlarged;
    }
}

void class_environment::set_edge_dims(const std::shared_ptr<class_mps> &MPS, const Tensor<Scalar, 4> &M) {
    if (side == "L") {
        long chiL = MPS->LB_left.dimension(0);
        block.resize(array3{chiL,chiL, M.dimension(0)});
        block.setZero();
        for (long i = 0; i < chiL; i++){
            block(i,i,2) = 1;
        }
    }
    if(side == "R"){
        long chiR = MPS->LB.dimension(0);
        block.resize(array3{chiR,chiR, M.dimension(1)});
        block.setZero();
        for (long i = 0; i < chiR; i++){
            block(i,i,0) = 1;
        }
    }
//    enlarge(MPS,M);
    size = 0;
}


void class_environment_var::enlarge(const std::shared_ptr<class_mps> &MPS, const Tensor<Scalar,4> &M){
    /*!< Contracts a site into the block. */
    Tensor<Scalar,4> block_enlarged;
    if (side == "L"){

        /*! # Left environment contraction
         * [      ]--0 0--[LB]--1 1--[  GA   ]--2
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
         * [      ]--1 0--[LB]--1  1--[GA conj ]--2
         */
        size++;
//        std::cout << "LB_left: " << MPS->LB_left.dimensions() << std::endl;
//        std::cout << "GA     : " << MPS->GA.dimensions() << std::endl;
//        std::cout << "M      : " << M.dimensions() << std::endl;

        block_enlarged =
                block.contract(asDiagonal(MPS->LB_left),    idx({0},{0}))
                        .contract(MPS->GA,                  idx({3},{1}))
                        .contract(M,                        idx({1,3},{0,2}))
                        .contract(M,                        idx({1,4},{0,2}))
                        .contract(asDiagonal(MPS->LB_left), idx({0},{0}))
                        .contract(MPS->GA.conjugate(),      idx({4,3},{1,0}))
                        .shuffle(array4{0,3,1,2});
        block = block_enlarged;

    }
    if (side == "R"){
        /*! # Right environment contraction
         *  1--[   GB    ]--2 0--[LB]--1  0--[      ]
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
         *  1--[ GB conj ]--2 0--[LB]--1  1--[      ]
        */

        size++;
        block_enlarged =
                block.contract(asDiagonal(MPS->LB),    idx({0},{1}))
                        .contract(MPS->GB,             idx({3},{2}))
                        .contract(M,                   idx({1,3},{1,2}))
                        .contract(M,                   idx({1,4},{1,2}))
                        .contract(asDiagonal(MPS->LB), idx({0},{1}))
                        .contract(MPS->GB.conjugate(), idx({4,3},{2,0}))
                        .shuffle(array4{0, 3, 1, 2});
        block = block_enlarged;
    }
}



void class_environment_var::set_edge_dims(const std::shared_ptr<class_mps> &MPS, const Tensor<Scalar, 4> &M) {
    if (side == "L") {
        long chiL = MPS->LB_left.dimension(0);
        block.resize(array4{chiL,chiL, M.dimension(0), M.dimension(0)});
        block.setZero();
        for (long i = 0; i < chiL; i++){
            block(i,i,2,2) = 1;
        }

    }
    if(side == "R"){
        long chiR = MPS->LB.dimension(0);
        block.resize(array4{chiR,chiR, M.dimension(1),M.dimension(1)});
        block.setZero();
        for (long i = 0; i < chiR; i++){
            block(i,i,0,0) = 1;
        }
    }
//    enlarge(MPS,M);
    size = 0;
}