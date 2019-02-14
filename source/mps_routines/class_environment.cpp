//
// Created by david on 2017-12-02.
//


#include "class_environment.h"
#include <mps_routines/class_mps_2site.h>

using namespace std;
using namespace Textra;
using Scalar = class_environment::Scalar;
void class_environment::enlarge(class_mps_2site & MPS, const Eigen::Tensor<Scalar,4> &M){
    /*!< Contracts a site into the block. */
    Eigen::Tensor<Scalar,3> block_enlarged;
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
//            position = MPS.MPS_A->get_position();
        block_enlarged =
                block.contract(MPS.A(),                  idx({0},{1}))
                        .contract(M,                      idx({1,2},{0,2}))
                        .contract(MPS.A().conjugate(),   idx({0,3},{1,0}))
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
//        position = MPS.MPS_B->get_position();
        block_enlarged =
                block.contract(MPS.B(),                idx({0},{2}))
                        .contract(M,                    idx({1,2},{1,2}))
                        .contract(MPS.B().conjugate(), idx({0,3},{2,0}))
                        .shuffle(array3{0,2,1});
        block = block_enlarged;
    }
}

void class_environment::set_edge_dims(class_mps_2site & MPS, const Eigen::Tensor<Scalar, 4> &M) {
    if (side == "L") {
        long chiA = MPS.chiA();
        block.resize(array3{chiA,chiA, M.dimension(0)});
        block.setZero();
        for (long i = 0; i < chiA; i++){
            block(i,i,M.dimension(0)-1) = 1;
        }
        position = 0;
    }
    if(side == "R"){
        long chiB = MPS.chiB();
        block.resize(array3{chiB,chiB, M.dimension(1)});
        block.setZero();
        for (long i = 0; i < chiB; i++){
            block(i,i,0) = 1;
        }
        position = 1;
    }
//    enlarge(MPS,M);
    size = 0;
}


void class_environment_var::enlarge(class_mps_2site & MPS, const Eigen::Tensor<Scalar,4> &M){
    /*!< Contracts a site into the block. */
    Eigen::Tensor<Scalar,4> block_enlarged;
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
//        position = MPS.MPS_A->get_position();
        block_enlarged =
                block.contract(MPS.A(),                    idx({0},{1}))
                        .contract(M,                        idx({1,3},{0,2}))
                        .contract(M,                        idx({1,4},{0,2}))
                        .contract(MPS.A().conjugate(),      idx({0,4},{1,0}))
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
//        position = MPS.MPS_B->get_position();
        block_enlarged =
                block.contract(MPS.B(),                idx({0},{2}))
                        .contract(M,                    idx({1,3},{1,2}))
                        .contract(M,                    idx({1,4},{1,2}))
                        .contract(MPS.B().conjugate(), idx({0,4},{2,0}))
                        .shuffle(array4{0, 3, 1, 2});
        block = block_enlarged;
    }
}



void class_environment_var::set_edge_dims(class_mps_2site & MPS, const Eigen::Tensor<Scalar, 4> &M) {
    if (side == "L") {
        long chiA = MPS.chiA();
        block.resize(array4{chiA,chiA, M.dimension(0), M.dimension(0)});
        block.setZero();
        for (long i = 0; i < chiA; i++){
            block(i,i, M.dimension(0)-1, M.dimension(0)-1) = 1;
        }
        position = 0;

    }
    if(side == "R"){
        long chiB = MPS.chiB();
        block.resize(array4{chiB,chiB, M.dimension(1),M.dimension(1)});
        block.setZero();
        for (long i = 0; i < chiB; i++){
            block(i,i,0,0) = 1;
        }
        position = 1;
    }
//    enlarge(MPS,M);
    size = 0;
}