//
// Created by david on 2017-12-02.
//


#include "class_environment.h"
#include <mps_state/class_mps_2site.h>

using namespace std;
using namespace Textra;
using Scalar = class_environment::Scalar;

bool class_environment::isReal() const {
    return Textra::isReal(block, "env " + side);
}
bool class_environment_var::isReal() const {
    return Textra::isReal(block, "env2" + side);
}


void class_environment::enlarge(const class_mps_2site & MPS, const Eigen::Tensor<Scalar,4> &MPO){
    if (side == "L"){
        enlarge(MPS.A(),MPO);
    }else if (side == "R"){
        enlarge(MPS.B(),MPO);
    }
}

void class_environment::enlarge(const Eigen::Tensor<Scalar,3> MPS, const Eigen::Tensor<Scalar,4> &MPO){
    /*!< Contracts a site into the block. */
    if(size == 0){set_edge_dims(MPS,MPO);}

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
                block.contract(MPS,                  idx({0},{1}))
                        .contract(MPO,                      idx({1,2},{0,2}))
                        .contract(MPS.conjugate(),   idx({0,3},{1,0}))
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
                block.contract(MPS,                idx({0},{2}))
                        .contract(MPO,                    idx({1,2},{1,2}))
                        .contract(MPS.conjugate(), idx({0,3},{2,0}))
                        .shuffle(array3{0,2,1});
        block = block_enlarged;
    }
}

void class_environment::set_edge_dims(const class_mps_2site & MPS, const Eigen::Tensor<Scalar, 4> &MPO) {
    if (side == "L") {
        set_edge_dims(MPS.chiA(),MPO.dimension(0));
    }
    if(side == "R"){
        set_edge_dims(MPS.chiB(),MPO.dimension(1));
    }
}

void class_environment::set_edge_dims(const Eigen::Tensor<Scalar,3> MPS, const Eigen::Tensor<Scalar, 4> &MPO) {
    if (side == "L") {
        set_edge_dims(MPS.dimension(1),MPO.dimension(0));
    }
    if(side == "R"){
        set_edge_dims(MPS.dimension(2),MPO.dimension(1));
    }
}

void class_environment::set_edge_dims(const int mpsDim, const int mpoDim) {
    if (side == "L") {
        block.resize(array3{mpsDim,mpsDim, mpoDim});
        block.setZero();
        for (long i = 0; i < mpsDim; i++){
            block(i,i,mpoDim-1) = 1;
        }
        position = 0;
    }
    if(side == "R"){
        block.resize(array3{mpsDim,mpsDim, mpoDim});
        block.setZero();
        for (long i = 0; i < mpsDim; i++){
            block(i,i,0) = 1;
        }
        position = 1;
    }
    size = 0;
}




void class_environment_var::enlarge(const class_mps_2site & MPS, const Eigen::Tensor<Scalar,4> &MPO){
    if (side == "L"){
        enlarge(MPS.A(),MPO);
    }else if (side == "R"){
        enlarge(MPS.B(),MPO);
    }
}

void class_environment_var::enlarge(const Eigen::Tensor<Scalar,3> MPS, const Eigen::Tensor<Scalar,4> &MPO){
    /*!< Contracts a site into the block. */
    if(size == 0){set_edge_dims(MPS,MPO);}
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
                block.contract(MPS,                    idx({0},{1}))
                        .contract(MPO,                        idx({1,3},{0,2}))
                        .contract(MPO,                        idx({1,4},{0,2}))
                        .contract(MPS.conjugate(),      idx({0,4},{1,0}))
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
                block.contract(MPS,                idx({0},{2}))
                        .contract(MPO,                    idx({1,3},{1,2}))
                        .contract(MPO,                    idx({1,4},{1,2}))
                        .contract(MPS.conjugate(), idx({0,4},{2,0}))
                        .shuffle(array4{0, 3, 1, 2});
        block = block_enlarged;
    }
}

void class_environment_var::set_edge_dims(const class_mps_2site & MPS, const Eigen::Tensor<Scalar, 4> &MPO) {
    if (side == "L") {
        set_edge_dims(MPS.chiA(),MPO.dimension(0));
    }
    if(side == "R"){
        set_edge_dims(MPS.chiB(),MPO.dimension(1));
    }
}

void class_environment_var::set_edge_dims(const Eigen::Tensor<Scalar,3> MPS, const Eigen::Tensor<Scalar, 4> &MPO) {
    if (side == "L") {
        set_edge_dims(MPS.dimension(1),MPO.dimension(0));
    }
    if(side == "R"){
        set_edge_dims(MPS.dimension(2),MPO.dimension(1));
    }
}

void class_environment_var::set_edge_dims(const int mpsDim, const int mpoDim) {
    block.resize(array4{mpsDim,mpsDim, mpoDim,mpoDim});
    block.setZero();
    if (side == "L") {
        for (long i = 0; i < mpsDim; i++){
            block(i,i, mpoDim-1, mpoDim-1) = 1;
        }
        position = 0;
    }
    if(side == "R"){

        for (long i = 0; i < mpsDim; i++){
            block(i,i,0,0) = 1;
        }
        position = 1;
    }
    size = 0;
}

