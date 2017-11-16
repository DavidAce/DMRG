//
// Created by david on 7/21/17.
//

#ifndef DMRG_CLASS_ENVIRONMENT_H
#define DMRG_CLASS_ENVIRONMENT_H

#include "general/n_tensor_extra.h"
#include "mps_routines/class_mps.h"
#include "sim_parameters/n_model.h"

using namespace Textra;
using namespace std;

/*! \brief Environment block och type Left or Right.
 *
 * # Left environment
 * This class contains the Left environment block as a rank-3 tensor. New sites are contracted into the left environment as
 * \f[
 * L \leftarrow L \Lambda^B_{n-1} \Gamma^A_n W \Lambda^B_{n-1} (\Gamma^A_n)^*
 * \f]
 * where \f$W\f$ is the rank-4 tensor Hamiltonian MPO.

 *
 * # Right environment
 * This class contains the Right environment block as a rank-3 tensor. New sites are contracted into the Right environment as
 * \f[
 * R \leftarrow R \Gamma^B_{n+1} \Lambda^B_{n+1} W (\Gamma^B_{n+1})^* \Lambda^B_{n+1}
 * \f]
 * where \f$W\f$ is the rank-4 tensor Hamiltonian MPO.
 */

enum class Side {L,R};
template<Side side>
class class_environment{
public:
    using Scalar = Model::Scalar;
private:

public:
    int size;                                               /*!< Number of particles that have been contracted into this left environment. */
    Textra::Tensor<3,Scalar> block;                         /*!< The environment block. */
    class_environment(){
        size = 0;
        block.resize(1, 1, 3);
        block.setZero();
        switch (side){
            case Side::L:
                block(0, 0, 0) = 1;
                break;
            case Side::R:
                block(0, 0, 2) = 1;
                break;
        }
    };

    void enlarge(const class_mps &MPS, const Textra::Tensor<4,Scalar> &W){
        /*!< Contracts a site into the block. */
        Textra::Tensor<3,Scalar> block_enlarged;
        switch (side){
            case (Side::L):
                /*! # Left environment contraction
                 * [      ]--0 0--[LB]--1 1--[ GA conj ]--2
                 * [      ]                      |
                 * [      ]                      0
                 * [      ]
                 * [      ]                      2
                 * [      ]                      |
                 * [ left ]--2            0--[   W    ]--1
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
                                .contract(W,                      idx<2>({1,2},{0,2}))
                                .contract(asDiagonal(MPS.L_tail), idx<1>({0},{0}))
                                .contract(MPS.GA,                 idx<2>({2,3},{0,1}))
                                .shuffle(array3{0,2,1});
                block = block_enlarged;
                break;
            case (Side::R):
                /*! # Right environment contraction
                 *  1--[ GB conj ]--2 0--[LB]--1  0--[      ]
                 *          |                        [      ]
                 *          0                        [      ]
                 *                                   [      ]
                 *          2                        [      ]
                 *          |                        [      ]
                 *   0--[   W    ]--1             2--[ right]
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
                                .contract(W,                  idx<2>({1,2},{1,2}))
                                .contract(asDiagonal(MPS.LB), idx<1>({0},{1}))
                                .contract(MPS.GB,             idx<2>({2,3},{0,2}))
                                .shuffle(array3{0,2,1});
                block = block_enlarged;
                break;
        }
    }

};


#endif //DMRG_CLASS_ENVIRONMENT_H
