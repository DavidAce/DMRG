//
// Created by david on 7/21/17.
//

#ifndef DMRG_CLASS_ENVIRONMENT_H
#define DMRG_CLASS_ENVIRONMENT_H

#include "general/n_tensor_extra.h"
#include "mps_routines/class_mps.h"
#include "sim_parameters/n_model.h"

//using namespace Textra;
//using namespace std;

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
    using Scalar = class_mps::Scalar;
    int size;                                       /*!< Number of particles that have been contracted into this left environment. */
    Textra::Tensor<Scalar,3> block;                         /*!< The environment block. */
    class_environment(){
        size = 0;
        block.resize(1, 1, 3);
        block.setZero();
        if constexpr(side == Side::L){
//            block(0, 0, 0) = 1;
            block(0, 0, 2) = 1;
        }
        else
        if constexpr(side== Side::R){
//            block(0, 0, 2) = 1;
            block(0, 0, 0) = 1;
        }
    };

    void enlarge(const class_mps &MPS, const Textra::Tensor<Scalar,4> &M);
};


template<Side side>
class class_environment_var{
public:
    using Scalar = class_mps::Scalar;
    int size;                                       /*!< Number of particles that have been contracted into this left environment. */
    Textra::Tensor<Scalar,4> block;                         /*!< The environment block. */
    class_environment_var(){
        size = 0;
        block.resize(1, 1, 3, 3) ;
        block.setZero();
        if constexpr(side == Side::L){
            block(0, 0, 2 ,2) = 1;
//            block(0, 0, 2 ,2) = 1;
//            block(0, 0, 2 ,0) = 1;
//            block(0, 0, 2 ,1) = 1;
        }
        else
        if constexpr(side== Side::R){
            block(0, 0, 0, 0) = 1;
//            block(0, 0, 0, 1) = 1;
//            block(0, 0, 0, 2) = 1;
        }
    };

    void enlarge(const class_mps &MPS, const Textra::Tensor<Scalar,4> &M);
};

#endif //DMRG_CLASS_ENVIRONMENT_H
