//
// Created by david on 7/21/17.
//

#ifndef DMRG_CLASS_ENVIRONMENT_H
#define DMRG_CLASS_ENVIRONMENT_H

#include <memory>
#include "general/nmspc_tensor_extra.h"
#include "sim_parameters/nmspc_model.h"

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

class class_mps;

class class_environment{
public:
    using Scalar = std::complex<double>;
    std::string side;
    int size;                                       /*!< Number of particles that have been contracted into this left environment. */
    Textra::Tensor<Scalar,3> block;                 /*!< The environment block. */
    explicit class_environment(std::string side_):side(std::move(side_)){
        size = 0;
        block.resize(1, 1, 3);
        block.setZero();
        if (side == "L"){
            block(0, 0, 2) = 1;
        }
        else
        if (side == "R"){
            block(0, 0, 0) = 1;
        }
        else{
            std::cerr << "Pass strings L or R to initialize environment" << std::endl;
            exit(1);
        }
    };

    void enlarge(const std::shared_ptr<class_mps> &MPS, const Textra::Tensor<Scalar,4> &M);
    void set_edge_dims(const std::shared_ptr<class_mps> &MPS, const Textra::Tensor<Scalar, 4> &M);
};



class class_environment_var{
public:
    using Scalar = std::complex<double>;
    int size;                                       /*!< Number of particles that have been contracted into this left environment. */
    std::string side;

    Textra::Tensor<Scalar,4> block;                         /*!< The environment block. */
    explicit class_environment_var(std::string side_):side(std::move(side_)){
        size = 0;
        block.resize(1, 1, 3, 3) ;
        block.setZero();
        if (side == "L"){
            block(0, 0, 2 ,2) = 1;
        }
        else
        if (side == "R"){
            block(0, 0, 0, 0) = 1;
        }
        else{
            std::cerr << "Pass strings L or R to initialize environment" << std::endl;
            exit(1);
        }
    };
    void enlarge(const std::shared_ptr<class_mps> &MPS, const Textra::Tensor<Scalar,4> &M);
    void set_edge_dims(const std::shared_ptr<class_mps> &MPS, const Textra::Tensor<Scalar, 4> &M);
};

#endif //DMRG_CLASS_ENVIRONMENT_H
