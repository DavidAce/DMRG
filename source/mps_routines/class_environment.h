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

class class_mps_2site;

class class_environment{
private:
    long position;
public:
    using Scalar = std::complex<double>;
    std::string side;
    unsigned long size = 0;                        /*!< Number of particles that have been contracted into this left environment. */
    Eigen::Tensor<Scalar,3> block;                 /*!< The environment block. */
    explicit class_environment(std::string side_):side(std::move(side_)){size = 0;};
    void enlarge(const std::unique_ptr<class_mps_2site> &MPS, const Eigen::Tensor<Scalar,4> &M);
    void set_edge_dims(const std::unique_ptr<class_mps_2site> &MPS, const Eigen::Tensor<Scalar, 4> &M);
    void set_position(const long position_){position = position_;}
    auto get_position() const {return position;}
};



class class_environment_var{
private:
    long position;
public:
    using Scalar = std::complex<double>;
    unsigned long size;                                       /*!< Number of particles that have been contracted into this left environment. */
    std::string side;
    Eigen::Tensor<Scalar,4> block;                         /*!< The environment block. */
    explicit class_environment_var(std::string side_):side(std::move(side_)){size = 0;};
    void enlarge(const std::unique_ptr<class_mps_2site> &MPS, const Eigen::Tensor<Scalar,4> &M);
    void set_edge_dims(const std::unique_ptr<class_mps_2site> &MPS, const Eigen::Tensor<Scalar, 4> &M);
    void set_position(const long position_){position = position_;}
    auto get_position() const {return position;}
};

#endif //DMRG_CLASS_ENVIRONMENT_H
