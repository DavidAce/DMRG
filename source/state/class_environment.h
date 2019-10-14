//
// Created by david on 7/21/17.
//

#ifndef DMRG_CLASS_ENVIRONMENT_H
#define DMRG_CLASS_ENVIRONMENT_H

#include <memory>

#include "general/nmspc_tensor_extra.h"
#include <spdlog/fmt/fmt.h>

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

class class_mps_site;

class class_environment{
public:
    using Scalar = std::complex<double>;
private:
    std::optional<size_t> position;
    bool edge_has_been_set = false;
    void enlarge      (const Eigen::Tensor<Scalar,3>  & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
public:
    std::string side;
    size_t sites = 0;                               /*!< Number of particles that have been contracted into this environment. */
    Eigen::Tensor<Scalar,3> block;                 /*!< The environment block. */
    explicit class_environment(std::string side_):side(side_){};
    class_environment(
            std::string side_,
            int mpsDim,
            int mpoDim)
            :side(side_)
    {
        set_edge_dims(mpsDim,mpoDim);
    }
    bool isReal () const;
    void enlarge      (const class_mps_site & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
    void set_edge_dims(const class_mps_site & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
    void set_edge_dims(const Eigen::Tensor<Scalar,3> & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
    void set_edge_dims(const int mpsDim, const int mpoDim);
    void set_position(const size_t position_){position = position_;}
    size_t get_position() const {
        if(position) {return position.value();}
        else{throw std::runtime_error(fmt::format("Position hasn't been set on environment side {}", side));}
    }

};



class class_environment_var{
public:
    using Scalar = std::complex<double>;
private:
    std::optional<size_t> position;
    bool edge_has_been_set = false;
    void enlarge(const Eigen::Tensor<Scalar,3>  & MPS, const Eigen::Tensor<Scalar,4> &MPO);
public:
    size_t sites = 0;                                      /*!< Number of particles that have been contracted into this left environment. */
    std::string side;
    Eigen::Tensor<Scalar,4> block;                         /*!< The environment block. */
    explicit class_environment_var(std::string side_):side(side_){};
    class_environment_var(
            std::string side_,
            int mpsDim,
            int mpoDim)
            :side(side_)
    {
        set_edge_dims(mpsDim,mpoDim);
    }
    bool isReal () const;
    void enlarge(const class_mps_site & MPS, const Eigen::Tensor<Scalar,4> &MPO);
    void set_edge_dims(const class_mps_site & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
    void set_edge_dims(const Eigen::Tensor<Scalar,3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO);
    void set_edge_dims(const int mpsDim, const int mpoDim);
    void set_position(const size_t position_){position = position_;}
    size_t get_position() const {
        if(position) {return position.value();}
        else{throw std::runtime_error("Position hasn't been set on environment var " + side);}
    }};

#endif //DMRG_CLASS_ENVIRONMENT_H
