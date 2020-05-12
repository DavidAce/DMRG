//
// Created by david on 7/21/17.
//

#pragma once

#include <memory>

#include "general/nmspc_tensor_extra.h"
#include <spdlog/fmt/fmt.h>

// using namespace Textra;
// using namespace std;

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
class class_mpo_base;

class class_environment_base {
    public:
    using Scalar = std::complex<double>;

    protected:
    std::optional<size_t> position;
    bool                  edge_has_been_set                                                                 = false;
    virtual void          enlarge(const Eigen::Tensor<Scalar, 3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO) = 0;

    public:
    std::string side;
    size_t      sites = 0; /*!< Number of particles that have been contracted into this environment. */
    explicit class_environment_base(std::string side_, size_t position_);
    explicit class_environment_base(std::string side_, const class_mps_site &MPS, const class_mpo_base &MPO);

    virtual bool isReal() const                                                        = 0;
    virtual bool hasNaN() const                                                        = 0;
    virtual void assertValidity() const                                                = 0;
    virtual void set_edge_dims(const class_mps_site &MPS, const class_mpo_base &MPO) = 0;
    void         set_position(const size_t position_) { position = position_; }
    size_t       get_position() const {
        if(position) {
            return position.value();
        } else {
            throw std::runtime_error(fmt::format("Position hasn't been set on environment side {}", side));
        }
    }
};

class class_environment final : public class_environment_base {
    private:
    void enlarge(const Eigen::Tensor<Scalar, 3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO) final;

    public:
    Eigen::Tensor<Scalar, 3> block; /*!< The environment block. */
    using class_environment_base::class_environment_base;
    explicit class_environment(std::string side_, const class_mps_site &MPS, const class_mpo_base &MPO);
    [[nodiscard]] class_environment enlarge(const class_mps_site &MPS, const class_mpo_base &MPO);

    bool isReal() const final;
    bool hasNaN() const final;
    void assertValidity() const final;
    void set_edge_dims(const class_mps_site &MPS, const class_mpo_base &MPO) final;
};

class class_environment_var final : public class_environment_base {
    private:
    void enlarge(const Eigen::Tensor<Scalar, 3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO) final;

    public:
    Eigen::Tensor<Scalar, 4> block; /*!< The environment block. */
    using class_environment_base::class_environment_base;
    explicit class_environment_var(std::string side_, const class_mps_site &MPS, const class_mpo_base &MPO);
    [[nodiscard]] class_environment_var enlarge(const class_mps_site &MPS, const class_mpo_base &MPO);

    bool isReal() const final;
    bool hasNaN() const final;
    void assertValidity() const final;
    void set_edge_dims(const class_mps_site &MPS, const class_mpo_base &MPO) final;
};

//
//
// class class_environment_var{
// public:
//    using Scalar = std::complex<double>;
// private:
//    std::optional<size_t> position;
//    bool edge_has_been_set = false;
//    void enlarge(const Eigen::Tensor<Scalar,3>  & MPS, const Eigen::Tensor<Scalar,4> &MPO);
//    void set_edge_dims(const class_mps_site & MPS, const class_model_base &MPO);
////    void set_edge_dims(const class_mps_site & MPS, const Eigen::Tensor<Scalar, 4> &MPO);
//    void set_edge_dims(const Eigen::Tensor<Scalar,3> &MPS, const Eigen::Tensor<Scalar, 4> &MPO);
//    void set_edge_dims(const int mpsDim, const int mpoDim);
// public:
//    size_t sites = 0;                                      /*!< Number of particles that have been contracted into this left environment. */
//    std::string side;
//    Eigen::Tensor<Scalar,4> block;                         /*!< The environment block. */
//    explicit class_environment_var(std::string side_, size_t position_):position(position_),side(side_){};
//
////    explicit class_environment_var(std::string side_):side(side_){};
////    class_environment_var(
////            std::string side_,
////            int mpsDim,
////            int mpoDim)
////            :side(side_)
////    {
////        set_edge_dims(mpsDim,mpoDim);
////    }
//    bool isReal () const;
////    void enlarge(const class_mps_site & MPS, const class_model_base &MPO);
////    void enlarge(const class_mps_site & MPS, const Eigen::Tensor<Scalar,4> &MPO);
//    class_environment_var enlarge(const class_mps_site & MPS, const class_model_base &MPO);
//
//
//    void set_position(const size_t position_){position = position_;}
//    size_t get_position() const {
//        if(position) {return position.value();}
//        else{throw std::runtime_error("Position hasn't been set on environment var " + side);}
//    }};
//
