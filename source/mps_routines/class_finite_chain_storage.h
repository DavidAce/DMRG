//
// Created by david on 2019-01-29.
//

#ifndef CLASS_FINITE_CHAIN_STORAGE_H
#define CLASS_FINITE_CHAIN_STORAGE_H


class class_finite_chain_storage {
public:
    using Scalar = std::complex<double>;

private:

    std::list<class_vidal_mps>  MPS_L;                                   /*!< A list of stored \f$ \Lambda^B \Gamma^A...  \f$-tensors. */
    std::list<class_vidal_mps>  MPS_R;                                   /*!< A list of stored \f$ \Gamma^B \Lambda^B...  \f$-tensors. */
    Eigen::Tensor<Scalar,1>     MPS_C;                                   //Current center bond matrix;
    std::list<class_environment> ENV_L;
    std::list<class_environment> ENV_R;
    std::list<class_environment_var> ENV2_L;
    std::list<class_environment_var> ENV2_R;
    std::list<std::unique_ptr<class_hamiltonian_base>> MPO_L;            /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::list<std::unique_ptr<class_hamiltonian_base>> MPO_R;            /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */

    unsigned long max_length = 0;                                                 /*!< The maximum length of the chain */

public:
    class_finite_chain_storage()=default;
    explicit class_finite_chain_storage(int max_length_):max_length(max_length_);                        /*!< The maximum length of the chain. */

    void set_max_length(int max_length_);                                        /*!< Sets the maximum length of the chain. */
    bool max_length_is_set         = false;

    const auto & get_MPS_L() const {return std::as_const(MPS_L);}
    const auto & get_MPS_R() const {return std::as_const(MPS_R);}
    const auto & get_MPS_C() const {return std::as_const(MPS_C);}
    const auto & get_MPO_L() const {return std::as_const(MPO_L);}
    const auto & get_MPO_R() const {return std::as_const(MPO_R);}
    const auto & get_ENV_L() const {return std::as_const(ENV_L);}
    const auto & get_ENV_R() const {return std::as_const(ENV_R);}
    const auto & get_ENV2_L()const {return std::as_const(ENV2_L);}
    const auto & get_ENV2_R()const {return std::as_const(ENV2_R);}
    auto & ref_MPS_L() {return MPS_L;}
    auto & ref_MPS_R() {return MPS_R;}
    auto & ref_MPS_C() {return MPS_C;}
    auto & ref_MPO_L() {return MPO_L;}
    auto & ref_MPO_R() {return MPO_R;}
    auto & ref_ENV_L() {return ENV_L;}
    auto & ref_ENV_R() {return ENV_R;}
    auto & ref_ENV2_L(){return ENV2_L;}
    auto & ref_ENV2_R(){return ENV2_R;}

};


#endif //DMRG_CLASS_FINITE_CHAIN_STORAGE_H
