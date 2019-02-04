//
// Created by david on 2019-01-29.
//

#ifndef CLASS_FINITE_CHAIN_STORAGE_H
#define CLASS_FINITE_CHAIN_STORAGE_H
#include <memory>
#include <variant>
#include <complex>
#include <general/nmspc_tensor_extra.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_finite_chain_state.h>
#include <model/class_hamiltonian_base.h>


/*!
 \class class_finite_chain_storage
 \brief Storage class for the finite chain and effective environments during finite DMRG
 The storage is always partitioned into two lists, corresponding to the two halves about the current position.
 The back and front of these lists correspond to the current state of the superblock's left and right parts, respectively,
 while "MPS_C" is the central \f$ \Lambda^A \f$.
 Once full, the particle content of the chain should match the length given in the constant max_length.
*/

class class_finite_chain_state {
public:
    using Scalar = std::complex<double>;

private:

    std::list<class_vidal_mps>                         MPS_L;   /*!< A list of stored \f$ \Lambda^B \Gamma^A...  \f$-tensors. */
    std::list<class_vidal_mps>                         MPS_R;   /*!< A list of stored \f$ \Gamma^B \Lambda^B...  \f$-tensors. */
    Eigen::Tensor<Scalar,1>                            MPS_C;   //Current center bond matrix;
    std::list<class_environment>                       ENV_L;
    std::list<class_environment>                       ENV_R;
    std::list<class_environment_var>                   ENV2_L;
    std::list<class_environment_var>                   ENV2_R;
    std::list<std::shared_ptr<class_hamiltonian_base>> MPO_L;     /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::list<std::shared_ptr<class_hamiltonian_base>> MPO_R;     /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */

    size_t max_sites = 0;                                                 /*!< The maximum length of the chain */
    size_t num_sweeps = 0;
    int direction  = -1;

public:
    class_finite_chain_state()=default;
    class_finite_chain_state(size_t max_sites_):max_sites(max_sites_){};
    void set_max_length(size_t max_length_);                                        /*!< Sets the maximum length of the chain. */


    bool max_length_is_set                          = false;
    bool all_mps_have_been_written_to_hdf5          = false;
    bool all_mpo_have_been_written_to_hdf5          = false;
    bool all_env_have_been_written_to_hdf5          = false;
    bool all_model_params_have_been_written_to_hdf5 = false;

    bool energy_has_been_measured   = false;
    bool variance_has_been_measured = false;
    bool entropy_has_been_measured  = false;
    void set_measured_true();
    void set_measured_false();



    size_t get_sweeps() const ;
    void   increment_sweeps() {num_sweeps++;}
    size_t reset_sweeps();


    size_t get_length() const;
    size_t get_position() const;


    int get_direction() const;
    void flip_direction();
    bool position_is_the_middle() const ;
    bool position_is_the_middle_any_direction() const ;
    bool position_is_the_left_edge() const ;
    bool position_is_the_right_edge() const ;

    const auto & get_MPS_L() const {return std::as_const(MPS_L);}
    const auto & get_MPS_R() const {return std::as_const(MPS_R);}
    const auto & get_MPS_C() const {return std::as_const(MPS_C);}
    const auto & get_MPO_L() const {return std::as_const(MPO_L);}
    const auto & get_MPO_R() const {return std::as_const(MPO_R);}
    const auto & get_ENV_L() const {return std::as_const(ENV_L);}
    const auto & get_ENV_R() const {return std::as_const(ENV_R);}
    const auto & get_ENV2_L()const {return std::as_const(ENV2_L);}
    const auto & get_ENV2_R()const {return std::as_const(ENV2_R);}
    auto & get_MPS_L() {return MPS_L;}
    auto & get_MPS_R() {return MPS_R;}
    auto & get_MPS_C() {return MPS_C;}
    auto & get_MPO_L() {return MPO_L;}
    auto & get_MPO_R() {return MPO_R;}
    auto & get_ENV_L() {return ENV_L;}
    auto & get_ENV_R() {return ENV_R;}
    auto & get_ENV2_L(){return ENV2_L;}
    auto & get_ENV2_R(){return ENV2_R;}

};


#endif //DMRG_CLASS_FINITE_CHAIN_STORAGE_H
