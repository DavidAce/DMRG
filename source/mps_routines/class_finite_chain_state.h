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

    int max_sites = 0;                                                 /*!< The maximum length of the chain */
    int num_sweeps = 0;
    int direction  = -1;

public:
    class_finite_chain_state()=default;
    explicit class_finite_chain_state(int max_sites_);

    void do_all_measurements();
    void set_max_sites(int max_sites_);                                        /*!< Sets the maximum length of the chain. */
    void clear();

    bool max_sites_is_set                       = false;
    bool mps_have_been_written_to_hdf5          = false;
    bool mpo_have_been_written_to_hdf5          = false;
    bool env_have_been_written_to_hdf5          = false;
    bool model_params_have_been_written_to_hdf5 = false;


    bool energy_has_been_measured      = false;
    bool variance_has_been_measured    = false;
    bool entropy_has_been_measured     = false;
    bool norm_has_been_measured        = false;
    bool parity_has_been_measured      = false;
    bool everything_has_been_measured  = false;
    bool has_been_measured();
    void set_measured_true();
    void set_measured_false();
    bool has_been_written();
    void set_written_true();
    void set_written_false();



    int    get_sweeps() const ;
    void   set_sweeps(int num_sweeps_) {num_sweeps = num_sweeps_;}
    void   increment_sweeps() {num_sweeps++;}
    int    reset_sweeps();

    void set_positions();
    size_t get_length() const;
    int get_position() const;


    int get_direction() const;
    void flip_direction();
    bool position_is_the_middle() const ;
    bool position_is_the_middle_any_direction() const ;
    bool position_is_the_left_edge() const ;
    bool position_is_the_right_edge() const ;
    bool position_is_any_edge() const ;
    bool position_is_at(int pos)const;
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


    Eigen::Tensor<Scalar,3> & get_G(size_t pos);
    Eigen::Tensor<Scalar,1> & get_L(size_t pos);
    Eigen::Tensor<Scalar,3>   get_A(size_t pos);
    Eigen::Tensor<Scalar,3>   get_B(size_t pos);
    Eigen::Tensor<Scalar,4>   get_theta(size_t pos);
    std::tuple<long,long,long>  get_dims(size_t pos);

    Eigen::Tensor<Scalar,3> & get_GA();
    Eigen::Tensor<Scalar,3> & get_GB();
    Eigen::Tensor<Scalar,1> & get_LA();
    Eigen::Tensor<Scalar,1> & get_LC();
    Eigen::Tensor<Scalar,1> & get_LB();
    Eigen::Tensor<Scalar,3>   get_A();
    Eigen::Tensor<Scalar,3>   get_B();
    Eigen::Tensor<Scalar,4>   get_theta();


    struct Measurements {
        int    length                                   = 0;
        std::vector<int> bond_dimensions                ;
        double norm                                     = 0;
        double energy_mpo, energy_per_site_mpo          = 0;
        double energy_per_site_ham                      = 0;
        double energy_variance_mpo                      = 0;
        double energy_variance_per_site_mpo             = 0;
        double spin_component_sx                                = 0;
        double spin_component_sy                                = 0;
        double spin_component_sz                                = 0;
        std::vector<double> spin_components                       ;
        std::vector<double> entanglement_entropies         ;
    } measurements;


};


#endif //DMRG_CLASS_FINITE_CHAIN_STORAGE_H
