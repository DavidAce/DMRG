//
// Created by david on 2019-01-29.
//

#ifndef CLASS_FINITE_CHAIN_STORAGE_H
#define CLASS_FINITE_CHAIN_STORAGE_H
#include <memory>
#include <complex>
#include <optional>
#include <general/nmspc_tensor_extra.h>
#include <mps_state/class_environment.h>
#include <mps_state/class_mps_2site.h>
#include <mps_state/class_finite_chain_state.h>
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
private:


    int num_sweeps   = 0;
    int direction    = 1;

public:
    using Scalar = std::complex<double>;
    class_finite_chain_state()=default;

    std::list<class_vidal_site>                        MPS_L;   /*!< A list of stored \f$ \Lambda^B \Gamma^A...  \f$-tensors. */
    std::list<class_vidal_site>                        MPS_R;   /*!< A list of stored \f$ \Gamma^B \Lambda^B...  \f$-tensors. */
    Eigen::Tensor<Scalar,1>                            MPS_C;   //Current center bond matrix;
    std::list<class_environment>                       ENV_L;
    std::list<class_environment>                       ENV_R;
    std::list<class_environment_var>                   ENV2_L;
    std::list<class_environment_var>                   ENV2_R;
    std::list<std::shared_ptr<class_hamiltonian_base>> MPO_L;     /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::list<std::shared_ptr<class_hamiltonian_base>> MPO_R;     /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */

    void do_all_measurements();



    int    get_sweeps()                         const ;
    void   set_sweeps(int num_sweeps_) {num_sweeps = num_sweeps_;}
    void   increment_sweeps() {num_sweeps++;}
    int    reset_sweeps();

    void set_positions();
    size_t get_length()                         const;
    size_t get_position()                       const;
    void flip_direction();
    int  get_direction()                        const;
    bool position_is_the_middle()               const;
    bool position_is_the_middle_any_direction() const;
    bool position_is_the_left_edge()            const;
    bool position_is_the_right_edge()           const;
    bool position_is_any_edge()                 const;
    bool position_is_at(size_t pos)             const;
    bool isReal()                               const;


//    std::tuple<long,long,long>  get_dims(size_t pos)                    const;
    const class_vidal_site       & get_MPS(size_t pos)                  const;
          class_vidal_site       & get_MPS(size_t pos);
    const class_hamiltonian_base & get_MPO(size_t pos)                  const;
          class_hamiltonian_base & get_MPO(size_t pos);
    const class_environment      & get_ENVL(size_t pos)                 const;
    const class_environment      & get_ENVR(size_t pos)                 const;
    const class_environment_var  & get_ENV2L(size_t pos)                const;
    const class_environment_var  & get_ENV2R(size_t pos)                const;
    const Eigen::Tensor<Scalar,3> & get_G(size_t pos)                   const;
    const Eigen::Tensor<Scalar,1> & get_L(size_t pos)                   const;
    Eigen::Tensor<Scalar,3> & get_G(size_t pos);
    Eigen::Tensor<Scalar,1> & get_L(size_t pos);
    Eigen::Tensor<Scalar,3>   get_A()                                   const;
    Eigen::Tensor<Scalar,3>   get_B()                                   const;
    Eigen::Tensor<Scalar,3>   get_A(size_t pos)                         const;
    Eigen::Tensor<Scalar,3>   get_B(size_t pos)                         const;
    Eigen::Tensor<Scalar,4>   get_theta()                               const;
    Eigen::Tensor<Scalar,4>   get_theta(size_t pos)                     const;

    //For multisite
    std::list<size_t>      active_sites;
    std::list<size_t>      activate_sites(long threshold);
    Eigen::DSizes<long,3>  active_dimensions() const;

    Eigen::Tensor<Scalar,3>   get_multitheta()    const;
    Eigen::Tensor<Scalar,4>   get_multimpo  ()    const;
    std::pair<std::reference_wrapper<const class_environment>     , std::reference_wrapper<const class_environment>>      get_multienv ()     const;
    std::pair<std::reference_wrapper<const class_environment_var> , std::reference_wrapper<const class_environment_var>>  get_multienv2()     const;

    Eigen::Tensor<Scalar,6>   get_multi_hamiltonian() const;
    Eigen::Tensor<Scalar,6>   get_multi_hamiltonian2() const;
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian_matrix() const;
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> get_multi_hamiltonian2_matrix() const;


    std::vector<double>  truncation_error;

    struct Measurements {
        std::optional<size_t>               length                                  = {};
        std::optional<size_t>               bond_dimension_midchain                 = {};
        std::optional<size_t>               bond_dimension_current                  = {};
        std::optional<std::vector<size_t>>  bond_dimensions                         = {};
        std::optional<double>               norm                                    = {};
        std::optional<double>               energy                                  = {};
        std::optional<double>               energy_per_site                         = {};
        std::optional<double>               energy_variance_mpo                     = {};
        std::optional<double>               energy_variance_per_site                = {};
        std::optional<double>               spin_component_sx                       = {};
        std::optional<double>               spin_component_sy                       = {};
        std::optional<double>               spin_component_sz                       = {};
        std::optional<std::vector<double>>  spin_components                         = {};
        std::optional<double>               entanglement_entropy_midchain           = {};
        std::optional<double>               entanglement_entropy_current            = {};
        std::optional<std::vector<double>>  entanglement_entropies                  = {};
    };
    mutable Measurements measurements;
    void unset_measurements() const;
    void do_all_measurements()const;
};


#endif //DMRG_CLASS_FINITE_CHAIN_STORAGE_H
