
#ifndef DMRG_CLASS_STORAGE_H
#define DMRG_CLASS_STORAGE_H
#include <map>
#include <list>
#include <memory>
#include <variant>
#include <complex>
#include <general/nmspc_tensor_extra.h>
#include <mps_routines/class_environment.h>
//#include <mps_routines/class_hamiltonian.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_finite_chain_storage.h>

#include <sim_parameters/nmspc_model.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_file.h>
#include <iostream>
#include <model/class_hamiltonian_base.h>
class class_superblock;





/*!
 \class class_finite_chain
 \brief Storage class for the finite chain and effective environments during finite DMRG
 The storage is always partitioned into two lists, corresponding to the two halves about the current position.
 The back and front of these lists correspond to the current state of the superblock's left and right parts, respectively,
 while "MPS_C" is the central \f$ \Lambda^A \f$.
 Once full, the particle content of the chain should match the length given in the constant max_length.
*/



class class_finite_chain {
public:
    using Scalar = std::complex<double>;


private:

    class_finite_chain_storage storage;

    std::list<class_vidal_mps>  MPS_L;                                   /*!< A list of stored \f$ \Lambda^B \Gamma^A...  \f$-tensors. */
    std::list<class_vidal_mps>  MPS_R;                                   /*!< A list of stored \f$ \Gamma^B \Lambda^B...  \f$-tensors. */
    Eigen::Tensor<Scalar,1>     MPS_C;                                   //Current center bond matrix;
    std::list<class_environment> ENV_L;
    std::list<class_environment> ENV_R;
    std::list<class_environment_var> ENV2_L;
    std::list<class_environment_var> ENV2_R;
    std::list<std::unique_ptr<class_hamiltonian_base>> MPO_L;            /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::list<std::unique_ptr<class_hamiltonian_base>> MPO_R;            /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */



    std::shared_ptr<class_superblock> superblock;
    std::shared_ptr<class_hdf5_file>  hdf5;
    SimulationType sim_type;
    std::string    sim_name;

    bool max_length_is_set         = false;
    bool superblock_is_set         = false;
    bool hdf5_file_is_set          = false;

    bool full_mps_has_been_written           = false;
    bool full_mpo_has_been_written           = false;
    bool hamiltonian_params_has_been_written = false;
    int direction = -1;
    int sweeps    = 0;
    unsigned long max_length = 0;                                                 /*!< The maximum length of the chain */
//    unsigned long current_length = 0;

public:

    void print_error_and_exit(int error_type);
    class_finite_chain()=default;
    explicit class_finite_chain(int max_length_,                         /*!< The maximum length of the chain. */
                                 std::shared_ptr<class_superblock> superblock_,
                                 std::shared_ptr<class_hdf5_file>  hdf5_,
                                 SimulationType sim_type,
                                 std::string    sim_name
    );

    void set_max_length(int max_length_);                                        /*!< Sets the maximum length of the chain. */
    void set_superblock(std::shared_ptr<class_superblock> superblock_);          /*!< Sets the pointer to a superblock */
    void update_current_length();
    int  insert();                                                               /*!< Store current MPS and environments indexed by their respective positions on the chain. */
    int  load();                                                                 /*!< Load MPS and environments according to current position. */
    void overwrite_local_MPS();                                                  /*!< Update the MPS stored at current position.*/
    void overwrite_local_MPO();                                                  /*!< Update the MPO stored at current position.*/
    void overwrite_local_ENV();                                                  /*!< Update the ENV stored at current position.*/
    int  move();                                                                 /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */

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

    void print_storage();                                                        /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    void print_storage_compact();                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    void print_hamiltonians();
    int reset_sweeps();
    int get_direction() const;
    int get_position() const;
    int get_length() const;
    int get_sweeps() const;
    bool position_is_the_middle();
    bool position_is_the_middle_any_direction();
    bool position_is_the_left_edge();
    bool position_is_the_right_edge();


    // Functions that operate on the whole chain
    void apply_mpo(const Eigen::Tensor<Scalar,4> mpo, const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge);
    void normalize_chain();
    void refresh_environments();
    void refresh_superblock();
    void set_parity_projected_mps(const Eigen::MatrixXcd paulimatrix);
    void set_parity_projected_mps2(const Eigen::MatrixXcd paulimatrix);
    void set_parity_projected_mps3(const Eigen::MatrixXcd paulimatrix);


    // Functions relating to HDF5 file storage
    void write_all_to_hdf5();
    void write_bond_matrices_to_hdf5();
    void write_2site_mps_to_hdf5();
    void write_2site_mpo_to_hdf5();
    void write_2site_env_to_hdf5();
    void write_2site_env2_to_hdf5();
    void write_full_mps_to_hdf5();
    void write_full_mpo_to_hdf5();
    void write_hamiltonian_params_to_hdf5();
    void write_entanglement_to_hdf5();
    void set_hdf5_file (std::shared_ptr<class_hdf5_file>  hdf5_);                /*!< Sets the pointer to an hdf5-file for storage */



};




#endif //DMRG_CLASS_STORAGE_H
