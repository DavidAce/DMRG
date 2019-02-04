
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
#include <mps_routines/class_finite_chain_state.h>

#include <sim_parameters/nmspc_model.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <IO/class_hdf5_file.h>
#include <iostream>
#include <model/class_hamiltonian_base.h>
class class_superblock;
class class_finite_chain_state;




/*!
 \class class_finite_chain_sweeper
 \brief Sweeper class for the finite chain and effective environments during finite DMRG
 The storage is always partitioned into two lists, corresponding to the two halves about the current position.
 The back and front of these lists correspond to the current state of the superblock's left and right parts, respectively,
 while "MPS_C" is the central \f$ \Lambda^A \f$.
 Once full, the particle content of the chain should match the length given in the constant max_length.
*/





class class_finite_chain_sweeper {
private:
    using Scalar         = std::complex<double>;
    using state_ptr      = std::shared_ptr<class_finite_chain_state>;
    using superblock_ptr = std::shared_ptr<class_superblock>;

public:

    void print_error_and_exit(int error_type);

    static int  insert_superblock_to_chain(state_ptr state, superblock_ptr superblock);                                /*!< Store current MPS and environments indexed by their respective positions on the chain. */
//    int  load();                                  /*!< Load MPS and environments according to current position. */

    static void copy_superblock_to_chain(state_ptr, superblock_ptr);
    static void copy_superblock_mps_to_chain(state_ptr, superblock_ptr);                   /*!< Update the MPS stored at current position.*/
    static void copy_superblock_mpo_to_chain(state_ptr, superblock_ptr);                   /*!< Update the MPO stored at current position.*/
    static void copy_superblock_env_to_chain(state_ptr, superblock_ptr);                   /*!< Update the ENV stored at current position.*/
    static int  move(state_ptr, superblock_ptr);                                  /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */



    static void print_state(state_ptr);                                                        /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    static void print_state_compact(state_ptr);                                                /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    static void print_hamiltonians(state_ptr);


};




#endif //DMRG_CLASS_STORAGE_H
