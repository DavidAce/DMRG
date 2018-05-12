
#ifndef DMRG_CLASS_STORAGE_H
#define DMRG_CLASS_STORAGE_H
#include <map>
#include <list>
#include <memory>
#include <variant>
#include <complex>
#include <general/nmspc_tensor_extra.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_hamiltonian.h>
#include <sim_parameters/nmspc_model.h>


class class_superblock;
class class_hdf5_file;


/*!
 \class class_environment_storage
 \brief Storage class for environments during finite DMRG
 During the infinite-DMRG part of the algorithm, the current state for each chain-length
 is stored in a `class_finite_chain_storage` object, for later use in the finite-DMRG stage.

 There are as many elements in storage as there are particles in the current environment, i.e. superblock->environment_size.



 The storage is always partitioned into two lists, corresponding to the two halves about the current position.
 The back and front of these lists correspond to the current state of the superblock's left and right parts, respectively,
 while "LA" is the central \f$ \Lambda^A \f$.

 Once full, the particle content of the chain should match the length given in the constant max_length.
 This can be checked with the following conditions:
 ENV_L.back().size + ENV_R.front().size == max_length - 2
 ENV_L.back().size + ENV_R.front().size == superblock->environment_size
 ENV_L.back().size == superblock->Lblock->size
 ENV_R.back().size == superblock->Rblock->size


*/



class class_finite_chain_storage {
public:
    using Scalar = std::complex<double>;
public:

    std::list<std::tuple<Textra::Tensor<Scalar,1>,Textra::Tensor<Scalar,3>>>  MPS_L;  /*!< A list of stored \f$ \Lambda^B \Gamma^A...  \f$-tensors. */
    std::list<std::tuple<Textra::Tensor<Scalar,3>,Textra::Tensor<Scalar,1>>>  MPS_R;  /*!< A list of stored \f$ \Gamma^B \Lambda^B...  \f$-tensors. */
//    std::listTextra::Tensor<Scalar,1>> MPS_C; /*!< A list of stored \f$ \Lambda^A\f$-tensors, i.e. center bond matrices. */
    Textra::Tensor<Scalar,1> LA;  //Current center bond matrix;
    std::list<class_environment> ENV_L;
    std::list<class_environment> ENV_R;
    std::list<class_environment_var> ENV2_L;
    std::list<class_environment_var> ENV2_R;
    std::list<class_hamiltonian> MPO_L;   /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */
    std::list<class_hamiltonian> MPO_R;   /*!< A list of stored Hamiltonian MPO tensors,indexed by chain position. */



private:
    std::shared_ptr<class_superblock> superblock;
    std::shared_ptr<class_hdf5_file>  hdf5;
public:
    int max_length = 0;                                             /*!< The maximum length of the chain */
//    int current_length = 0;
//    int position_L;                                                 /*!< The current position of \f$\Gamma^A\f$ w.r.t the full chain. */
//    int position_R;                                                 /*!< The current position of \f$\Gamma^B\f$ w.r.t the full chain. */
    long current_length = 0;
    bool max_length_is_set = false;
    bool superblock_is_set = false;
    bool hdf5_file_is_set  = false;

    void print_error_and_exit(int error_type);
public:


    class_finite_chain_storage()=default;
    explicit class_finite_chain_storage(int max_length_,                                   /*!< The maximum length of the chain. */
                                 std::shared_ptr<class_superblock> superblock_,
                                 std::shared_ptr<class_hdf5_file>  hdf5_
    );

    template<typename T>
    void insert_two_in_the_middle(std::list<T>& list, T& objLeft, T&objRight){
        auto it = list.begin();
        std::advance(it, std::distance(list.begin(), list.end())/2);
        list.insert(it, {objLeft,objRight});
    }



    void set_max_length(int max_length_);                                       /*!< Sets the maximum length of the chain. */
    void set_superblock(std::shared_ptr<class_superblock> superblock_);     /*!< Sets the pointer to a superblock */
    void set_hdf5_file (std::shared_ptr<class_hdf5_file>  hdf5_);           /*!< Sets the pointer to an hdf5-file for storage */
    void update_current_length();
    int  insert();                        /*!< Store current MPS and environments indexed by their respective positions on the chain. */
//    int  insert_edges();
    int  load();                                /*!< Load MPS and environments according to current position. */
    void overwrite_MPS();                 /*!< Update the MPS stored at current position.*/
    int  move(int &direction, int &sweep);    /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
    void print_storage();                                                   /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
    void print_storage_compact();                                                   /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */

    void print_hamiltonian_energies();

};




#endif //DMRG_CLASS_STORAGE_H
