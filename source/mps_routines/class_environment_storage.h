
#ifndef DMRG_CLASS_STORAGE_H
#define DMRG_CLASS_STORAGE_H

#include "general/nmspc_tensor_extra.h"
#include "class_environment.h"
#include "sim_parameters/nmspc_model.h"
#include <map>
#include <list>
#include <memory>
#include <variant>

//using namespace Textra;
//using namespace std;
class class_superblock;
class class_hdf5_file;


/*!
 \class class_environment_storage
 \brief Storage class for environments during finite DMRG
 During the infinite-DMRG part of the algorithm, the current state for each chain-length
 is stored in a `class_storage` object, for later use in the finite-DMRG stage.

 Specifically, this class stores objects from the current superblock:
 - The left and right environment blocks `Lblock` and `Rblock`, both type `Tensor3`, and
 - the MPS tensors \f$\Gamma^{A,B}=\f$ `MPS.GA, MPS.GB` and \f$\Lambda^{A,B}=\f$ `MPS.LA, MPS.LB`, required for computing \f$\Theta\f$, following
 the notation used in [Frank Pollmann's notes](http://quantumtensor.pks.mpg.de/wp-content/uploads/2016/06/notes_1.pdf):
 \f$\Theta= \Lambda^B\Gamma^A\Lambda^A\Gamma^B\Lambda^B\f$. <br> Here \f$\Gamma^{A,B}\f$ are a rank-3 tensors of type `Tensor3`  and \f$\Lambda^{A,B}\f$
 are rank-1 tensors of type `Tensor1`.

 All objects are stored and indexed by their position relative to the final chain of length `max_length`.
*/


class class_environment_storage {
public:
    using Scalar = std::complex<double>;
public:
    std::list<std::tuple<Textra::Tensor<Scalar,3>,Textra::Tensor<Scalar,1>, Textra::Tensor<Scalar,3>>>  MPS_L;  /*!< A list of stored \f$ \Gamma^A \Gamma^B...  \f$-tensors with corresponding,  \f$ \Lambda^A \Lambda^B...  \f$ and left environments  indexed by chain position. */
    std::list<std::tuple<Textra::Tensor<Scalar,3>,Textra::Tensor<Scalar,1>, Textra::Tensor<Scalar,3>>>  MPS_R;  /*!< A list of stored \f$ \Gamma^A \Gamma^B...  \f$-tensors with corresponding,  \f$ \Lambda^A \Lambda^B...  \f$ and right environments indexed by chain position. */
    std::list<Textra::Tensor<Scalar,4>>  MPO_L; /*!< A list of stored Hamiltonian MPO tensors,  indexed by chain position. */
    std::list<Textra::Tensor<Scalar,4>>  MPO_R; /*!< A list of stored Hamiltonian MPO tensors,  indexed by chain position. */
    std::list<Textra::Tensor<Scalar,4>>  block_Sq_L; /*!< A list of stored Hamiltonian MPO tensors,  indexed by chain position. */
    std::list<Textra::Tensor<Scalar,4>>  block_Sq_R; /*!< A list of stored Hamiltonian MPO tensors,  indexed by chain position. */

private:
    std::shared_ptr<class_superblock> superblock;
    std::shared_ptr<class_hdf5_file>  hdf5;
public:
    int max_length = 0;                                             /*!< The maximum length of the chain */
    int current_length = 0;
    int position_L;                                                 /*!< The current position of \f$\Gamma^A\f$ w.r.t the full chain. */
    int position_R;                                                 /*!< The current position of \f$\Gamma^B\f$ w.r.t the full chain. */
    bool max_length_is_set = false;
    bool superblock_is_set = false;
    bool hdf5_file_is_set  = false;

    void print_error_and_exit(int error_type);
public:


    class_environment_storage(){};
    explicit class_environment_storage(int max_length_,                                   /*!< The maximum length of the chain. */
                                 std::shared_ptr<class_superblock> superblock_,
                                 std::shared_ptr<class_hdf5_file>  hdf5_
    );

    template<typename T>
    void insert_two_in_the_middle(std::list<T>& list, T& objLeft, T&objRight){
        auto it = list.begin();
        std::advance(it, std::distance(list.begin(), list.end())/2);
        list.insert(it, {objLeft,objRight});
    }



    void set_length(int max_length_);                                       /*!< Sets the maximum length of the chain. */
    void set_superblock(std::shared_ptr<class_superblock> superblock_);     /*!< Sets the pointer to a superblock */
    void set_hdf5_file (std::shared_ptr<class_hdf5_file>  hdf5_);           /*!< Sets the pointer to an hdf5-file for storage */
    int  insert();                        /*!< Store current MPS and environments indexed by their respective positions on the chain. */
    int  load();                                /*!< Load MPS and environments according to current position. */
    void overwrite_MPS();                 /*!< Update the MPS stored at current position.*/
    int  move(int &direction, int &sweep);    /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment. Turn direction around if the edge is reached. */
    void print_storage();                                                   /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
};




#endif //DMRG_CLASS_STORAGE_H
