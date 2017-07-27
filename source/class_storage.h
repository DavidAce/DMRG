
#ifndef FINITE_DMRG_EIGEN_CLASS_STORAGE_H
#define FINITE_DMRG_EIGEN_CLASS_STORAGE_H

#include <n_tensor_extra.h>
#include <map>

using namespace Textra;
using namespace std;
class class_superblock;
class class_environment_L;
class class_environment_R;

/*!
 # Storage Class

 During the infinite-DMRG part of the algorithm, the current state for each chain-length
 is stored in a `class_storage` object, for later use in the finite-DMRG stage.

 Specifically, this class stores objects from the current superblock:
 - The left and right environment blocks `Lblock` and `Rblock`, both type `Textra::Tensor3`, and
 - the MPS tensors \f$\Gamma^{A,B}=\f$ `MPS.GA, MPS.GB` and \f$\Lambda^{A,B}=\f$ `MPS.LA, MPS.LB`, required for computing \f$\Theta\f$, following
 the notation used in [Frank Pollmann's notes](http://quantumtensor.pks.mpg.de/wp-content/uploads/2016/06/notes_1.pdf):
 \f$\Theta= \Lambda^B\Gamma^A\Lambda^A\Gamma^B\Lambda^B\f$. <br> Here \f$\Gamma^{A,B}\f$ are a rank-3 tensors of type `Textra::Tensor3`  and \f$\Lambda^{A,B}\f$
 are rank-1 tensors of type `Textra::Tensor1`.

 All objects are stored and indexed by their position relative to the final chain of length `max_length`.
*/


class class_storage {
public:

    std::map<int, Tensor3> G_list;                                  /*!< A list of stored \f$\Gamma\f$-tensors,  indexed by chain position. */
    std::map<int, Tensor1> L_list;                                  /*!< A list of stored \f$\Lambda\f$-tensors, indexed by chain position. */
    std::map<int,class_environment_L> Lblock_list;                  /*!< A list of stored Left block environments,  indexed by chain position. */
    std::map<int,class_environment_R> Rblock_list;                  /*!< A list of stored Right block environments, indexed by chain position. */

    const int max_length;                                           /*!< The maximum length of the chain */
    int position_L = 0;                                             /*!< The current position of \f$\Gamma^A\f$ w.r.t the full chain. */
    int position_R = max_length - 1;                                /*!< The current position of \f$\Gamma^B\f$ w.r.t the full chain. */

    class_storage(int max_length_                                   /*!< The maximum length of the chain. */
                  ):max_length(max_length_){};

    void store_insert(const class_superblock &superblock);          /*!< Store current MPS and environments indexed by their respective positions on the chain. */
    void load(class_superblock &superblock);                        /*!< Load MPS and environments according to current position. */
    void overwrite_MPS(const class_superblock &superblock);         /*!< Update the MPS stored at current position.*/
    void move(class_superblock &superblock, const int direction);   /*!< Move current position to the left (`direction=1`) or right (`direction=-1`), and store the **newly enlarged** environment.  */
    void print_storage();                                           /*!< Print the tensor dimensions for all \f$\Gamma\f$-tensors. */
};


#endif //FINITE_DMRG_EIGEN_CLASS_STORAGE_H
