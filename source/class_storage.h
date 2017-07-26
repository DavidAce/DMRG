
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
 - the MPS tensors \f$\Gamma^{A,B}=\f$ `MPS.G.A, MPS.G.B` and \f$\Lambda^{A,B}=\f$ `MPS.L.A, MPS.L.B`, required for computing \f$\Theta\f$, following
 the notation used in [Frank Pollmann's notes](http://quantumtensor.pks.mpg.de/wp-content/uploads/2016/06/notes_1.pdf):
 \f$\Theta= \Lambda^B\Gamma^A\Lambda^A\Gamma^B\Lambda^B.\f$.  Here \f$\Gamma^{A,B}\f$ are a rank-3 tensors of type `Textra::Tensor3`  and \f$\Lambda^{A,B}\f$
 are rank-1 tensors of `Textra::Tensor1`.

*/


class class_storage {
private:
    template <typename list_type, typename inType>
    void insert_middle(list_type &target_list, const inType &elem){
        typename list_type::iterator it = target_list.begin();
        std::advance(it, std::distance(target_list.begin(), target_list.end())/2);
        target_list.insert(it, elem);
    }
    template <typename list_type, typename inType>
    void replace(list_type &target_list, const inType &elem, const int at){
        typename list_type::iterator it = target_list.begin();
        std::advance(it, std::distance(target_list.begin(), target_list.begin()+at));
        *it = elem;
    }
public:

    std::map<int, Tensor3> G_list;                               /*!< Detailed description after the member */
    std::map<int, Tensor1> L_list;                               /*!< Detailed description after the member */
    std::map<int,class_environment_L> Lblock_list;               /*!< Detailed description after the member */
    std::map<int,class_environment_R> Rblock_list;               /*!< Detailed description after the member */

    const int max_length;
    int position_L = 0;
    int position_R = max_length - 1;

    class_storage(int L):max_length(L){
    };

    void store_insert(const class_superblock &superblock);
    void load(class_superblock &superblock);
    void overwrite(const class_superblock &superblock);
    void move(class_superblock &superblock, const int direction);
    void print_storage();
};


#endif //FINITE_DMRG_EIGEN_CLASS_STORAGE_H
