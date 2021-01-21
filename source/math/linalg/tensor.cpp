#include "tensor.h"
#include <general/nmspc_tensor_extra.h>


template<typename Scalar, auto rank>
Eigen::Tensor<Scalar, rank> linalg::tensor::mirror(const Eigen::Tensor<Scalar, rank> &tensor){
    /*
     Returns a mirrored tensor

     Example: Starting with A
            0 1 2 3
            | | | |
           [  A   ]
           | | | |
           4 5 6 7

     returns

            3 2 1 0
            | | | |
           [  A   ]
           | | | |
           7 6 5 4

     This is useful for compaitibility with kronecker products which gives results indexed right to left.

    */
    if constexpr(rank <= 2) return tensor;
    else{
        std::array<Eigen::Index, rank> shf_idx{};
        for(size_t i = 0; i < static_cast<size_t>(rank); i++) {
            shf_idx[i] = static_cast<Eigen::Index>(i);
        }
        std::reverse(shf_idx.begin(), shf_idx.begin()+rank/2);
        std::reverse(shf_idx.begin()+rank/2, shf_idx.end());
        return tensor.shuffle(shf_idx);
    }
}

template Eigen::Tensor<linalg::cplx, 0> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 0> &tensor);
template Eigen::Tensor<linalg::cplx, 1> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 1> &tensor);
template Eigen::Tensor<linalg::cplx, 2> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 2> &tensor);
template Eigen::Tensor<linalg::cplx, 3> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 3> &tensor);
template Eigen::Tensor<linalg::cplx, 4> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 4> &tensor);
template Eigen::Tensor<linalg::cplx, 5> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 5> &tensor);
template Eigen::Tensor<linalg::cplx, 6> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 6> &tensor);
template Eigen::Tensor<linalg::cplx, 7> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 7> &tensor);
template Eigen::Tensor<linalg::cplx, 8> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 8> &tensor);
template Eigen::Tensor<linalg::cplx, 9> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 9> &tensor);
template Eigen::Tensor<linalg::cplx, 10> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 10> &tensor);
template Eigen::Tensor<linalg::cplx, 11> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 11> &tensor);
template Eigen::Tensor<linalg::cplx, 12> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 12> &tensor);
template Eigen::Tensor<linalg::cplx, 13> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 13> &tensor);
template Eigen::Tensor<linalg::cplx, 14> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 14> &tensor);
template Eigen::Tensor<linalg::cplx, 15> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 15> &tensor);
template Eigen::Tensor<linalg::cplx, 16> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 16> &tensor);
template Eigen::Tensor<linalg::cplx, 17> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 17> &tensor);
template Eigen::Tensor<linalg::cplx, 18> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 18> &tensor);
template Eigen::Tensor<linalg::cplx, 19> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 19> &tensor);
template Eigen::Tensor<linalg::cplx, 20> linalg::tensor::mirror(const Eigen::Tensor<linalg::cplx, 20> &tensor);
template Eigen::Tensor<linalg::real, 0> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 0> &tensor);
template Eigen::Tensor<linalg::real, 1> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 1> &tensor);
template Eigen::Tensor<linalg::real, 2> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 2> &tensor);
template Eigen::Tensor<linalg::real, 3> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 3> &tensor);
template Eigen::Tensor<linalg::real, 4> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 4> &tensor);
template Eigen::Tensor<linalg::real, 5> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 5> &tensor);
template Eigen::Tensor<linalg::real, 6> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 6> &tensor);
template Eigen::Tensor<linalg::real, 7> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 7> &tensor);
template Eigen::Tensor<linalg::real, 8> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 8> &tensor);
template Eigen::Tensor<linalg::real, 9> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 9> &tensor);
template Eigen::Tensor<linalg::real, 10> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 10> &tensor);
template Eigen::Tensor<linalg::real, 11> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 11> &tensor);
template Eigen::Tensor<linalg::real, 12> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 12> &tensor);
template Eigen::Tensor<linalg::real, 13> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 13> &tensor);
template Eigen::Tensor<linalg::real, 14> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 14> &tensor);
template Eigen::Tensor<linalg::real, 15> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 15> &tensor);
template Eigen::Tensor<linalg::real, 16> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 16> &tensor);
template Eigen::Tensor<linalg::real, 17> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 17> &tensor);
template Eigen::Tensor<linalg::real, 18> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 18> &tensor);
template Eigen::Tensor<linalg::real, 19> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 19> &tensor);
template Eigen::Tensor<linalg::real, 20> linalg::tensor::mirror(const Eigen::Tensor<linalg::real, 20> &tensor);






//template<typename Scalar, auto rank, auto npair>
//Eigen::Tensor<Scalar, rank-2*npair> linalg::tensor::trace(const Eigen::Tensor<Scalar, rank> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, npair> & idx_pair, bool mirror) {
//    /*
//     * Returns the partial trace of a tensor
//     * Note that a the tensor given here may be mirrored!
//     */
//    static_assert(rank >= 2*npair, "Rank must be large enough");
//    if constexpr(npair == 1){
//        auto idx_pair_r = idx_pair[0];
//        if(mirror) idx_pair_r = mirror_idx_pair<rank>(idx_pair[0]); // Find the corresponding indices if this tensor was mirrored
//
//        // Collect indices and dimensions traced
//        auto idx1 = static_cast<size_t>(idx_pair_r.first);
//        auto idx2 = static_cast<size_t>(idx_pair_r.second);
//        std::array<size_t,2> idx_tr{idx1,idx2};
//        std::array<Eigen::Index,2> dim_tr{tensor.dimension(idx1), tensor.dimension(idx2)};
//
//        if(dim_tr[0] != dim_tr[1]) throw std::runtime_error("Traced dimensions must be equal size");
//        Eigen::Tensor<Scalar,1> id(dim_tr[0]);
//        id.setConstant(1.0);
//        Eigen::Tensor<Scalar,rank-2> result = Textra::asDiagonal(id).contract(tensor, Textra::idx({1ul,0ul},{idx1,idx2}));
//        return result;
//    }else if constexpr(npair == 2){
//        std::array<long,2> pair1{idx_pair[1].first, idx_pair[1].second};
//        std::array<long,2> pair0{idx_pair[0].first, idx_pair[0].second};
//        pair0[0] -= std::count_if(pair1.begin(),pair1.end(),[&pair0](auto i){return i < pair0[0];});
//        pair0[1] -= std::count_if(pair1.begin(),pair1.end(),[&pair0](auto i){return i < pair0[1];});
//        auto res1 = linalg::tensor::trace(tensor,Textra::idx({pair1[0]},{pair1[1]}));
//        return linalg::tensor::trace(res1,Textra::idx({pair0[0]},{pair0[1]}));
//    }else if constexpr(npair == 3){
//        std::array<long,2> pair2{idx_pair[2].first,idx_pair[2].second};
//        std::array<long,2> pair1{idx_pair[1].first,idx_pair[1].second};
//        std::array<long,2> pair0{idx_pair[0].first,idx_pair[0].second};
//        pair1[0] -= std::count_if(pair2.begin(),pair2.end(),[&pair1](auto i){return i < pair1[0];});
//        pair1[1] -= std::count_if(pair2.begin(),pair2.end(),[&pair1](auto i){return i < pair1[1];});
//        pair0[0] -= std::count_if(pair2.begin(),pair2.end(),[&pair0](auto i){return i < pair0[0];});
//        pair0[1] -= std::count_if(pair2.begin(),pair2.end(),[&pair0](auto i){return i < pair0[1];});
//        auto res = linalg::tensor::trace(tensor,Textra::idx({pair2[0]},{pair2[1]}));
//        return linalg::tensor::trace(tensor,Textra::idx({pair0[0], pair1[0]}, {pair0[1], pair1[1]}));
//    }else
//        throw std::runtime_error("Trace not implemented");
//}
//
//template Eigen::Tensor<linalg::cplx, 0> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 2> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 2> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 4> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 4> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 6> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 6> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 8> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 8> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 10> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 0> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 4> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  2> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 2> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 6> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  2> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 4> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 8> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  2> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 6> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 10> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 2> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 8> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 12> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 2> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 0> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 6> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  3> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 2> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 8> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  3> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 4> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 10> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 3> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 6> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 12> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 3> & idx_pair, bool);
//template Eigen::Tensor<linalg::cplx, 8> linalg::tensor::trace(const Eigen::Tensor<linalg::cplx, 14> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 3> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 0> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 2> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 2> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 4> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 4> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 6> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 6> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 8> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 8> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 10> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 1> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 0> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 4> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  2> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 2> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 6> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  2> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 4> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 8> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  2> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 6> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 10> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 2> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 8> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 12> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 2> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 0> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 6> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  3> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 2> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 8> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>,  3> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 4> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 10> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 3> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 6> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 12> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 3> & idx_pair, bool);
//template Eigen::Tensor<linalg::real, 8> linalg::tensor::trace(const Eigen::Tensor<linalg::real, 14> &tensor, const std::array<Eigen::IndexPair<Eigen::Index>, 3> & idx_pair, bool);