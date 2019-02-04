//
// Created by david on 2019-01-29.
//

#ifndef DMRG_CLASS_FINITE_CHAIN_OPS_H
#define DMRG_CLASS_FINITE_CHAIN_OPS_H
#include <memory>
#include <general/nmspc_tensor_extra.h>
class class_finite_chain_state;
class class_superblock;

class class_finite_chain_ops{
public:
    using Scalar         = std::complex<double>;
    using state_ptr      = std::shared_ptr<class_finite_chain_state>;
    using superblock_ptr = std::shared_ptr<class_superblock>;
    // Functions that operate on the whole chain
    static void apply_mpo(state_ptr state,const Eigen::Tensor<Scalar,4> mpo, const Eigen::Tensor<Scalar,3> Ledge, const Eigen::Tensor<Scalar,3> Redge);
    static void normalize_chain(state_ptr state);
    static void refresh_environments(state_ptr state);
    static void refresh_superblock(state_ptr state, superblock_ptr superblock);
    static void set_parity_projected_mps(state_ptr state,const Eigen::MatrixXcd paulimatrix);
};



#endif //DMRG_CLASS_FINITE_CHAIN_OPS_H
