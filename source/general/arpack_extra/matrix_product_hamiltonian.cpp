//
// Created by david on 2018-10-30.
//

#include "matrix_product_hamiltonian.h"
#include <general/nmspc_tensor_extra.h>

template<class T>
void DenseHamiltonianProduct<T>::MultMv(T* theta_in_, T* theta_out_) {
    t_mul.tic();
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       Lblock_map(Lblock,shape_theta4[1], shape_theta4[1], shape_mpo4[0] );
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       Rblock_map(Rblock,shape_theta4[3], shape_theta4[3], shape_mpo4[1] );
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       HA_map    (HA, shape_mpo4);
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       HB_map    (HB, shape_mpo4);
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       theta_in  (theta_in_, shape_theta4);
    Eigen::TensorMap<Eigen::Tensor<T, 4>>             theta_out (theta_out_, shape_theta4);



    //Best yet! The sparcity of the effective hamiltonian (Lblock HA HB Rblock) is about 58% nonzeros.
    //L have shown this to be the fastest contraction ordering
    theta_out = Lblock_map
            .contract(theta_in,    Textra::idx({0},{1}))
            .contract(HA_map ,     Textra::idx({1,2},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(HB_map ,     Textra::idx({3,1},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(Rblock_map,  Textra::idx({1,3},{0,2}))
            .shuffle(Textra::array4{1,0,2,3});
    counter++;
    t_mul.toc();
}

template class DenseHamiltonianProduct<std::complex<double>>;
template class DenseHamiltonianProduct<double>;

