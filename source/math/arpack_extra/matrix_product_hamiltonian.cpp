#include "matrix_product_hamiltonian.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>

template<typename T>
void DenseHamiltonianProduct<T>::MultAx(T* mps_in_, T* mps_out_) {
    if(counter == 0) omp = std::make_shared<OMP>(num_threads);
    t_mul.tic();
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       envL_map(Lblock,shape_mps[1], shape_mps[1], shape_mpo[0] );
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       envR_map(Rblock,shape_mps[2], shape_mps[2], shape_mpo[1] );
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       mpo_map    (mpo, shape_mpo);
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       mps_map_in  (mps_in_, shape_mps);
    Eigen::TensorMap<Eigen::Tensor<T, 3>>             mps_map_out (mps_out_, shape_mps);

    mps_map_out.device(omp->dev) =
        envL_map
            .contract(mps_map_in, Textra::idx({0},{1}))
            .contract(mpo_map,    Textra::idx({1,2},{0,2}))
            .contract(envR_map,   Textra::idx({1,2},{0,2}))
            .shuffle(Textra::array3{1,0,2});

    //Best yet! I have shown this to be the fastest contraction ordering
//    theta_out.device(omp->dev) = Lblock_map
//        .contract(theta_in,    Textra::idx({0},{1}))
//        .contract(HA_map ,     Textra::idx({1,2},{0,2}))//  idx({1,2,3},{0,4,5}))
//        .contract(HB_map ,     Textra::idx({3,1},{0,2}))//  idx({1,2,3},{0,4,5}))
//        .contract(Rblock_map,  Textra::idx({1,3},{0,2}))
//        .shuffle(Textra::array4{1,0,2,3});
    counter++;
    t_mul.toc();
}

//Explicit instantiations
template class DenseHamiltonianProduct<double>;
template class DenseHamiltonianProduct<std::complex<double>>;

