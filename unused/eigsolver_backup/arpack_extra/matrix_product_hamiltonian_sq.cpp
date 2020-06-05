#include "matrix_product_hamiltonian_sq.h"
#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_tensor_omp.h>

template<typename T>
void DenseHamiltonianSqProduct<T>::MultAx(T* theta_in_, T* theta_out_) {
    if(counter == 0) omp = std::make_shared<OMP>(num_threads);
    t_mul.tic();
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       env2L(env2L_ptr,shape_theta[1], shape_theta[1], shape_mpo[0], shape_mpo[0] );
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       env2R(env2R_ptr,shape_theta[2], shape_theta[2], shape_mpo[1], shape_mpo[1]  );
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       mpo    (mpo_ptr, shape_mpo);
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       theta_in  (theta_in_, shape_theta);
    Eigen::TensorMap<Eigen::Tensor<T, 3>>             theta_out (theta_out_, shape_theta);

    double log2spin  = std::log2(shape_theta[0]);
    double log2chiL  = std::log2(shape_theta[1]);
    double log2chiR  = std::log2(shape_theta[2]);

    if (log2spin >= std::max(log2chiL, log2chiR)){
        if (log2chiL > log2chiR){
            Eigen::Tensor<Scalar,3> theta = theta_in.shuffle(Textra::array3{1,0,2});
            theta_out.device(omp->dev) =
                theta
                    .contract(env2L, Textra::idx({0}, {0}))
                    .contract(mpo  , Textra::idx({0,3}, {2,0}))
                    .contract(env2R, Textra::idx({0,3}, {0,2}))
                    .contract(mpo  , Textra::idx({2,1,4}, {2,0,1}))
                    .shuffle(Textra::array3{2,0,1});
        }

        else{
            Eigen::Tensor<Scalar,3> theta = theta_in.shuffle(Textra::array3{2,0,1});
            theta_out.device(omp->dev) =
                theta
                    .contract(env2R, Textra::idx({0}, {0}))
                    .contract(mpo  , Textra::idx({0,3}, {2,1}))
                    .contract(env2L, Textra::idx({0,3}, {0,2}))
                    .contract(mpo  , Textra::idx({2,4,1}, {2,0,1}))
                    .shuffle(Textra::array3{2,1,0});
        }

    }else{
        Eigen::Tensor<Scalar,3> theta = theta_in.shuffle(Textra::array3{1,0,2});
        theta_out.device(omp->dev) =
            theta
                .contract(env2L, Textra::idx({0}, {0}))
                .contract(mpo  , Textra::idx({0,3}, {2,0}))
                .contract(mpo  , Textra::idx({4,2}, {2,0}))
                .contract(env2R, Textra::idx({0,2,3}, {0,2,3}))
                .shuffle(Textra::array3{1, 0, 2});
    }

    Eigen::Tensor<Scalar,3> theta_add = theta_out + 100*theta_out;
    theta_out = theta_add;
    counter++;
    t_mul.toc();

}

//Explicit instantiations
template class DenseHamiltonianSqProduct<double>;
template class DenseHamiltonianSqProduct<std::complex<double>>;

