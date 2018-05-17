//
// Created by david on 2018-05-08.
//

#include <general/class_arpack_custom_products.h>
#include <general/nmspc_tensor_extra.h>
#include <Eigen/Core>


template<typename T,Form form>
void DenseMatrixProduct<T,form>::MultMv(T* x_in, T* x_out) {
    Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>> A_mat(A,n,n);
    if constexpr(form == Form::GENERAL) {
        if (side == Side::R) {
            Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>> x_vec_in  (x_in, n);
            Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>> x_vec_out (x_out, n);
            x_vec_out.noalias() = A_mat * x_vec_in;
        }
        if (side == Side::L) {
            Eigen::Map<Eigen::Matrix<T,1,Eigen::Dynamic>> x_vec_in  (x_in, n);
            Eigen::Map<Eigen::Matrix<T,1,Eigen::Dynamic>> x_vec_out (x_out, n);
            x_vec_out.noalias() = x_vec_in * A_mat;
        }
    }
    if constexpr (form == Form::SYMMETRIC){
        Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>> x_vec_in  (x_in, n);
        Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>> x_vec_out (x_out, n);
        x_vec_out.noalias() = A_mat.template selfadjointView<Eigen::Lower>() * x_vec_in;
    }
    counter++;
}

template class DenseMatrixProduct<std::complex<double>, Form::GENERAL>;
template class DenseMatrixProduct<std::complex<double>, Form::SYMMETRIC>;
template class DenseMatrixProduct<double, Form::GENERAL>;
template class DenseMatrixProduct<double, Form::SYMMETRIC>;






template<class T>
void DenseHamiltonianProduct<T>::MultMv(T* theta_in_, T* theta_out_) {
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       Lblock_map(Lblock,shape_theta4[1], shape_theta4[1], shape_mpo4[0] );
    Eigen::TensorMap<Eigen::Tensor<const T, 3>>       Rblock_map(Rblock,shape_theta4[3], shape_theta4[3], shape_mpo4[1] );
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       HA_map    (HA, shape_mpo4);
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       HB_map    (HB, shape_mpo4);
    Eigen::TensorMap<Eigen::Tensor<const T, 4>>       theta_in (theta_in_, shape_theta4);
    Eigen::TensorMap<Eigen::Tensor<T, 4>>             theta_out(theta_out_, shape_theta4);



    //Best yet! The sparcity of the effective hamiltonian (Lblock HA HB Rblock) is about 58% nonzeros.
    //I have shown this to be the fastest contraction ordering
    theta_out = Lblock_map
            .contract(theta_in,    Textra::idx({0},{1}))
            .contract(HA_map ,     Textra::idx({1,2},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(HB_map ,     Textra::idx({3,1},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(Rblock_map,  Textra::idx({1,3},{0,2}))
            .shuffle(Textra::array4{1,0,2,3});
    counter++;
}

template class DenseHamiltonianProduct<std::complex<double>>;
template class DenseHamiltonianProduct<double>;

