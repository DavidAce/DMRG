//
// Created by david on 2018-05-04.
//
//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACK
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <mps_routines/class_custom_contraction.h>
using Scalar = std::complex<double>;

template<class T>
void class_custom_contraction<T>::MultMv(T* theta_in_, T* theta_out_) {
    Eigen::TensorMap<Textra::Tensor<const T, 4>> theta_in (theta_in_, shape4);
    Eigen::TensorMap<Textra::Tensor<T, 1>>       theta_out(theta_out_, shape1);
    //Best yet! The sparcity of the effective hamiltonian (Lblock HA Rblock) is about 58% nonzeros.

    theta_out = Lblock
            .contract(theta_in, Textra::idx({0},{1}))
            .contract(HA ,     Textra::idx({1,2},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(HB ,     Textra::idx({3,1},{0,2}))//  idx({1,2,3},{0,4,5}))
            .contract(Rblock,  Textra::idx({1,3},{0,2}))
            .shuffle(Textra::array4{1,0,2,3})
            .reshape(shape1);
    counter++;

//    Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>> mymat(mytens.data(), shape1[0],shape1[0]);
//    Eigen::TensorMap<Textra::Tensor<const T, 1>> theta_in(theta_in_, shape1[0]);
//    Eigen::TensorMap<Textra::Tensor<T, 1>>       theta_out(theta_out_, shape1[0]);
//    theta_out = mytens.contract(theta_in, Textra::idx({1},{0}));

//
//    Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>> theta_in(theta_in_, shape1[0]);
//    Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>> theta_out(theta_out_, shape1[0]);
//    theta_out = mymat * theta_in;
//


//    Eigen::TensorMap<Textra::Tensor<const T, 4>> theta_in (theta_in_, shape4);
//    Eigen::TensorMap<Textra::Tensor<T, 1>>       theta_out(theta_out_, shape1);
//    //Best yet! The sparcity of the effective hamiltonian (Lblock HA Rblock) is about 58% nonzeros.
//
//    theta_out = theta_in
//            .contract(Lblock, Textra::idx({1},{0}))
//            .contract(HA ,     Textra::idx({4,0},{0,2}))//  idx({1,2,3},{0,4,5}))
//            .contract(HB ,     Textra::idx({3,0},{0,2}))//  idx({1,2,3},{0,4,5}))
//            .contract(Rblock,  Textra::idx({0,3},{0,2}))
//            .shuffle(Textra::array4{1,0,2,3})
//            .reshape(shape1);
//
//


//
//    Eigen::TensorMap<Textra::Tensor<const T, 4>> theta_in (theta_in_, shape4);
//    Eigen::TensorMap<Textra::Tensor<T, 1>>       theta_out(theta_out_, shape1);
//    //Best yet! The sparcity of the effective hamiltonian (Lblock HA Rblock) is about 58% nonzeros.
//
//    theta_out = Lblock
//            .contract(theta_in,    Textra::idx({0},{1}))
//            .contract(Rblock ,     Textra::idx({4},{0}))//  idx({1,2,3},{0,4,5}))
//            .contract(HA ,         Textra::idx({1,2},{0,2}))//  idx({1,2,3},{0,4,5}))
//            .contract(HB,          Textra::idx({4,3,1},{0,1,2}))
//            .shuffle(Textra::array4{2,0,3,1})
//            .reshape(shape1);
//



}

//template class class_custom_contraction<double>;
template class class_custom_contraction<Scalar>;
