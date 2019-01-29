//
// Created by david on 2019-01-17.
//

#ifndef CLASS_MPO_H
#define CLASS_MPO_H

#include <general/nmspc_tensor_extra.h>

class class_mpo{

public:
    using Scalar = std::complex<double>;

    static std::tuple<
            Eigen::Tensor<class_mpo::Scalar,4>,
            Eigen::Tensor<class_mpo::Scalar,3>,
            Eigen::Tensor<class_mpo::Scalar,3>>
    pauli_mpo(const Eigen::MatrixXcd paulimatrix);


    static std::tuple<
            Eigen::Tensor<class_mpo::Scalar,4>,
            Eigen::Tensor<class_mpo::Scalar,3>,
            Eigen::Tensor<class_mpo::Scalar,3>>
    parity_selector_mpo(const Eigen::MatrixXcd paulimatrix, const int sector = 1);




};


#endif //CLASS_MPO_H
