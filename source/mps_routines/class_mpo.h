//
// Created by david on 2019-01-17.
//

#ifndef CLASS_MPO_H
#define CLASS_MPO_H

#include <general/nmspc_tensor_extra.h>

class class_mpo{
public:
    using Scalar = std::complex<double>;
    static Eigen::Tensor<Scalar,4> parity(const Eigen::Matrix2cd & paulimatrix);
    static Eigen::Tensor<Scalar,4> parity(const Eigen::Matrix3cd & paulimatrix);
};


#endif //CLASS_MPO_H
