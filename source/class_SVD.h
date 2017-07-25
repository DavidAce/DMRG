//
// Created by david on 7/18/17.
//

#ifndef FINITE_DMRG_EIGEN_CLASS_SVD_H
#define FINITE_DMRG_EIGEN_CLASS_SVD_H

#include <n_tensor_extra.h>
using namespace Textra;

class class_SVD {
private:
    double renorm;
public:
    double SVDThreshold;
    long chi;

    class_SVD(){};
    std::tuple<Tensor3, Tensor1, Tensor3> decompose(const Eigen::TensorRef<Tensor2> theta2, const double SVDThreshold , const long d, const long chi, const long chia, const long chib);
};


#endif //FINITE_DMRG_EIGEN_CLASS_SVD_H
