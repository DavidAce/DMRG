//
// Created by david on 2017-10-04.
//

#ifndef DMRG_CLASS_SVD_H
#define DMRG_CLASS_SVD_H

#include "nmspc_tensor_extra.h"
#include <Eigen/SVD>

template<typename Scalar>
class class_SVD{
private:
    double SVDThreshold         = 1e-12;
    double truncation_error     = 0;
    int chi                     = 0;
    using MatrixType  = Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    Eigen::BDCSVD<MatrixType> SVD;
public:

    class_SVD(){
        setThreshold(SVDThreshold);
    }

    double get_truncation_error();
    void setThreshold(double newThreshold);
    Eigen::Tensor<Scalar, 2> pseudo_inverse(const Eigen::Tensor<Scalar,2> &tensor);

    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,2> &tensor);

    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,2> &tensor, const long chi_max);

    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,2> &tensor, long d, long chiL, long chi_max, long chiR);

    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,2> &tensor);

    std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,4> &tensor, long chi_max);

    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,4> &tensor);
};

#endif //DMRG_CLASS_SVD_H
