//
// Created by david on 2017-10-04.
//
#include <Eigen/SVD>
#include "class_svd_wrapper.h"

double class_SVD::get_truncation_error(){
    return truncation_error;
}

void class_SVD::setThreshold(double newThreshold) {
    SVDThreshold = newThreshold;
}

class_SVD::TensorTuple<2,1,2> class_SVD::decompose(const Textra::Tensor2d &tensor) {
    Eigen::BDCSVD<Eigen::MatrixXd> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(Textra::Tensor2d_to_MatrixXd(tensor), Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    return std::make_tuple
            (
                    Textra::MatrixXd_to_Tensor2d(SVD.matrixU()),
                    Textra::MatrixXd_to_Tensor2d(SVD.singularValues().normalized()),
                    Textra::MatrixXd_to_Tensor2d(SVD.matrixV().transpose())
            );
}

class_SVD::TensorTuple<2,1,2> class_SVD::decompose(const Textra::Tensor2d &tensor, const long chi_max) {
    Eigen::BDCSVD<Eigen::MatrixXd> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(Textra::Tensor2d_to_MatrixXd(tensor), Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    long chi = std::min(SVD.rank(), chi_max);
    return std::make_tuple
            (
                Textra::MatrixXd_to_Tensor2d(SVD.matrixU().leftCols(chi)),
                Textra::MatrixXd_to_Tensor2d(SVD.singularValues().head(chi).normalized()),
                Textra::MatrixXd_to_Tensor2d(SVD.matrixV().leftCols(chi).transpose())
            );
}

class_SVD::TensorTuple<3,1,3> class_SVD::schmidt(const Textra::Tensor2d &tensor, long d, long chiL, long chi_max, long chiR) {
    Eigen::BDCSVD<Eigen::MatrixXd> SVD;
    SVD.setThreshold(SVDThreshold);
    SVD.compute(Textra::Tensor2d_to_MatrixXd(tensor), Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi            = std::min(SVD.rank(),chi_max);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).sum();
    return std::make_tuple
            (
                Textra::Matrix_to_Tensor<3, double>(SVD.matrixU().leftCols(chi), {d, chiL, chi}),
                Textra::Matrix_to_Tensor<1, double>(SVD.singularValues().head(chi).normalized(), { chi }),
                Textra::Matrix_to_Tensor<3, double>(SVD.matrixV().leftCols(chi), { d, chiR, chi }).shuffle(Textra::array3{ 0, 2, 1 })
            );
}

