////
//// Created by david on 2017-11-03.
////
//
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "class_eig_wrapper.h"

class_eig::TensorPair <2,1,class_eig::cplx>  class_eig::solve_gen(const Textra::Tensor2d &tensor) {
    Eigen::Map<const Eigen::MatrixXd> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::EigenSolver<Eigen::MatrixXd> eig_solver;
    eig_solver.compute(mat);
    return std::make_pair(
            Textra::Matrix_to_Tensor<2,class_eig::cplx>(eig_solver.eigenvectors(), tensor.dimensions()),
            Textra::Matrix_to_Tensor<2,class_eig::cplx>(eig_solver.eigenvalues() , {eig_solver.eigenvalues().size()} )
    );
}

class_eig::TensorPair <2,1,class_eig::real>  class_eig::solve_sym(const Textra::Tensor2d &tensor) {
    Eigen::Map<const Eigen::MatrixXd> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver;
    eig_solver.compute(mat);
    return std::make_pair(
            Textra::Matrix_to_Tensor<2,class_eig::real>(eig_solver.eigenvectors(), tensor.dimensions()),
            Textra::Matrix_to_Tensor<2,class_eig::real>(eig_solver.eigenvalues() , {eig_solver.eigenvalues().size()} )
    );
}

