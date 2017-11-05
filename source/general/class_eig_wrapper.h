//
// Created by david on 2017-11-03.
//

#ifndef DMRG_CLASS_EIG_WRAPPER_H
#define DMRG_CLASS_EIG_WRAPPER_H

#include <Util/SelectionRule.h>
#include "n_tensor_extra.h"
#include <SymEigsSolver.h>
#include <GenEigsSolver.h>

class class_eig {
private:
    double eigThreshold = 1e-12;
    int eigMaxIter = 10000;
    template<typename Scalar>
    using MatVecPair = std::pair<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> ,
                                 Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>;

    template<long rank1, long rank2, typename Scalar>
    using TensorPair = std::pair<Textra::Tensor<rank1,Scalar>, Textra::Tensor<rank2,Scalar>>;

public:
    using real = double;
    using cplx = std::complex<double>;
    using selectionRule = Spectra::SELECT_EIGENVALUE;
    enum inputForm{SYM, GEN};
    enum direction {L,R}; //Left or right eigenvectors

    TensorPair <2,1,cplx> solve_gen(const Textra::Tensor2d & tensor);
    TensorPair <2,1,real> solve_sym(const Textra::Tensor2d & tensor);

    template<inputForm form, selectionRule rule, direction dir>
    auto SpectraSelector(Eigen::MatrixXd mat, int num){
        if constexpr (dir == direction::L){mat.transposeInPlace();}
        if constexpr(form = inputForm::GEN){
            int ncv =  3;
            Spectra::DenseGenMatProd<double> op(mat);
            return Spectra::GenEigsSolver<double, rule, Spectra::DenseGenMatProd<double>> (&op, num, ncv);
        }
        if constexpr(form = inputForm::SYM){
            int ncv =  2;
            Spectra::DenseSymMatProd<double> op(mat);
            return Spectra::SymEigsSolver<double, rule, Spectra::DenseSymMatProd<double>> (&op, num, ncv);
        }
    };

    template<typename Scalar, inputForm form, selectionRule rule, direction dir>
    MatVecPair<Scalar> solve_to_matrix(const Textra::Tensor2d &tensor, const int num) {
        Eigen::Map<const Eigen::MatrixXd> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
        auto eigs = SpectraSelector<form, rule,dir>(mat,num);
        eigs.init();
        eigs.compute(eigMaxIter, eigThreshold, rule);
        if(eigs.info() != Spectra::SUCCESSFUL){
            std::cout << "Eigenvalue solver failed." << '\n';
        }
        std::cout << "Eigenvalue: " << eigs.eigenvalues()<< std::endl;
        return std::make_pair(eigs.eigenvectors(),
                              eigs.eigenvalues());

    }

    template<typename Scalar, inputForm form, selectionRule rule, direction dir>
    TensorPair<2,1,Scalar> solve_to_tensor(const Textra::Tensor2d &tensor, const int num) {
        auto[eigvec,eigval] = solve_to_matrix<Scalar, form, rule, dir>(tensor, num);
        return std::make_pair(
                Textra::Matrix_to_Tensor<2,Scalar>(eigvec, {eigvec.rows(),eigvec.cols()}),
                Textra::Matrix_to_Tensor<1,Scalar>(eigval, {eigval.size()})
        );
    }

    TensorPair<2,2,real> squareRoot_sym(const Textra::Tensor2d & tensor){
        Eigen::Map<const Eigen::MatrixXd> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver;
        eig_solver.compute(mat);
        return std::make_pair(
                Textra::Matrix_to_Tensor<2,real>(eig_solver.operatorSqrt()       , tensor.dimensions()),
                Textra::Matrix_to_Tensor<2,real>(eig_solver.operatorInverseSqrt(), tensor.dimensions())
        );
    };
};


#endif //DMRG_CLASS_EIG_WRAPPER_H
