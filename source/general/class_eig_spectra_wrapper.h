//
// Created by david on 2017-11-03.
//

#ifndef DMRG_CLASS_EIG_WRAPPER_H
#define DMRG_CLASS_EIG_WRAPPER_H

#include "n_tensor_extra.h"
#include <SymEigsSolver.h>
#include <GenEigsSolver.h>

using namespace Textra;

enum class eig_Form{SYM, GEN};
enum class eig_Side {L,R}; //Left or right eigenvectors
using      eig_Mode = Spectra::SELECT_EIGENVALUE;

template<typename Scalar>
class class_eig {
private:
    double eigThreshold = 1e-12;
    int eigMaxIter = 10000;
    int ncv        = 20;

    using MatMapR = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>;
    using MatMapL = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;
    using MatrixType = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    template <long rank1, long rank2>
    using TensorPair = std::pair<Textra::Tensor<rank1,Scalar>,
                                 Textra::Tensor<rank2,Scalar>>;
public:

    void setThreshold(double new_eigThreshold){  eigThreshold   = new_eigThreshold;    }
    void setMaxIter(double new_eigMaxIter)    {  eigMaxIter     = new_eigMaxIter;      }
    void setncv(double new_ncv)               {  ncv            = new_ncv;             }

    Textra::Tensor<2,Scalar>  squareRoot(const Textra::Tensor<2,Scalar> &tensor);
    TensorPair<2,2> squareRoot_w_inverse(const Textra::Tensor<2,Scalar> &tensor);
    TensorPair<2,1> solve_gen(const Textra::Tensor<2,Scalar> &tensor);
    TensorPair<2,1> solve_sym(const Textra::Tensor<2,Scalar> &tensor);

    template<eig_Form form, eig_Mode mode, eig_Side side>
    auto solve_dominant(const Textra::Tensor<2, Scalar> &tensor, const int num, Textra::array2 size);

};


//Definitions

template<typename Scalar>
Textra::Tensor<2,Scalar> class_eig<Scalar>::squareRoot(const Textra::Tensor<2,Scalar> &tensor){
    MatMapR mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::SelfAdjointEigenSolver<MatrixType> eig_solver;
    eig_solver.compute(mat);
    return Textra::Matrix_to_Tensor<2,Scalar>(eig_solver.operatorSqrt().template cast<Scalar>(), tensor.dimensions());
};


template<typename Scalar>
typename class_eig<Scalar>:: template TensorPair<2,2> class_eig<Scalar>::squareRoot_w_inverse(const Textra::Tensor<2,Scalar> &tensor){
    Eigen::Map<const MatrixType> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::SelfAdjointEigenSolver<MatrixType> eig_solver;
    eig_solver.compute(mat);
    return std::make_pair(
            Textra::Matrix_to_Tensor<2,Scalar>(eig_solver.operatorSqrt().template cast<Scalar>()       , tensor.dimensions()),
            Textra::Matrix_to_Tensor<2,Scalar>(eig_solver.operatorInverseSqrt().template cast<Scalar>(), tensor.dimensions()));
};



template<typename Scalar>
typename class_eig<Scalar>:: template TensorPair<2,1>  class_eig<Scalar>::solve_gen(const Textra::Tensor<2,Scalar> &tensor) {
    Eigen::Map<const MatrixType> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::EigenSolver<MatrixType> eig_solver;
    eig_solver.compute(mat);
    return std::make_pair (
            Textra::Matrix_to_Tensor<2,std::complex<double>>(eig_solver.eigenvectors(), tensor.dimensions()),
            Textra::Matrix_to_Tensor<1,std::complex<double>>(eig_solver.eigenvalues() , {eig_solver.eigenvalues().size()} ));
}

template<typename Scalar>
typename class_eig<Scalar>:: template TensorPair<2,1>  class_eig<Scalar>::solve_sym(const Textra::Tensor<2,Scalar> &tensor) {
    Eigen::Map<const MatrixType> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::SelfAdjointEigenSolver<MatrixType> eig_solver;
    eig_solver.compute(mat);
    return std::make_pair(
            Textra::Matrix_to_Tensor<2,double>(eig_solver.eigenvectors(), tensor.dimensions()).template cast<Scalar>(),
            Textra::Matrix_to_Tensor<1,double>(eig_solver.eigenvalues() , {eig_solver.eigenvalues().size()}).template cast<Scalar>());
}


template<typename Scalar>
template<eig_Form form, eig_Mode mode, eig_Side side>
auto class_eig<Scalar>::solve_dominant(const Textra::Tensor<2, Scalar> &tensor, const int num, Textra::array2 size) {
    int ncv_max = 40;
    int ncv = std::min((int) std::sqrt(tensor.size()), ncv_max);
    ncv = std::max(ncv,num+2);
    if constexpr(side == eig_Side::R){
        MatMapR mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
        if constexpr(form == eig_Form::SYM){
            Spectra::DenseSymMatProd<double> op(mat);
            Spectra::SymEigsSolver<double, mode, Spectra::DenseSymMatProd<double>>eigs (&op, num, ncv);
            eigs.init();
            eigs.compute(eigMaxIter, eigThreshold, mode);
            if(eigs.info() != Spectra::SUCCESSFUL){std::cout << "SYM R Eigenvalue solver failed." << '\n'; exit(1);}
            return std::make_pair(
                    Textra::Matrix_to_Tensor<2,double>(eigs.eigenvectors(), size),
                    Textra::Matrix_to_Tensor<1,double>(eigs.eigenvalues(), {eigs.eigenvalues().size()}));
        }
        if constexpr(form == eig_Form::GEN){
            Spectra::DenseGenMatProd<double> op(mat);
            Spectra::GenEigsSolver<double, mode, Spectra::DenseGenMatProd<double>> eigs (&op, num, ncv);
            eigs.init();
            std::cout << mat << std::endl;
            eigs.compute(eigMaxIter, eigThreshold, mode);
            if(eigs.info() != Spectra::SUCCESSFUL){std::cout << "GEN R Eigenvalue solver failed." << '\n';exit(1);}
            return std::make_pair(
                    Textra::Matrix_to_Tensor<2,std::complex<double>>(eigs.eigenvectors(), size),
                    Textra::Matrix_to_Tensor<1,std::complex<double>>(eigs.eigenvalues(), {eigs.eigenvalues().size()}));
        }
    }
    if constexpr(side == eig_Side::L){
        MatMapL mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
        if constexpr(form == eig_Form::SYM){
            Spectra::DenseSymMatProd<double> op(mat);
            Spectra::SymEigsSolver<double, mode, Spectra::DenseSymMatProd<double>>eigs (&op, num, ncv);
            eigs.init();
            eigs.compute(eigMaxIter, eigThreshold, mode);
            if(eigs.info() != Spectra::SUCCESSFUL){std::cout << "SYM L Eigenvalue solver failed." << '\n';exit(1);}
            return std::make_pair(
                    Textra::Matrix_to_Tensor<2,double>(eigs.eigenvectors(), size),
                    Textra::Matrix_to_Tensor<1,double>(eigs.eigenvalues(), {eigs.eigenvalues().size()}));
        }
        if constexpr(form == eig_Form::GEN){
            Spectra::DenseGenMatProd<double> op(mat.transpose());
            Spectra::GenEigsSolver<double, mode, Spectra::DenseGenMatProd<double>>eigs (&op, num, ncv);
            eigs.init();
            eigs.compute(eigMaxIter, eigThreshold, mode);
            if(eigs.info() != Spectra::SUCCESSFUL){std::cout << "GEN R Eigenvalue solver failed." << '\n';exit(1);}
            return std::make_pair(
                    Textra::Matrix_to_Tensor<2,std::complex<double>>(eigs.eigenvectors(), size),
                    Textra::Matrix_to_Tensor<1,std::complex<double>>(eigs.eigenvalues(), {eigs.eigenvalues().size()}));
        }
    }
}


#endif //DMRG_CLASS_EIG_WRAPPER_H
