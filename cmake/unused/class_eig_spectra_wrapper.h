//
// Created by david on 2017-11-03.
//


#ifndef DMRG_CLASS_EIG_WRAPPER_H
#define DMRG_CLASS_EIG_WRAPPER_H

#ifdef MKL_AVAILABLE
#define  EIGEN_USE_MKL_ALL
#endif

#include "general/nmspc_tensor_extra.h"
#include <SymEigsSolver.h>
#include <GenEigsSolver.h>
#include <experimental/filesystem>

//using namespace Textra;

namespace Spectra{
enum class Form{SYMMETRIC, GENERAL};
enum class Side {L,R}; //Left or right eigenvectors
using      Ritz = Spectra::SELECT_EIGENVALUE;
}

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
    using TensorPair = std::pair<Textra::Tensor<Scalar,rank1>,
                                 Textra::Tensor<Scalar,rank2>>;
public:

    void setThreshold(double new_eigThreshold){  eigThreshold   = new_eigThreshold;    }
    void setMaxIter(double new_eigMaxIter)    {  eigMaxIter     = new_eigMaxIter;      }
    void setncv(double new_ncv)               {  ncv            = new_ncv;             }

    Textra::Tensor<Scalar,2>  squareRoot(const Textra::Tensor<Scalar,2> &tensor);
    TensorPair<2,2> squareRoot_w_inverse(const Textra::Tensor<Scalar,2> &tensor);
    TensorPair<2,1> solve_gen(const Textra::Tensor<Scalar,2> &tensor);
    TensorPair<2,1> solve_sym(const Textra::Tensor<Scalar,2> &tensor);

    template<Spectra::Form form, Spectra::Ritz mode, Spectra::Side side = Spectra::Side::R>
    auto solve_dominant(const Textra::Tensor<Scalar,2> &tensor, const int num, Textra::array2 size);

};


//Definitions

template<typename Scalar>
Textra::Tensor<Scalar,2> class_eig<Scalar>::squareRoot(const Textra::Tensor<Scalar,2> &tensor){
    MatMapR mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::SelfAdjointEigenSolver<MatrixType> eig_solver;
    eig_solver.compute(mat);
    return Textra::Matrix_to_Tensor<2,Scalar>(eig_solver.operatorSqrt().template cast<Scalar>(), tensor.dimensions());
};


template<typename Scalar>
typename class_eig<Scalar>:: template TensorPair<2,2> class_eig<Scalar>::squareRoot_w_inverse(const Textra::Tensor<Scalar,2> &tensor){
    Eigen::Map<const MatrixType> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::SelfAdjointEigenSolver<MatrixType> eig_solver;
    eig_solver.compute(mat);
    return std::make_pair(
            Textra::Matrix_to_Tensor<Scalar,2>(eig_solver.operatorSqrt().template cast<Scalar>()       , tensor.dimensions()),
            Textra::Matrix_to_Tensor<Scalar,2>(eig_solver.operatorInverseSqrt().template cast<Scalar>(), tensor.dimensions()));
};



template<typename Scalar>
typename class_eig<Scalar>:: template TensorPair<2,1>  class_eig<Scalar>::solve_gen(const Textra::Tensor<Scalar,2> &tensor) {
    Eigen::Map<const MatrixType> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::EigenSolver<MatrixType> eig_solver;
    eig_solver.compute(mat);
    return std::make_pair (
            Textra::Matrix_to_Tensor<std::complex<double>,2>(eig_solver.eigenvectors(), tensor.dimensions()),
            Textra::Matrix_to_Tensor<std::complex<double>,1>(eig_solver.eigenvalues() , {eig_solver.eigenvalues().size()} ));
}

template<typename Scalar>
typename class_eig<Scalar>:: template TensorPair<2,1>  class_eig<Scalar>::solve_sym(const Textra::Tensor<Scalar,2> &tensor) {
    Eigen::Map<const MatrixType> mat(tensor.data(), tensor.dimension(0), tensor.dimension(1));
    Eigen::SelfAdjointEigenSolver<MatrixType> eig_solver;
    eig_solver.compute(mat);
    return std::make_pair(
            Textra::Matrix_to_Tensor<double,2>(eig_solver.eigenvectors(), tensor.dimensions()).template cast<Scalar>(),
            Textra::Matrix_to_Tensor<double,1>(eig_solver.eigenvalues() , {eig_solver.eigenvalues().size()}).template cast<Scalar>());
}


template<typename Scalar>
template<Spectra::Form form, Spectra::Ritz mode, Spectra::Side side>
auto class_eig<Scalar>::solve_dominant(const Textra::Tensor<Scalar,2> &tensor, const int num, Textra::array2 outsize) {
    int ncv_max = 80;
    int nev = 1;
    int n = (int) tensor.dimension(0);
    int ncv = std::max(n/2, 4*nev);
    ncv = std::min(ncv, ncv_max);

//
//    int ncv_max = 40;
//    int ncv = std::min((int) std::sqrt(tensor.size()), ncv_max);
//    ncv = std::max(ncv,num+2);
    if constexpr(side == Spectra::Side::R){
        MatMapR mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
        if constexpr(form == Spectra::Form::SYMMETRIC){
            Spectra::DenseSymMatProd<double> op(mat);
            Spectra::SymEigsSolver<double, mode, Spectra::DenseSymMatProd<double>>eigs (&op, num, ncv);
            eigs.init();
            eigs.compute(eigMaxIter, eigThreshold, mode);
            if(eigs.info() != Spectra::SUCCESSFUL){std::cout << "SYMMETRIC R Eigenvalue solver failed." << '\n'; exit(1);}
            return std::make_pair(
                    Textra::Matrix_to_Tensor<double,2>(eigs.eigenvectors(), outsize),
                    Textra::Matrix_to_Tensor<double,1>(eigs.eigenvalues(), {eigs.eigenvalues().size()}));
        }
        if constexpr(form == Spectra::Form::GENERAL){
            Spectra::DenseGenMatProd<double> op(mat);
            Spectra::GenEigsSolver<double, mode, Spectra::DenseGenMatProd<double>> eigs (&op, num, ncv);
            eigs.init();
            eigs.compute(eigMaxIter, eigThreshold, mode);
            if(eigs.info() != Spectra::SUCCESSFUL){std::cout << "GENERAL R Eigenvalue solver failed." << '\n';exit(1);}
            return std::make_pair(
                    Textra::Matrix_to_Tensor<std::complex<double>,2>(eigs.eigenvectors(), outsize),
                    Textra::Matrix_to_Tensor<std::complex<double>,1>(eigs.eigenvalues(), {eigs.eigenvalues().size()}));
        }
    }
    if constexpr(side == Spectra::Side::L){
        MatMapL mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
        if constexpr(form == Spectra::Form::SYMMETRIC){
            Spectra::DenseSymMatProd<double> op(mat);
            Spectra::SymEigsSolver<double, mode, Spectra::DenseSymMatProd<double>>eigs (&op, num, ncv);
            eigs.init();
            eigs.compute(eigMaxIter, eigThreshold, mode);
            if(eigs.info() != Spectra::SUCCESSFUL){std::cout << "SYMMETRIC L Eigenvalue solver failed." << '\n';exit(1);}
            return std::make_pair(
                    Textra::Matrix_to_Tensor<double,2>(eigs.eigenvectors(), outsize),
                    Textra::Matrix_to_Tensor<double,1>(eigs.eigenvalues(), {eigs.eigenvalues().size()}));
        }
        if constexpr(form == Spectra::Form::GENERAL){
            Spectra::DenseGenMatProd<double> op(mat.transpose());
            Spectra::GenEigsSolver<double, mode, Spectra::DenseGenMatProd<double>>eigs (&op, num, ncv);
            eigs.init();
            eigs.compute(eigMaxIter, eigThreshold, mode);
            if(eigs.info() != Spectra::SUCCESSFUL){std::cout << "GENERAL R Eigenvalue solver failed." << '\n';exit(1);}
            return std::make_pair(
                    Textra::Matrix_to_Tensor<std::complex<double>,2>(eigs.eigenvectors(), outsize),
                    Textra::Matrix_to_Tensor<std::complex<double>,1>(eigs.eigenvalues(), {eigs.eigenvalues().size()}));
        }
    }
}


#endif //DMRG_CLASS_EIG_WRAPPER_H
