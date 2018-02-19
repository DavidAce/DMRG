//
// Created by david on 2017-11-14.
//

#ifndef DMRG_CLASS_EIG_ARPACK_WRAPPER_H
#define DMRG_CLASS_EIG_ARPACK_WRAPPER_H
#ifdef MKL_AVAILABLE
#define  EIGEN_USE_MKL_ALL
#endif

#include <map>
#include <complex>
#include "ardsnsym.h"
#include "ardscomp.h"
#include "ardgcomp.h"
#include "ardssym.h"
#include <Eigen/Core>


namespace arpack{
    enum class Form{SYMMETRIC, GENERAL, COMPLEX};  // Real Symmetric, Real General or Complex General
    enum class Side {L,R};           //Left or right eigenvectors
    enum class Ritz {LA,SA,LM,SM,LR,SR,LI,SI,BE}; //choice of eigenvalue. LA is largest algebraic, and so on.
}

class class_arpackpp_wrapper {
private:
    double eigThreshold = 1e-12;
    int    eigMaxIter   = 1000;
    int    ncv_max = 10;

    using  MapType = std::map<arpack::Ritz, std::string>;
    MapType RitzToString;


    template<typename T>
    using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;
    template<typename T>
    using  VectorType = Eigen::Matrix<T,Eigen::Dynamic, 1>;



    template <bool return_eigenvectors = true,typename Derived>
    auto retrieve_solution(Derived & solution, int nev) {
        using namespace std::complex_literals;
        if constexpr(return_eigenvectors){
            solution.FindEigenvectors();
            using T_vec = decltype(solution.Eigenvector(0, 0));
            using T_val = decltype(solution.Eigenvalue(0));
            int rows = solution.GetN();
            int cols = std::min(nev, solution.GetNev());
            MatrixType<T_vec> eigvecs(rows, cols);
            VectorType<T_val> eigvals(cols);
            for (int i = 0; i < cols; i++) {
                eigvals(i) = solution.Eigenvalue(0);
                double phase = std::arg(solution.Eigenvector(i, 0));
                for (int j = 0; j < rows; j++) {
                    eigvecs(j, i) = solution.Eigenvector(i, j);
                    if constexpr(std::is_same_v<T_vec, std::complex<double>>){
                        std::complex<double> inv_phase = -1.0i * phase;
                        eigvecs(j, i) *= std::exp(inv_phase);
                    }
                }
            }
            return std::make_pair(eigvecs, eigvals);
        }
        if constexpr (return_eigenvectors == false){
            solution.FindEigenvalues();
            using T_val = decltype(solution.Eigenvalue(0));
            int cols = solution.GetNev();
            VectorType<T_val> eigvals(cols);
            for (int i = 0; i < cols; i++) {
                eigvals(i) = solution.Eigenvalue(0);
            }
            return eigvals;
        }
    }

public:
    class_arpackpp_wrapper() {
        RitzToString = {
                {arpack::Ritz::LA, "LA"},
                {arpack::Ritz::SA, "SA"},
                {arpack::Ritz::LM, "LM"},
                {arpack::Ritz::SM, "SM"},
                {arpack::Ritz::LR, "LR"},
                {arpack::Ritz::SR, "SR"},
                {arpack::Ritz::LI, "LI"},
                {arpack::Ritz::SI, "SI"},
                {arpack::Ritz::BE, "BE"}
        };

    }



    template<arpack::Form form, arpack::Ritz ritz, arpack::Side side = arpack::Side::R, bool return_eigenvectors = true, typename Scalar>
    auto solve_dominant(Scalar *data, int rows, int cols, const int nev, Scalar *residp =  NULL) {
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>> inputdata(data, rows,cols );
        assert(rows == cols && "Input matrix is not square. Error in Arpack eigenvalue solver.");
        int n   = rows;
        int ncv = std::max(n/2, 4*nev);
        ncv = std::min(ncv, ncv_max);
        Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> internal_data;
        if constexpr(side == arpack::Side::R){
            internal_data = inputdata;
        }else{
            internal_data = inputdata.transpose();
        }


        if constexpr(form == arpack::Form::SYMMETRIC && std::is_same_v<Scalar, double>){
            ARdsSymMatrix<double> matrix(n, internal_data.data());
            ARluSymStdEig<double> eigs(nev, matrix, RitzToString.at(ritz), ncv, eigThreshold, eigMaxIter, residp);
            return retrieve_solution<return_eigenvectors>(eigs,nev);
        }
        if constexpr(form == arpack::Form::GENERAL && std::is_same_v<Scalar, double>){
            int nev_temp = nev == 1 ? 2 : nev;
            ARdsNonSymMatrix<double,double> matrix(n, internal_data.data());
            ARluNonSymStdEig<double> eigs(nev_temp, matrix, RitzToString.at(ritz), ncv, eigThreshold, eigMaxIter,residp);
            return retrieve_solution<return_eigenvectors>(eigs,nev);
        }

        if constexpr(form == arpack::Form::COMPLEX && std::is_same_v<Scalar, std::complex<double>>){
            int nev_temp = nev == 1 ? 2 : nev;
            ARdsNonSymMatrix<std::complex<double>,double> matrix(n, internal_data.data());
            ARluCompStdEig<double> eigs(nev_temp, matrix, RitzToString.at(ritz), ncv, eigThreshold, eigMaxIter,residp);
//            std::cout << "rows: " << rows << std::endl
//                      << "cols: " << cols << std::endl
//                      << "ncv:  " << ncv  << std::endl
//                      << "Ritz: " << RitzToString.at(ritz) <<std::endl;
            return retrieve_solution<return_eigenvectors>(eigs,nev);

        }
    }
};

#endif //DMRG_CLASS_EIG_ARPACK_WRAPPER_H
