//
// Created by david on 2017-11-14.
//

#ifndef DMRG_CLASS_EIG_ARPACK_WRAPPER_H
#define DMRG_CLASS_EIG_ARPACK_WRAPPER_H
#include <map>
#include <complex>
#include "n_tensor_extra.h"
#include "ardsnsym.h"
#include "ardscomp.h"
#include "ardgcomp.h"
#include "ardssym.h"

//using namespace Textra;
//using namespace std::complex_literals;

namespace arpack{
    enum class Form{SYMMETRIC, GENERAL, COMPLEX};  // Real Symmetric, Real General or Complex General
    enum class Side {L,R};           //Left or right eigenvectors
    enum class Ritz {LA,SA,LM,SM,LR,SR,LI,SI,BE}; //choice of eigenvalue. LA is largest algebraic, and so on.
}


class class_eig_arpack {
private:
    double eigThreshold = 1e-14;
    int    eigMaxIter   = 10000;
    using  MapType = std::map<arpack::Ritz, std::string>;
    MapType RitzToString;

    template <bool return_eigenvectors = true,typename Derived>
    auto retrieve_solution(Derived & solution, int nev) {
        using namespace std::complex_literals;
        if constexpr(return_eigenvectors){
            solution.FindEigenvectors();
            using T_vec = decltype(solution.Eigenvector(0, 0));
            using T_val = decltype(solution.Eigenvalue(0));
            int rows = solution.GetN();
            int cols = std::min(nev, solution.GetNev());
            Textra::Tensor<T_vec, 2> eigvecs(rows, cols);
            Textra::Tensor<T_val, 1> eigvals(cols);
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
            Textra::Tensor<T_val, 1> eigvals(cols);
            for (int i = 0; i < cols; i++) {
                eigvals(i) = solution.Eigenvalue(0);
            }
            return eigvals;
        }
    }

public:
    class_eig_arpack() {
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
    auto solve_dominant(const Textra::Tensor<Scalar,2> &tensor, const int nev) {
        assert(tensor.dimension(0) == tensor.dimension(1) &&
               "Input matrix is not square. Error in Arpack eigenvalue solver.");
        int ncv_max = 80;
        int n = (int)tensor.dimension(0);
        int ncv = std::max(n/2, 4*nev);
        ncv = std::min(ncv, ncv_max);
        Textra::Tensor<Scalar,2> internal_tensor;
        if constexpr(side == arpack::Side::R){
            internal_tensor = tensor;
        }else{
            internal_tensor = tensor.shuffle(Textra::array2{1,0});
        }


        if constexpr(form == arpack::Form::SYMMETRIC && std::is_same_v<Scalar, double>){
            ARdsSymMatrix<double> matrix(n, internal_tensor.data());
            ARluSymStdEig<double> eigs(nev, matrix, RitzToString.at(ritz), n, eigThreshold, eigMaxIter);
            return retrieve_solution<return_eigenvectors>(eigs,nev);
        }
        if constexpr(form == arpack::Form::GENERAL && std::is_same_v<Scalar, double>){
            int nev_temp = nev == 1 ? 2 : nev;
            ARdsNonSymMatrix<double,double> matrix(n, internal_tensor.data());
            ARluNonSymStdEig<double> eigs(nev_temp, matrix, RitzToString.at(ritz), n, eigThreshold, eigMaxIter);
            return retrieve_solution<return_eigenvectors>(eigs,nev);
        }

        if constexpr(form == arpack::Form::COMPLEX && std::is_same_v<Scalar, std::complex<double>>){
            ARdsNonSymMatrix<std::complex<double>,double> matrix(n, internal_tensor.data());
            ARluCompStdEig<double> eigs(nev, matrix, RitzToString.at(ritz), n, eigThreshold, eigMaxIter);
            return retrieve_solution<return_eigenvectors>(eigs,nev);
        }
    }
};

#endif //DMRG_CLASS_EIG_ARPACK_WRAPPER_H
