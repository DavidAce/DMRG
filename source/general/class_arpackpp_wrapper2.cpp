//
// Created by david on 2018-05-06.
//

#include "class_arpackpp_wrapper2.h"
#include <algorithm>
#include <arpackpp/ardsnsym.h>
#include <arpackpp/ardscomp.h>
#include <arpackpp/ardgcomp.h>
#include <arpackpp/ardssym.h>
#include <assert.h>


template<typename Scalar, Form form>
class_arpackpp_wrapper2<Scalar, form>::class_arpackpp_wrapper2() {
    RitzToString = {
            {arpackpp::Ritz::LA, "LA"},
            {arpackpp::Ritz::SA, "SA"},
            {arpackpp::Ritz::LM, "LM"},
            {arpackpp::Ritz::SM, "SM"},
            {arpackpp::Ritz::LR, "LR"},
            {arpackpp::Ritz::SR, "SR"},
            {arpackpp::Ritz::LI, "LI"},
            {arpackpp::Ritz::SI, "SI"},
            {arpackpp::Ritz::BE, "BE"}
    };

}
template<typename Scalar, Form form>
const std::vector<Scalar> & class_arpackpp_wrapper2<Scalar, form>::get_eigvecs() const{
    if(eigvecs_found) {
        return eigvecs;
    }else{
        std::cerr << "Eigenvectors haven't been computed yet. Exiting " << std::endl;
        exit(1);
    }
}

template<typename Scalar, Form form>
const std::vector<Scalar> & class_arpackpp_wrapper2<Scalar, form>::get_eigvals() const {
    if(eigvals_found) {
        return eigvals;
    }else{
        std::cerr << "Eigenvalues haven't been computed yet. Exiting " << std::endl;
        exit(1);
    }
}


template<typename Scalar, Form form>
void class_arpackpp_wrapper2<Scalar, form>::dephase() {
    if constexpr (std::is_same_v<Scalar, std::complex<double>>) {
        using namespace std::complex_literals;
        for (int i = 0; i < nev_found; i++) {
            auto begin = eigvecs.begin() + i * rows;
            auto end = begin + rows;
            Scalar inv_phase = -1.0i * std::arg(eigvecs[i * rows]);
            std::transform(begin, end, begin,
                           std::bind1st(std::multiplies<Scalar>(), inv_phase));
        }
    }
}


template<typename Scalar, Form form>
void class_arpackpp_wrapper2<Scalar,form>::eig(Scalar *data,
                                          Ritz ritz,
                                          const int n,
                                          const int nev,
                                          const int ncv,
                                          bool compute_eigv){
    Scalar * residp = NULL;
    if constexpr (form == Form::REAL_SYMMETRIC && std::is_same_v<Scalar, double>){
        ARdsSymMatrix<Scalar> matrix(n, data);
        ARluSymStdEig<Scalar> eigs(nev, matrix, RitzToString.at(ritz), ncv, eigThreshold, eigMaxIter, residp);
        retrieve_solution(eigs,nev,compute_eigv);
    }

    if  constexpr(form == Form::REAL_GENERAL && std::is_same_v<Scalar, double>){
        int nev_temp = nev == 1 ? 2 : nev;
        ARdsNonSymMatrix<Scalar,Scalar> matrix(n, data);
        ARluNonSymStdEig<Scalar> eigs(nev_temp, matrix, RitzToString.at(ritz), ncv, eigThreshold, eigMaxIter,residp);
        retrieve_solution(eigs,nev,compute_eigv);
    }

    if constexpr (form == Form::COMPLEX_GENERAL && std::is_same_v<Scalar, std::complex<double>>){
        ARdsNonSymMatrix<std::complex<double>,double> matrix(n, data);
        ARluCompStdEig<double> eigs(nev, matrix, RitzToString.at(ritz), ncv, eigThreshold, eigMaxIter,residp);
        retrieve_solution(eigs,nev, compute_eigv);
    }


}

template<typename Scalar, Form form>
void class_arpackpp_wrapper2<Scalar,form>::eig(Scalar *data,
                                               std::string ritz,
                                               const int n,
                                               const int nev,
                                               const int ncv,
                                               bool compute_eigv)
{
    std::string ritz_str;
    Ritz ritz_enum;
    for (auto it = RitzToString.begin(); it != RitzToString.end(); ++it) {
        if (it->second == ritz){
            ritz_str = it->second;
            ritz_enum = it->first;
        }
    }
    assert(ritz_str == ritz);
    eig(data,ritz_enum, n, nev, ncv,compute_eigv);
};


template class class_arpackpp_wrapper2<double, Form::REAL_SYMMETRIC> ;
template class class_arpackpp_wrapper2<double, Form::REAL_GENERAL> ;
template class class_arpackpp_wrapper2<std::complex<double>, Form::COMPLEX_GENERAL> ;
