//
// Created by david on 2020-06-04.
//

#include "solution.h"



void eig::solution::build_eigvecs_cplx(){
    bool build_eigvecsR_cplx =  eigvecsR_cplx.empty() and not eigvecsR_imag.empty() and eigvecsR_real.size() == eigvecsR_imag.size();
    bool build_eigvecsL_cplx =  eigvecsL_cplx.empty() and not eigvecsL_imag.empty() and eigvecsL_real.size() == eigvecsL_imag.size();

    if(build_eigvecsR_cplx){
        eigvecsR_cplx.resize(eigvecsR_real.size());
        for(size_t i = 0; i < eigvecsR_real.size(); i++)
            eigvecsR_cplx[i] = std::complex<double>(eigvecsR_real[i], eigvecsR_imag[i]);
        eigvecsR_real.clear();
        eigvecsR_imag.clear();
    }
    if(build_eigvecsL_cplx){
        eigvecsL_cplx.resize(eigvecsL_real.size());
        for(size_t i = 0; i < eigvecsL_real.size(); i++)
            eigvecsL_cplx[i] = std::complex<double>(eigvecsL_real[i], eigvecsL_imag[i]);
        eigvecsL_real.clear();
        eigvecsL_imag.clear();
    }

}

void eig::solution::build_eigvals_cplx(){
    bool build_cplx =  eigvals_cplx.empty() and not eigvals_imag.empty() and eigvals_real.size() == eigvals_imag.size();
    if(build_cplx){
        eigvals_cplx.resize(eigvals_real.size());
        for(size_t i = 0; i < eigvals_real.size(); i++)
            eigvals_cplx[i] = std::complex<double>(eigvals_real[i], eigvals_imag[i]);
        eigvals_real.clear();
        eigvals_imag.clear();
    }
}