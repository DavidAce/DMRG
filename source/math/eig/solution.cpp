#include "solution.h"

void eig::solution::build_eigvecs_cplx() {
    bool build_eigvecsR_cplx = eigvecsR_cplx.empty() and (not eigvecsR_real.empty() or not eigvecsR_imag.empty());
    bool build_eigvecsL_cplx = eigvecsL_cplx.empty() and (not eigvecsL_real.empty() or not eigvecsL_imag.empty());

    if(build_eigvecsR_cplx) {
        eigvecsR_cplx.resize(std::max(eigvecsR_real.size(), eigvecsR_imag.size()));
        for(size_t i = 0; i < eigvecsR_cplx.size(); i++) {
            if(not eigvecsR_real.empty() and not eigvecsR_imag.empty() and i < eigvecsR_real.size() and i < eigvecsR_imag.size())
                eigvecsR_cplx[i] = std::complex<double>(eigvecsR_real[i], eigvecsR_imag[i]);
            else if(not eigvecsR_real.empty() and i < eigvecsR_real.size())
                eigvecsR_cplx[i] = std::complex<double>(eigvecsR_real[i], 0.0);
            else if(not eigvecsR_imag.empty() and i < eigvecsR_imag.size())
                eigvecsR_cplx[i] = std::complex<double>(0.0, eigvecsR_imag[i]);
        }
        eigvecsR_real.clear();
        eigvecsR_imag.clear();
    }
    if(build_eigvecsL_cplx) {
        eigvecsL_cplx.resize(std::max(eigvecsL_real.size(), eigvecsL_imag.size()));
        for(size_t i = 0; i < eigvecsL_cplx.size(); i++) {
            if(not eigvecsL_real.empty() and not eigvecsL_imag.empty() and i < eigvecsL_real.size() and i < eigvecsL_imag.size())
                eigvecsL_cplx[i] = std::complex<double>(eigvecsL_real[i], eigvecsL_imag[i]);
            else if(not eigvecsL_real.empty() and i < eigvecsL_real.size())
                eigvecsL_cplx[i] = std::complex<double>(eigvecsL_real[i], 0.0);
            else if(not eigvecsL_imag.empty() and i < eigvecsL_imag.size())
                eigvecsL_cplx[i] = std::complex<double>(0.0, eigvecsL_imag[i]);
        }
        eigvecsL_real.clear();
        eigvecsL_imag.clear();
    }
}

void eig::solution::build_eigvecs_real() {
    bool build_eigvecsR_real = eigvecsR_real.empty() and not eigvecsR_cplx.empty();
    bool build_eigvecsL_real = eigvecsL_real.empty() and not eigvecsL_cplx.empty();

    if(build_eigvecsR_real) {
        eigvecsR_real.resize(eigvecsR_cplx.size());
        for(size_t i = 0; i < eigvecsR_real.size(); i++) {
            if(std::imag(eigvecsR_cplx[i]) > 1e-12) throw std::runtime_error("Error building real eigvecR: Nonzero imaginary part");
            eigvecsR_real[i] = std::real(eigvecsR_cplx[i]);
        }
        eigvecsR_cplx.clear();
    }
    if(build_eigvecsL_real) {
        eigvecsL_real.resize(eigvecsL_cplx.size());
        for(size_t i = 0; i < eigvecsL_real.size(); i++) {
            if(std::imag(eigvecsL_cplx[i]) > 1e-12) throw std::runtime_error("Error building real eigvecL: Nonzero imaginary part");
            eigvecsL_real[i] = std::real(eigvecsL_cplx[i]);
        }
        eigvecsL_cplx.clear();
    }
}

void eig::solution::build_eigvals_cplx() {
    bool build_cplx = eigvals_cplx.empty() and not eigvals_imag.empty() and eigvals_real.size() == eigvals_imag.size();
    if(build_cplx) {
        eigvals_cplx.resize(eigvals_real.size());
        for(size_t i = 0; i < eigvals_real.size(); i++) eigvals_cplx[i] = std::complex<double>(eigvals_real[i], eigvals_imag[i]);
        eigvals_real.clear();
        eigvals_imag.clear();
    }
}

void eig::solution::build_eigvals_real() {
    bool build_real = (eigvals_real.empty() or eigvals_imag.empty()) and not eigvals_cplx.empty();
    if(build_real) {
        eigvals_real.resize(eigvals_cplx.size());
        eigvals_imag.resize(eigvals_cplx.size());
        for(size_t i = 0; i < eigvals_cplx.size(); i++) {
            eigvals_real[i] = std::real(eigvals_cplx[i]);
            eigvals_imag[i] = std::imag(eigvals_cplx[i]);
        }
        eigvals_cplx.clear();
    }
}