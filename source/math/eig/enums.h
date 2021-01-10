#pragma once
#include <complex>
#include <tools/common/log.h>

namespace eig {
    using real      = double;
    using cplx      = std::complex<double>;
    using size_type = long;
    inline std::shared_ptr<spdlog::logger> log;

    // Enums
    enum class Form { SYMM, NSYM };                         // Symmetric or non-symmetric problems (complex symmetric are assumed Hermitian)
    enum class Type { REAL, CPLX };                         // Real or complex, i.e. double or std::complex<double> matrix
    enum class Side { L, R, LR };                           // Left, right or both eigenvectors (for nsym problems)
    enum class Ritz { LA, SA, LM, SM, LR, SR, LI, SI, BE }; // Choice of eigenvalue. LA is largest algebraic, and so on.
    enum class Storage { DENSE, SPARSE, TENSOR };           // Eigen Dense or Sparse, or std::vector for container
    enum class Shinv { ON, OFF };
    enum class Dephase { ON, OFF };
    enum class Vecs { ON, OFF };

    inline Ritz stringToRitz(std::string_view ritzstring) {
        if(ritzstring == "LA") return Ritz::LA;
        if(ritzstring == "SA") return Ritz::SA;
        if(ritzstring == "LM") return Ritz::LM;
        if(ritzstring == "SM") return Ritz::SM;
        if(ritzstring == "LR") return Ritz::LR;
        if(ritzstring == "SR") return Ritz::SR;
        if(ritzstring == "LI") return Ritz::LI;
        if(ritzstring == "SI") return Ritz::SI;
        if(ritzstring == "BE") return Ritz::BE;
        throw std::runtime_error("Wrong ritz string: " + std::string(ritzstring));
    }

    inline std::string_view RitzToString(Ritz ritz) {
        if(ritz == Ritz::LA) return "LA";
        if(ritz == Ritz::SA) return "SA";
        if(ritz == Ritz::LM) return "LM";
        if(ritz == Ritz::SM) return "SM";
        if(ritz == Ritz::LR) return "LR";
        if(ritz == Ritz::SR) return "SR";
        if(ritz == Ritz::LI) return "LI";
        if(ritz == Ritz::SI) return "SI";
        if(ritz == Ritz::BE) return "BE";
        throw std::runtime_error("Wrong ritz enum");
    }

}