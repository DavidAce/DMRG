#pragma once
#include <complex>

namespace eig {
    using real      = double;
    using cplx      = std::complex<double>;
    using size_type = long;

    // Enums
    enum class Lib { ARPACK, PRIMME }; // Choose the underlying library
    enum class Form { SYMM, NSYM };    // Symmetric or non-symmetric problems (complex symmetric are assumed Hermitian)
    enum class Type { REAL, CPLX };    // Real or complex, i.e. double or std::complex<double> matrix
    enum class Side { L, R, LR };      // Left, right or both eigenvectors (for nsym problems)
    enum class Ritz {
        LA,
        SA,
        LM,
        SM,
        LR,
        SR,
        LI,
        SI,
        BE,
        primme_smallest,
        primme_largest,
        primme_closest_geq,
        primme_closest_leq,
        primme_closest_abs,
        primme_largest_abs
    };                                            // Choice of eigenvalue. LA is largest algebraic, and so on.
    enum class Storage { DENSE, SPARSE, MPS }; // Eigen Dense or Sparse, or std::vector for container
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
        if(ritzstring == "primme_smallest") return Ritz::primme_smallest;
        if(ritzstring == "primme_largest") return Ritz::primme_largest;
        if(ritzstring == "primme_closest_geq") return Ritz::primme_closest_geq;
        if(ritzstring == "primme_closest_leq") return Ritz::primme_closest_leq;
        if(ritzstring == "primme_closest_abs") return Ritz::primme_closest_abs;
        if(ritzstring == "primme_largest_abs") return Ritz::primme_largest_abs;
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
        if(ritz == Ritz::primme_smallest) return "primme_smallest";
        if(ritz == Ritz::primme_largest) return "primme_largest";
        if(ritz == Ritz::primme_closest_geq) return "primme_closest_geq";
        if(ritz == Ritz::primme_closest_leq) return "primme_closest_leq";
        if(ritz == Ritz::primme_closest_abs) return "primme_closest_abs";
        if(ritz == Ritz::primme_largest_abs) return "primme_largest_abs";
        throw std::runtime_error("Wrong ritz enum");
    }

}