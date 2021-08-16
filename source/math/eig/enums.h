#pragma once
#include <complex>
#include <optional>
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

    enum class PrimmeMethod{
        PRIMME_DEFAULT_METHOD,
        PRIMME_DYNAMIC,
        PRIMME_DEFAULT_MIN_TIME,
        PRIMME_DEFAULT_MIN_MATVECS,
        PRIMME_Arnoldi,
        PRIMME_GD,
        PRIMME_GD_plusK,
        PRIMME_GD_Olsen_plusK,
        PRIMME_JD_Olsen_plusK,
        PRIMME_RQI,
        PRIMME_JDQR,
        PRIMME_JDQMR,
        PRIMME_JDQMR_ETol,
        PRIMME_STEEPEST_DESCENT,
        PRIMME_LOBPCG_OrthoBasis,
        PRIMME_LOBPCG_OrthoBasis_Window,
    };

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

    inline PrimmeMethod stringToMethod(std::optional<std::string> methodstring) {
        if(not methodstring.has_value()) return PrimmeMethod::PRIMME_DEFAULT_MIN_MATVECS;
        if(methodstring.value() == "PRIMME_DEFAULT_METHOD") return PrimmeMethod::PRIMME_DEFAULT_METHOD;
        if(methodstring.value() == "PRIMME_DYNAMIC") return PrimmeMethod::PRIMME_DYNAMIC;
        if(methodstring.value() == "PRIMME_DEFAULT_MIN_TIME") return PrimmeMethod::PRIMME_DEFAULT_MIN_TIME;
        if(methodstring.value() == "PRIMME_DEFAULT_MIN_MATVECS") return PrimmeMethod::PRIMME_DEFAULT_MIN_MATVECS;
        if(methodstring.value() == "PRIMME_Arnoldi") return PrimmeMethod::PRIMME_Arnoldi;
        if(methodstring.value() == "PRIMME_GD") return PrimmeMethod::PRIMME_GD;
        if(methodstring.value() == "PRIMME_GD_plusK") return PrimmeMethod::PRIMME_GD_plusK;
        if(methodstring.value() == "PRIMME_GD_Olsen_plusK") return PrimmeMethod::PRIMME_GD_Olsen_plusK;
        if(methodstring.value() == "PRIMME_JD_Olsen_plusK") return PrimmeMethod::PRIMME_JD_Olsen_plusK;
        if(methodstring.value() == "PRIMME_RQI") return PrimmeMethod::PRIMME_RQI;
        if(methodstring.value() == "PRIMME_JDQR") return PrimmeMethod::PRIMME_JDQR;
        if(methodstring.value() == "PRIMME_JDQMR") return PrimmeMethod::PRIMME_JDQMR;
        if(methodstring.value() == "PRIMME_JDQMR_ETol") return PrimmeMethod::PRIMME_JDQMR_ETol;
        if(methodstring.value() == "PRIMME_STEEPEST_DESCENT") return PrimmeMethod::PRIMME_STEEPEST_DESCENT;
        if(methodstring.value() == "PRIMME_LOBPCG_OrthoBasis") return PrimmeMethod::PRIMME_LOBPCG_OrthoBasis;
        if(methodstring.value() == "PRIMME_LOBPCG_OrthoBasis_Window") return PrimmeMethod::PRIMME_LOBPCG_OrthoBasis_Window;
        throw std::runtime_error("Wrong method string: " + methodstring.value());
    }
    inline std::string_view MethodToString(std::optional<PrimmeMethod> method) {
        if(not method.has_value()) return "PRIMME_DEFAULT_MIN_MATVECS";
        if(method.value() == PrimmeMethod::PRIMME_DEFAULT_METHOD) return "PRIMME_DEFAULT_METHOD";
        if(method.value() == PrimmeMethod::PRIMME_DYNAMIC) return "PRIMME_DYNAMIC";
        if(method.value() == PrimmeMethod::PRIMME_DEFAULT_MIN_TIME) return "PRIMME_DEFAULT_MIN_TIME";
        if(method.value() == PrimmeMethod::PRIMME_DEFAULT_MIN_MATVECS) return "PRIMME_DEFAULT_MIN_MATVECS";
        if(method.value() == PrimmeMethod::PRIMME_Arnoldi) return "PRIMME_Arnoldi";
        if(method.value() == PrimmeMethod::PRIMME_GD) return "PRIMME_GD";
        if(method.value() == PrimmeMethod::PRIMME_GD_plusK) return "PRIMME_GD_plusK";
        if(method.value() == PrimmeMethod::PRIMME_GD_Olsen_plusK) return "PRIMME_GD_Olsen_plusK";
        if(method.value() == PrimmeMethod::PRIMME_JD_Olsen_plusK) return "PRIMME_JD_Olsen_plusK";
        if(method.value() == PrimmeMethod::PRIMME_RQI) return "PRIMME_RQI";
        if(method.value() == PrimmeMethod::PRIMME_JDQR) return "PRIMME_JDQR";
        if(method.value() == PrimmeMethod::PRIMME_JDQMR) return "PRIMME_JDQMR";
        if(method.value() == PrimmeMethod::PRIMME_JDQMR_ETol) return "PRIMME_JDQMR_ETol";
        if(method.value() == PrimmeMethod::PRIMME_STEEPEST_DESCENT) return "PRIMME_STEEPEST_DESCENT";
        if(method.value() == PrimmeMethod::PRIMME_LOBPCG_OrthoBasis) return "PRIMME_LOBPCG_OrthoBasis";
        if(method.value() == PrimmeMethod::PRIMME_LOBPCG_OrthoBasis_Window) return "PRIMME_LOBPCG_OrthoBasis_Window";
        return "PRIMME_DEFAULT_MIN_MATVECS";
    }
}