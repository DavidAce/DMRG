#include "solution.h"
namespace eig {

    template<typename Scalar, Side side>
    std::vector<Scalar> &solution::get_eigvecs() const {
        static_assert(side != Side::LR and "Cannot get both L/R eigvecs simultaneusly");
        if constexpr(std::is_same_v<Scalar, real>) {
            build_eigvecs_real();
            if constexpr(side == Side::R) return eigvecsR_real;
            if constexpr(side == Side::L) return eigvecsL_real;
            //            if constexpr(side == Side::LR) return std::pair(eigvecsL_real, eigvecsR_real);
        }
        if constexpr(std::is_same_v<Scalar, cplx>) {
            build_eigvecs_cplx();
            if constexpr(side == Side::R) return eigvecsR_cplx;
            if constexpr(side == Side::L) return eigvecsL_cplx;
            //            if constexpr(side == Side::LR) return std::pair(eigvecsL_cplx, eigvecsR_cplx);
        }
    }

    template std::vector<eig::real> &solution::get_eigvecs<eig::real, Side::L>() const;
    template std::vector<eig::cplx> &solution::get_eigvecs<eig::cplx, Side::L>() const;
    template std::vector<eig::real> &solution::get_eigvecs<eig::real, Side::R>() const;
    template std::vector<eig::cplx> &solution::get_eigvecs<eig::cplx, Side::R>() const;

    template<typename Scalar>
    std::vector<Scalar> &solution::get_eigvecs(Side side) const {
        if(side == Side::R) return get_eigvecs<Scalar, Side::R>();
        if(side == Side::L) return get_eigvecs<Scalar, Side::L>();
        throw std::runtime_error("Cannot return both L and R eigenvectors");
    }
    template std::vector<eig::real> &solution::get_eigvecs<eig::real>(Side side) const;
    template std::vector<eig::cplx> &solution::get_eigvecs<eig::cplx>(Side side) const;

    template<typename Scalar, Form form, Side side>
    std::vector<Scalar> &solution::get_eigvecs() const {
        if constexpr(std::is_same<real, Scalar>::value) return get_eigvecs<form, Type::REAL, side>();
        if constexpr(std::is_same<cplx, Scalar>::value) return get_eigvecs<form, Type::CPLX, side>();
    }
    template std::vector<eig::real> &solution::get_eigvecs<eig::real, Form::SYMM, Side::L>() const;
    template std::vector<eig::real> &solution::get_eigvecs<eig::real, Form::SYMM, Side::R>() const;
    template std::vector<eig::cplx> &solution::get_eigvecs<eig::cplx, Form::SYMM, Side::L>() const;
    template std::vector<eig::cplx> &solution::get_eigvecs<eig::cplx, Form::SYMM, Side::R>() const;
    template std::vector<eig::cplx> &solution::get_eigvecs<eig::cplx, Form::NSYM, Side::L>() const;
    template std::vector<eig::cplx> &solution::get_eigvecs<eig::cplx, Form::NSYM, Side::R>() const;

    template<typename Scalar>
    std::vector<Scalar> &solution::get_eigvals() const {
        if constexpr(std::is_same_v<Scalar, real>) {
            build_eigvals_real();
            return eigvals_real;
        }
        if constexpr(std::is_same_v<Scalar, cplx>) {
            build_eigvals_cplx();
            return eigvals_cplx;
        }
    }

    template std::vector<eig::real> &solution::get_eigvals<eig::real>() const;
    template std::vector<eig::cplx> &solution::get_eigvals<eig::cplx>() const;

    const std::vector<double> &solution::get_resnorms() const { return meta.residual_norms; }

    void solution::reset() {
        eigvals_real.clear();
        eigvals_imag.clear();
        eigvals_cplx.clear();
        eigvecsR_real.clear();
        eigvecsR_imag.clear();
        eigvecsL_real.clear();
        eigvecsL_imag.clear();
        eigvecsR_cplx.clear();
        eigvecsL_cplx.clear();
        meta = Meta();
    }

    bool solution::eigvecs_are_real() const { return meta.form == Form::SYMM and meta.type == Type::REAL; }

    bool solution::eigvals_are_real() const { return meta.form == Form::SYMM; }

    std::type_index solution::get_eigvecs_type() const {
        if(eigvecs_are_real())
            return typeid(real);
        else
            return typeid(cplx);
    }

    void solution::build_eigvecs_cplx() const {
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

    void solution::build_eigvecs_real() const {
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

    void solution::build_eigvals_cplx() const {
        bool build_cplx = eigvals_cplx.empty() and not eigvals_imag.empty() and eigvals_real.size() == eigvals_imag.size();
        if(build_cplx) {
            eigvals_cplx.resize(eigvals_real.size());
            for(size_t i = 0; i < eigvals_real.size(); i++) eigvals_cplx[i] = std::complex<double>(eigvals_real[i], eigvals_imag[i]);
            eigvals_real.clear();
            eigvals_imag.clear();
        }
    }

    void solution::build_eigvals_real() const {
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
}