#pragma once
#include "common.h"
#include <Eigen/Core>
#include "math/float.h"

namespace linalg::matrix {
    template<auto N, auto M>
    constexpr auto multiply() {
        if constexpr(N == -1 or M == -1)
            return -1;
        else
            return N * M;
    }

    template<typename DA, typename DB>
    using KroneckerResultType = Eigen::Matrix<cplx_or_real<typename DA::Scalar, typename DB::Scalar>, multiply<DA::RowsAtCompileTime, DB::RowsAtCompileTime>(),
                                              multiply<DA::ColsAtCompileTime, DB::ColsAtCompileTime>(), Eigen::ColMajor>;

    template<typename DerivedA, typename DerivedB>
    extern KroneckerResultType<DerivedA, DerivedB> kronecker(const Eigen::PlainObjectBase<DerivedA> &A, const Eigen::PlainObjectBase<DerivedB> &B,
                                                             bool mirror = false);

    template<typename DerivedA, typename DerivedB>
    auto kronecker(const Eigen::EigenBase<DerivedA> &A, const Eigen::EigenBase<DerivedB> &B, bool mirror = false) {
        if constexpr(is_PlainObject<DerivedA>::value and is_PlainObject<DerivedB>::value)
            return kronecker(A, B, mirror);
        else { return kronecker(A.derived().eval(), B.derived().eval(), mirror); }
    }

    template<typename T>
    std::string to_string(const Eigen::EigenBase<T> &m, const Eigen::IOFormat &f = Eigen::IOFormat(4, 0, ", ", "\n", "  [", "]")) {
        using Scalar = typename T::Scalar;
        if constexpr(std::is_floating_point_v<Scalar>) {
            std::stringstream ss;
            ss << m.derived().format(f);
            return ss.str();
        } else if constexpr(is_std_complex_v<Scalar>) {
            using Inner = typename Scalar::value_type;
            if constexpr(std::is_floating_point_v<Inner>) {
                std::stringstream ss;
                ss << m.derived().format(f);
                return ss.str();
            }
#if defined(USE_QUADMATH)
            else if constexpr(std::is_same_v<Inner, __float128>) {
                std::string s = "[";
                std::string rbuf;
                std::string ibuf;
                auto        bufsize = std::max<size_t>(64, static_cast<size_t>(f.precision) + 34);
                rbuf.resize(bufsize);
                ibuf.resize(bufsize);
                s.reserve(bufsize * m.derived().size());
                std::string fmt = "%." + std::to_string(f.precision) + "Qg";
                for(long i = 0; i < m.derived().rows(); ++i) {
                    s += "[";
                    for(long j = 0; j < m.derived().cols(); ++j) {
                        quadmath_snprintf(rbuf.data(), bufsize, fmt.data(), m.derived()(i, j).real());
                        quadmath_snprintf(ibuf.data(), bufsize, fmt.data(), m.derived()(i, j).imag());
                        s += "(" + rbuf + "," + ibuf + "i)" + f.coeffSeparator;
                    }
                    s += "]";
                    if(i + 1 != m.derived().rows()) s += '\n';
                }
                s += "]";
                return s;
            }
#endif
            else {
                throw std::runtime_error("to_string: this complex type is not implemented");
            }

        }
#if defined(USE_QUADMATH)
        else if constexpr(std::is_same_v<Scalar, __float128>) {
            std::string s = "[";
            std::string buf;
            auto        bufsize = std::max<size_t>(64, static_cast<size_t>(f.precision) + 34);
            buf.resize(bufsize);
            s.reserve(bufsize * m.derived().size());
            std::string fmt = "%." + std::to_string(f.precision) + "Qf";
            for(long i = 0; i < m.derived().rows(); ++i) {
                s += "[";
                for(long j = 0; j < m.derived().cols(); ++j) {
                    quadmath_snprintf(buf.data(), bufsize, fmt.data(), m.derived()(i, j));
                    s += buf + f.coeffSeparator;
                }
                s += "]";
                if(i + 1 != m.derived().rows()) s += '\n';
            }
            s += "]";
            return s;
        }
#endif
        else
            throw std::runtime_error("to_string: type is not implemented");
    }
    template<typename T>
    std::string to_string(const Eigen::EigenBase<T> &m, int prec, int flags, const std::string &sep = ", ") {
        Eigen::IOFormat f(prec, flags, sep, "\n", "  [", "]");
        return to_string(m, f);
    }
}
