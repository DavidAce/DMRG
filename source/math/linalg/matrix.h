#pragma once
#include "common.h"
#include <Eigen/Core>

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
        else {
            return kronecker(A.derived().eval(), B.derived().eval(), mirror);
        }
    }

    template<typename T>
    std::string to_string(const Eigen::EigenBase<T> &m, const Eigen::IOFormat &f = Eigen::IOFormat(4, 0, ", ", "\n", "  [", "]")) {
        std::stringstream ss;
        ss << m.derived().format(f);
        return ss.str();
    }

    template<typename T>
    std::string to_string(const Eigen::EigenBase<T> &m, int prec = 1, int flags = 0, const std::string &sep = ", " ) {
        Eigen::IOFormat f(prec,flags,sep, "\n", "  [", "]");
        return to_string(m, f);
    }
}
