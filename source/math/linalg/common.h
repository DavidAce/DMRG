#pragma once
#include <math/tenx/fwd_decl.h>
namespace Eigen {
    template<typename Idx>
    struct IndexPair;
}

namespace linalg {
    using cplx = std::complex<double>;
    using real = double;
    template<typename T>
    using EigenMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    template<typename Derived>
    using is_PlainObject = std::is_base_of<Eigen::PlainObjectBase<std::decay_t<Derived>>, std::decay_t<Derived>>;

    template<typename TA, typename TB>
    using cplx_or_real = typename std::conditional<std::is_same_v<TA, cplx> or std::is_same_v<TB, cplx>, cplx, real>::type;

    template<typename T>
    struct is_std_complex : public std::false_type {};
    template<typename T>
    struct is_std_complex<std::complex<T>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_std_complex_v = is_std_complex<T>::value;

}