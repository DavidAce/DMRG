#pragma once
#include "math/tenx/fwd_decl.h"
#include "math/float.h"

namespace Eigen {
    template<typename Idx>
    struct IndexPair;
}

namespace linalg {
    template<typename T>
    using EigenMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    template<typename Derived>
    using is_PlainObject = std::is_base_of<Eigen::PlainObjectBase<std::decay_t<Derived>>, std::decay_t<Derived>>;

    template<typename TA, typename TB>
    using cplx_or_real = typename std::conditional<std::is_same_v<TA, cx64> or std::is_same_v<TB, cx64>, cx64, fp64>::type;

    template<typename T>
    struct is_std_complex : public std::false_type {};
    template<typename T>
    struct is_std_complex<std::complex<T>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_std_complex_v = is_std_complex<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_value_type : public std::false_type {};
    template<typename T>
    struct has_value_type<T, std::void_t<typename T::value_type>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_value_type_v = has_value_type<T>::value;

}