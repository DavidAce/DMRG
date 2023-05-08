#pragma once
#ifndef DMRG_FLOAT_H
    #define DMRG_FLOAT_H
    #include <complex>
// This header defines the floating point types used throughout this program

using real = double;
using cplx = std::complex<real>;

    #if defined(USE_QUADMATH)
        #include <quadmath.h>
using real_t = __float128;
using cplx_t = std::complex<real_t>;
    #else
using real_t = long double;
using cplx_t = std::complex<real_t>;
    #endif

template<typename Lhs, typename Rhs>
bool cmp_t(Lhs lhs, Rhs rhs) {
    if constexpr(std::is_same_v<Lhs, cplx_t> or std::is_same_v<Rhs, cplx_t>)
        return cplx_t(lhs) == cplx_t(rhs);
    else if constexpr(std::is_same_v<Lhs, real_t> or std::is_same_v<Rhs, real_t>)
        return real_t(lhs) == real_t(rhs);
    else
        return lhs == rhs;
}

template<typename T>
auto abs_t(T val) {
    if constexpr(std::is_arithmetic_v<T>) return std::abs(val);
    #if defined(USE_QUADMATH)
    else if constexpr(std::is_same_v<T, cplx_t>) {
        __complex128 x;
        __real__ x = val.real();
        __imag__ x = val.imag();
        return cabsq(x);
    } else if constexpr(std::is_same_v<T, real_t>)
        return fabsq(val);
    #else
    static_assert(std::is_arithmetic_v<T>);
    #endif
}

#endif