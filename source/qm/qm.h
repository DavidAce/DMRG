#pragma once

#include <complex>
namespace qm {
    using cplx = std::complex<double>;
    using real = double;
    constexpr cplx imp(0.0, 1.0);
    constexpr cplx imn(0.0, -1.0);
}