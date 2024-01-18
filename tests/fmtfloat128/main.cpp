
#if defined(USE_QUADMATH)
    #include "io/fmt_custom.h"
    #include "math/f128.h"
    #include <fmt/format.h>



int main() {
    using real_t = __float128;
    using cplx_t = std::complex<real_t>;
    auto   v1  = f128_t("-3721.1244523527272214135626241245632");
    real_t v2 = f128_t("-3721.1244523527272214135626241245632");
    real_t v3 = 3;
    real_t v4 = 3.14;
    auto v5  = cplx_t(3.14, 2.3);
    cplx_t v6 = 0;
    cplx_t v7 = v1.value() * v6;
    fmt::print("v1: {:.7f}\n", v1);
    fmt::print("v2: {:.20f}\n", f128_t(v2));
    fmt::print("v3: {:>20.7f}\n", f128_t(v3));
    fmt::print("v4: {:<20.7f}\n", f128_t(v4));
    fmt::print("v4.string(5,8,f,<): {}\n", f128_t(v4).string(5,8, 'f', "<"));
    fmt::print("v5: {:>.5f}{:+>.5f}i\n", c128_t(v5).real(), c128_t(v5).imag());
    fmt::print("v6: {:>.5f}{:+>.5f}i\n", c128_t(v6).real(), c128_t(v6).imag());
    fmt::print("v7: {:>.5f}{:+>.5f}i\n", c128_t(v7).real(), c128_t(v7).imag());

    return 0;
}
#else
int main() { return 0; }
#endif