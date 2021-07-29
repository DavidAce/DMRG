#include <complex>
#include <functional>
#include <math/hash.h>
#include <tid/tid.h>

namespace hash {
    template<typename T>
    struct is_std_complex : public std::false_type {};
    template<typename T>
    struct is_std_complex<std::complex<T>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_std_complex_v = is_std_complex<T>::value;

    inline void hash_combine([[maybe_unused]] std::size_t &seed) {}

    template<typename T, typename... Rest>
    inline void hash_combine(std::size_t &seed, const T &v, Rest... rest) {
        std::hash<T> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        hash_combine(seed, rest...);
    }

    template<typename T>
    std::size_t hash_buffer(const T *v, unsigned long size, std::size_t seed) {
        auto        t_hash = tid::tic_token("hash");
        std::size_t h      = seed;
        if constexpr(is_std_complex_v<T>) {
            for(unsigned long idx = 0; idx < size; idx++) hash_combine(h, v[idx].real(), v[idx].imag());
            return h;
        } else {
            for(unsigned long idx = 0; idx < size; idx++) hash_combine(h, v[idx]);
            return h;
        }
    }

    template std::size_t hash_buffer(const double *v, unsigned long size, std::size_t seed);
    template std::size_t hash_buffer(const std::complex<double> *v, unsigned long size, std::size_t seed);
}
