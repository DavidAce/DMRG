#pragma once
#include <complex>
#include <optional>
#include <pcg_random.hpp>
#include <vector>

namespace rnd {

    enum class dist { uniform, normal, lognormal };
    constexpr std::string_view enum2sv(const dist &d);
    constexpr dist             sv2enum(std::string_view d);

    namespace internal {
        // The random number engine
        inline pcg64 rng;

    }
    // Random functions
    void                        seed(std::optional<long> n = std::nullopt);
    extern int                  uniform_integer_01();
    extern double               uniform_double_01();
    extern double               uniform_double_box(double min, double max);
    extern double               uniform_double_box(double halfwidth);
    extern std::complex<double> uniform_complex_in_unit_circle();
    extern std::complex<double> uniform_complex_on_unit_circle();
    extern std::complex<double> uniform_complex_box(double real_min, double real_max, double imag_min, double imag_max);
    extern std::complex<double> uniform_complex_slice(double radius_max, double angle_min, double angle_max);
    extern double               normal(double mean, double std);
    extern double               log_normal(double mean, double std);
    extern std::vector<int>     random_with_replacement(const std::vector<int> &indata);
    extern std::vector<int>     random_with_replacement(const std::vector<int> &indata, size_t num_choose);
    extern double               gaussian_truncated(double lowerLimit, double upperLimit, double mean, double std);

    template<typename T>
    T uniform_integer_box(T min, T max);

    template<typename T>
    std::vector<T> uniform_unit_n_sphere(size_t n);

    template<typename T>
    void shuffle(T &list);

    extern double              random(dist d, double mean, double width);
    extern double              random(std::string_view distribution, double mean, double width);
    extern std::vector<double> random(std::string_view distribution, double mean, double width, size_t num);
    extern std::vector<double> random(std::string_view distribution, double mean, double width, const std::vector<double> &weights);

}
