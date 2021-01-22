//
// Created by david on 2016-07-24.
//

#pragma once
#include <complex>
#include <optional>
#include <random>



namespace rnd {


    // Random functions
    extern void                 seed(std::optional<long> n = std::nullopt);
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
    extern std::vector<int>     random_with_replacement(const std::vector<int> & indata);
    extern std::vector<int>     random_with_replacement(const std::vector<int> & indata, size_t num_choose);
    extern double               gaussian_truncated(double lowerLimit, double upperLimit, double mean, double std);

    template<typename T>
    T uniform_integer_box(T min, T max);

    template<typename T>
    std::vector<T> uniform_unit_n_sphere(size_t n);

    template<typename T>
    void shuffle(T & list);

}
