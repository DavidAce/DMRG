//
// Created by david on 2016-07-24.
//

#pragma once
#include "rnd/pcg_random.hpp"
#include <Eigen/Core>
#include <complex>
#include <optional>
#include <random>
namespace rnd {

    namespace internal {
        // Make a random number engine
        //    std::mt19937 rng;
        inline pcg64 rng;
        // Commonly used distributions
        inline std::uniform_int_distribution<int>      rand_int_01(0, 1);
        inline std::uniform_real_distribution<double>  rand_double_01(0.0,1.0);
        inline std::uniform_real_distribution<double>  rand_double_0_2pi(0,2.0*M_PI);
        inline std::normal_distribution<double>        normal_double_01(0.0,1.0);
    }

    // Random functions
    extern void                 seed(std::optional<long> n = std::nullopt);
    extern int                  uniform_integer_01();
    extern double               uniform_double_01();
    extern double               uniform_double_box(double min, double max);
    extern double               uniform_double_box(double width);
    extern std::complex<double> uniform_complex_in_unit_circle();
    extern std::complex<double> uniform_complex_on_unit_circle();
    extern std::complex<double> uniform_complex_box(double real_min,double real_max, double imag_min, double imag_max);
    extern std::complex<double> uniform_complex_slice(double radius_max, double angle_min, double angle_max);
    extern double               normal(double mean, double std);
    extern double               log_normal(double mean, double std);
    extern Eigen::ArrayXd       random_with_replacement(const Eigen::ArrayXd &indata);
    extern Eigen::ArrayXd       random_with_replacement(const Eigen::ArrayXd &indata, int num_choose);
    extern double               gaussian_truncated(double lowerLimit, double upperLimit, double mean, double std);

    template<typename T, typename = std::enable_if<std::is_integral_v<T>>>
    T uniform_integer_box(const T min, const T max) {
        std::uniform_int_distribution<T> rand_int(std::min(min, max), std::max(min, max));
        return rand_int(internal::rng);
    }

    template<typename T>
    std::vector<T> uniform_unit_n_sphere(int n) {
        std::vector<T> arr;
        double         norm = 0.0;
        for(int i = 0; i < n; i++) {
            if constexpr(std::is_same<T, std::complex<double>>::value) {
                double re   = internal::normal_double_01(internal::rng);
                double im   = internal::normal_double_01(internal::rng);
                T      cplx = T(1.0, 0.0) * re + T(0.0, 1.0) * im;
                arr.push_back(cplx);
                norm += re * re + im * im;
            } else {
                arr.push_back(internal::normal_double_01(internal::rng));
                norm += std::abs(arr[i] * arr[i]);
            }
        }

        norm = std::sqrt(norm);
        for(int i = 0; i < n; i++) { arr[i] /= norm; }
        return arr;
    }

}
