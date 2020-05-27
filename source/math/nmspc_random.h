//
// Created by david on 2016-07-24.
//

#pragma once

#include <random>
#include <iostream>
#include <complex>
#include <Eigen/Core>
#include "pcg/pcg_random.hpp"
#include <optional>

namespace rn{

    namespace internal {
        // Make a random number engine
        //    extern std::mt19937 rng;
        extern pcg64 rng;
        // Commonly used distributions
        extern std::uniform_int_distribution<int>     rand_int_01;
        extern std::uniform_real_distribution<double> rand_double_01;
        extern std::uniform_real_distribution<double> rand_double_0_2pi;
        extern std::normal_distribution<double>       normal_double_01;

    }

    //Random functions
    extern void seed(std::optional<long> n = std::nullopt);
    extern int uniform_integer_01();
    extern double  uniform_double_01();
    extern double  uniform_double_box(const double min, const double max);
    extern std::complex<double>  uniform_complex_1();
    extern double  normal(const double mean, const double std);
    extern double  log_normal(const double mean, const double std);
    extern Eigen::ArrayXd random_with_replacement(const Eigen::ArrayXd & indata);
    extern Eigen::ArrayXd random_with_replacement(const Eigen::ArrayXd & indata, const int num_choose);
    extern double gaussian_truncated(const double lowerLimit, const double upperLimit, const double mean, const double std) ;

    template<typename T, typename = std::enable_if<std::is_integral_v<T>>>
    T uniform_integer_box(const T min, const T max){
        std::uniform_int_distribution<T>  rand_int(std::min(min,max),std::max(min,max));
        return rand_int(internal::rng);
    }


    template<typename T>
    std::vector<T> uniform_unit_n_sphere(int n){
        std::vector<T> arr;
        double norm = 0.0;
        for (int i = 0; i < n; i++){
            if constexpr(std::is_same<T,std::complex<double>>::value){
                double re = internal::normal_double_01(internal::rng);
                double im = internal::normal_double_01(internal::rng);
                T cplx = T(1.0,0.0)*re + T(0.0,1.0)*im;
                arr.push_back(cplx);
                norm += re*re + im*im;
            }else{
                arr.push_back(internal::normal_double_01(internal::rng));
                norm += std::abs(arr[i]*arr[i]);

            }
        }

        norm = std::sqrt(norm);
        for (int i = 0; i < n; i++){
            arr[i] /= norm;
        }
        return arr;

    }





}
