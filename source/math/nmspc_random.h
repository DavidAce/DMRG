//
// Created by david on 2016-07-24.
//

#pragma once

#include <random>
#include <iostream>
#include <complex>
#include <Eigen/Core>
#include "pcg/pcg_random.hpp"


namespace rn{
    //Random functions
//    extern std::mt19937 rng;


    // Make a random number engine
    extern pcg32 rng;

    extern void seed(unsigned long n);
    extern int uniform_integer_1();
    extern int uniform_integer(const int min, const int max);
    extern double  uniform_double_1();
    extern double  uniform_double(const double min, const double max);
    extern std::complex<double>  uniform_complex_1();
    extern double  normal(const double mean, const double std);
    extern double  log_normal(const double mean, const double std);
    extern Eigen::ArrayXd random_with_replacement(const Eigen::ArrayXd & indata);
    extern Eigen::ArrayXd random_with_replacement(const Eigen::ArrayXd & indata, const int num_choose);
    extern double gaussian_truncated(const double lowerLimit, const double upperLimit, const double mean, const double std) ;

    template<typename T>
    inline std::vector<T> uniform_unit_n_sphere(int n){
        std::vector<T> arr;
        std::normal_distribution<double> distribution(0.0,1.0);
        double norm = 0.0;
        for (int i = 0; i < n; i++){
            if constexpr(std::is_same<T,std::complex<double>>::value){
                double re = distribution(rng);
                double im = distribution(rng);
                T cplx = T(1.0,0.0)*re + T(0.0,1.0)*im;
                arr.push_back(cplx);
                norm += re*re + im*im;
            }else{
                arr.push_back(distribution(rng));
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
