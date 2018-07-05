//
// Created by david on 2016-07-24.
//

#ifndef WL_NMSPC_RANDOM_NUMBERS_H
#define WL_NMSPC_RANDOM_NUMBERS_H
#include <random>
#include <iostream>
#include <complex>
#include <Eigen/Core>

namespace rn{
    //typedef std::mt19937 RNGType;
    //RNGType rng;
    //Random functions
    extern std::mt19937 rng;
    extern void seed(unsigned long n);

    inline int __attribute__((hot)) uniform_integer_1(){
        std::uniform_int_distribution<>  rand_int(0, 1);
        return rand_int(rng);
    }

    inline int __attribute__((hot)) uniform_integer(const int min, const int max){
        std::uniform_int_distribution<>  rand_int(min,max);
        return rand_int(rng);
    }

    inline double __attribute__((hot)) uniform_double_1(){
        std::uniform_real_distribution<>  rand_real(0,1);
        return rand_real(rng);
    }

    inline double __attribute__((hot)) uniform_double(const double min, const double max){
        std::uniform_real_distribution<>  rand_real(std::min(min,max),std::max(min,max));
        return rand_real(rng);
    }


    inline std::complex<double> __attribute__((hot)) uniform_complex_1(){
        std::uniform_real_distribution<>  rand_real(0,2.0*M_PI);
        return std::polar(1.0,rand_real(rng));
    }

    template<typename T>
    inline std::vector<T> uniform_unit_n_sphere(int n){
        std::vector<T> arr;
        std::normal_distribution<double> distribution(0.0,1.0);
        double norm = 0.0;
        for (int i = 0; i < n; i++){
            if constexpr(std::is_same<T,std::complex<double>>::value){
                double re = distribution(rng);
                double im = distribution(rng);
                std::complex<double> cplx = 1.0*re + 1.0i*im;
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


    extern Eigen::ArrayXd random_with_replacement(const Eigen::ArrayXd & indata);
    extern Eigen::ArrayXd random_with_replacement(const Eigen::ArrayXd & indata, const int num_choose);

    extern double gaussian_truncated(const double lowerLimit, const double upperLimit, const double mean, const double std) ;
}

#endif //WL_NMSPC_RANDOM_NUMBERS_H
