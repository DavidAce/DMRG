
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_general_nmspc_random_numbers.h:

Program Listing for File nmspc_random_numbers.h
===============================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_general_nmspc_random_numbers.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/general/nmspc_random_numbers.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2016-07-24.
   //
   
   #ifndef NMSPC_RANDOM_NUMBERS_H
   #define NMSPC_RANDOM_NUMBERS_H
   #include <random>
   #include <iostream>
   #include <complex>
   #include <Eigen/Core>
   
   namespace rn{
       //Random functions
       extern std::mt19937 rng;
       extern void seed(unsigned long n);
       extern int __attribute__((hot)) uniform_integer_1();
       extern int __attribute__((hot)) uniform_integer(const int min, const int max);
       extern double __attribute__((hot)) uniform_double_1();
       extern double __attribute__((hot)) uniform_double(const double min, const double max);
       extern std::complex<double> __attribute__((hot)) uniform_complex_1();
       extern double __attribute__((hot)) log_normal(const double mean, const double std);
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
   
   #endif //NMSPC_RANDOM_NUMBERS_H
