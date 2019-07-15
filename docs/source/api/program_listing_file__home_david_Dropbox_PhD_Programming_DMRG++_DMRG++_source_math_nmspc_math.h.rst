
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_math_nmspc_math.h:

Program Listing for File nmspc_math.h
=====================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_math_nmspc_math.h>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/math/nmspc_math.h``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   #ifndef TRAINING_FUNCS_H
   #define TRAINING_FUNCS_H
   
   
   #include <iostream>
   #include <iterator>
   #include <vector>
   #include <iostream>
   #include <cmath>
   
   #include <functional>
   #include <memory>
   #include <utility>
   //#include <gsl/gsl_errno.h>
   //#include <gsl/gsl_integration.h>
   #include <numeric>
   #include <Eigen/Core>
   
   namespace math
   
   {
       template<typename T1, typename T2>
       inline auto mod(const T1 x, const T2 y)
       {
           return (x % y + y) % y;
       }
   
   
       template <typename T1, typename T2>
       inline std::vector<T2> LinSpaced(T1 num, T2 min, T2 max )
   
       {
           static_assert(std::is_integral<T1>::value && "math::LinSpaced -- Given type is not integral!");
           Eigen::Array<T2, Eigen::Dynamic, 1> temp =  Eigen::Array<T2, Eigen::Dynamic, 1> :: LinSpaced(num, min, max);
           return std::vector<T2> (temp.data(), temp.data() + temp.size());
       }
   
   
       template<typename Input, typename From, typename To>
       auto prod(const Input &in, const From from, const To to)
       {
           return std::accumulate(in.data() + from, in.data()+to,1,std::multiplies<>());
       }
   
   }
   
   
   
   #endif //TRAINING_FUNCS_H
