
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
   #include <list>
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
   
       template<typename T1, typename T2, typename T3>
       std::vector<T3> range(T1 first, T2 last, T3 step){
           if (step == 0) throw std::runtime_error("Range cannot have step size zero");
           if (first > last and step > 0 ) return range(first,last,-step);
           if (first < last and step < 0 ) return range(first,last,-step);
           if (first == last) return std::vector<T3>{first};
           T3 current = first;
           std::vector<T3> vec;
           size_t num_steps = std::abs(int((last-first+step) / step));
           if(num_steps > 1000000) throw std::runtime_error("Too many steps");
           while(current <= last){
               vec.push_back(current);
               current += step;
           }
           return vec;
       }
   
       template<typename T1, typename T2, typename T3>
       std::list<T3> range_list(T1 first, T2 last, T3 step){
           std::list<T3> vec2list;
           for(auto &item : range (first,last,step)){vec2list.emplace_back(item);}
           return vec2list;
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
