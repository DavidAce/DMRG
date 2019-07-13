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

/*!
 *  \namespace Math
 *  \brief Small convenience-type math functions like modulo and numerical integration using GSL.
 *  \tableofcontents
 */
namespace math

{
    /*! \brief MatLab-style modulo operator
    *   \param x first number
    *   \param y second number
    *   \return modulo of x and y. Example,  <code> mod(-0.5,10)  = 9.5 </code>, instead of <code> -0.5 </code>  as given by x%y.
    */
    template<typename T1, typename T2>
    inline auto mod(const T1 x, const T2 y)
    {
        return (x % y + y) % y;
    }


    /*! \brief MatLab-style linearly spaced array
    *   \param num number of linearly spaced values
    *   \param min minimum value in range
    *   \param max maximum value in range
    *   \return std::vector<T2>. Example,  <code> Linspaced(5,1,5) </code> gives a std::vector<int>: <code> [1,2,3,4,5] </code>
    */
    template <typename T1, typename T2>
    inline std::vector<T2> LinSpaced(T1 num, T2 min, T2 max )

    {
        static_assert(std::is_integral<T1>::value && "math::LinSpaced -- Given type is not integral!");
        Eigen::Array<T2, Eigen::Dynamic, 1> temp =  Eigen::Array<T2, Eigen::Dynamic, 1> :: LinSpaced(num, min, max);
        return std::vector<T2> (temp.data(), temp.data() + temp.size());
    }


    /*! \brief Product operator for containers such as vector
    *   \param in a vector, array or any 1D container with "<code> .data() </code>" method.
    *   \param from first element to multiply
    *   \param to last element to multiply
    *   \return std::vector<T2>. Example, let <code> my_vector = {1,2,3,4}</code>. Then <code> prod(my_vector,0,3) = 24 </code>.
    */
    template<typename Input, typename From, typename To>
    auto prod(const Input &in, const From from, const To to)
    {
        return std::accumulate(in.data() + from, in.data()+to,1,std::multiplies<>());
    }

}



#endif //TRAINING_FUNCS_H
