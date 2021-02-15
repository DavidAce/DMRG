#pragma once

#include <cmath>
#include <complex>
#include <iterator>
#include <numeric>
#include <vector>

/*!
 *  \namespace num
 *  \brief Small convenience-type num functions, like modulo
 *  \tableofcontents
 */

namespace num {


    /*! \brief MatLab-style modulo operator
     *   \param x first number
     *   \param y second number
     *   \return modulo of x and y. Example, <code> mod(7,2)  = 1 </code> but <code> mod(-0.5,10)  = 9.5 </code>, instead of <code> -0.5 </code>  as given by
     * x%y.
     */
    template<typename T>
    inline T mod(const T x, const T y) {
        if constexpr(std::is_integral_v<T>)
            return (x % y + y) % y;
        else
            return std::fmod((std::fmod(x, y) + y), y);
    }

    /*! \brief Python-style range generator, i.e. not-including "last"
     *   \return Range of T's. Example, <code> range(0,8,2) </code> gives a std::vector<int>: <code> [0,2,4,6] </code>
     */

    namespace internal{
        template<typename TA,typename TB>
        using int_or_dbl =
        typename std::conditional<std::is_floating_point_v<TA> or std::is_floating_point_v<TB>, double, int>::type;
    }

    template<typename T = int, typename T1, typename T2, typename T3 = internal::int_or_dbl<T1,T2>>
    std::vector<T> range(T1 first, T2 last, T3 step = static_cast<T3>(1)) {
        if(step == 0) throw std::runtime_error("Range cannot have step size zero");
        if constexpr(std::is_signed_v<T3>){
            if(static_cast<T3>(first) > static_cast<T3>(last) and step > 0) return range<T>(first, last, -step);
            if(static_cast<T3>(first) < static_cast<T3>(last) and step < 0) return range<T>(first, last, -step);
        }else{
            if(static_cast<T3>(first) > static_cast<T3>(last)) throw std::runtime_error("Range of unsigned step type cannot have first > last");
        }
        if(static_cast<T3>(first) == static_cast<T3>(last)) return {};
        if(static_cast<T3>(first) + step == static_cast<T3>(last)) return std::vector<T>{static_cast<T>(first)};

        auto num_steps = static_cast<size_t>(
            std::abs<double>((static_cast<double>(last) - static_cast<double>(first) + static_cast<double>(step)) / static_cast<double>(step)));
        if(num_steps > 10000000) throw std::runtime_error("Too many steps for range");

        std::vector<T> vec;
        vec.reserve(num_steps);
        for(T3 current = static_cast<T3>(first); current < static_cast<T3>(last); current += step) vec.emplace_back(current);
        return vec;
    }

    /*! \brief MatLab-style linearly spaced array
     *   \param num number of linearly spaced values
     *   \param a first value in range
     *   \param b last value in range
     *   \return std::vector<T2>. Example,  <code> Linspaced(5,1,5) </code> gives a std::vector<int>: <code> [1,2,3,4,5] </code>
     */
    inline std::vector<double> LinSpaced(std::size_t N, double a, double b) {
        double              h = (b - a) / static_cast<double>(N - 1);
        std::vector<double> xs(N);
        double              val = a;
        for(auto &x : xs) {
            x = val;
            val += h;
        }
        return xs;
    }

    inline std::vector<double> LogSpaced(std::size_t N, double a, double b, double base = 10.0) {
        if(a <= 0) throw std::range_error("a must be positive");
        if(b <= 0) throw std::range_error("b must be positive");
        double              loga   = std::log(a) / std::log(base);
        double              logb   = std::log(b) / std::log(base);
        double              h      = (logb - loga) / static_cast<double>(N - 1);
        double              factor = std::pow(base, h);
        double              val    = std::pow(base, loga);
        std::vector<double> xs(N);
        for(auto &x : xs) {
            x = val;
            val *= factor;
        }
        return xs;
    }

    /*! \brief Product operator for containers such as vector
     *   \param in a vector, array or any 1D container with "<code> .data() </code>" method.
     *   \param from first element to multiply
     *   \param to last element to multiply
     *   \return std::vector<T2>. Example, let <code> my_vector = {1,2,3,4}</code>. Then <code> prod(my_vector,0,3) = 24 </code>.
     */
    template<typename Input, typename From, typename To>
    auto prod(const Input &in, const From from, const To to) {
        return std::accumulate(in.data() + from, in.data() + to, 1, std::multiplies<>());
    }

    /*! \brief Checks if multiple values are equal to each other
     *   \param args any number of values
     *   \return bool, true if all args are equal
     */
    template<typename First, typename... T>
    bool all_equal(First &&first, T &&...t) {
        return ((first == t) && ...);
    }

    template<typename R, typename T>
    R next_power_of_two(T val) {
        return static_cast<R>(std::pow<long>(2, static_cast<long>(std::ceil(std::log2(std::real(val))))));
    }

}
