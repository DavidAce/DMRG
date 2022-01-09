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
#if defined(NDEBUG)
    static constexpr bool ndebug = true;
#else
    static constexpr bool ndebug = false;
#endif
    namespace internal {
        template<typename T>
        struct is_reference_wrapper : std::false_type {};

        template<typename T>
        struct is_reference_wrapper<std::reference_wrapper<T>> : std::true_type {};

        template<typename T>
        inline constexpr bool is_reference_wrapper_v = is_reference_wrapper<T>::value;
    }

    // Safe integer comparison functions from C++20

    template<class T, class U>
    [[nodiscard]] constexpr bool cmp_equal(T t, U u) noexcept {
        if constexpr(internal::is_reference_wrapper_v<T>)
            return cmp_equal(t.get(), u);
        else if constexpr(internal::is_reference_wrapper_v<U>)
            return cmp_equal(t, u.get());
        else if constexpr(std::is_floating_point_v<T> or std::is_floating_point_v<U>)
            return t == u;
        else {
            using UT = std::make_unsigned_t<T>;
            using UU = std::make_unsigned_t<U>;
            if constexpr(std::is_signed_v<T> == std::is_signed_v<U>)
                return t == u;
            else if constexpr(std::is_signed_v<T>)
                return t < 0 ? false : UT(t) == u;
            else
                return u < 0 ? false : t == UU(u);
        }
    }

    template<class T, class U>
    [[nodiscard]] constexpr bool cmp_not_equal(T t, U u) noexcept {
        return !cmp_equal(t, u);
    }

    template<class T, class U>
    [[nodiscard]] constexpr bool cmp_less(T t, U u) noexcept {
        if constexpr(internal::is_reference_wrapper_v<T>)
            return cmp_less(t.get(), u);
        else if constexpr(internal::is_reference_wrapper_v<U>)
            return cmp_less(t, u.get());
        else if constexpr(std::is_floating_point_v<T> or std::is_floating_point_v<U>)
            return t < u;
        else {
            using UT = std::make_unsigned_t<T>;
            using UU = std::make_unsigned_t<U>;
            if constexpr(std::is_signed_v<T> == std::is_signed_v<U>)
                return t < u;
            else if constexpr(std::is_signed_v<T>)
                return t < 0 ? true : UT(t) < u;
            else
                return u < 0 ? false : t < UU(u);
        }
    }

    template<class T, class U>
    [[nodiscard]] constexpr bool cmp_greater(T t, U u) noexcept {
        return cmp_less(u, t);
    }

    template<class T, class U>
    [[nodiscard]] constexpr bool cmp_less_equal(T t, U u) noexcept {
        return !cmp_greater(t, u);
    }

    template<class T, class U>
    [[nodiscard]] constexpr bool cmp_greater_equal(T t, U u) noexcept {
        return !cmp_less(t, u);
    }

    /*! \brief MatLab-style modulo operator
     *   \param x first number
     *   \param y second number
     *   \return modulo of x and y. Example, <code> mod(7,2)  = 1 </code> but <code> mod(-0.5,10)  = 9.5 </code>, instead of <code> -0.5 </code>  as given by
     * x%y.
     */
    template<typename T>
    [[nodiscard]] inline T mod(const T x, const T y) {
        if constexpr(!ndebug)
            if(y == 0) throw("num::mod(x,y): divisor y == 0");
        if constexpr(std::is_integral_v<T>) {
            if constexpr(std::is_unsigned_v<T>)
                return x >= y ? x % y : x;
            else {
                return x >= y ? x % y : (x < 0 ? (x % y + y) % y : x);
            }

        }
        //            return (x % y + y) % y;
        else
            return std::fmod((std::fmod(x, y) + y), y);
    }

    /*! \brief Similar to mod but faster for use with periodic boundary condition
     *   \param x first number
     *   \param y second number
     *   \return modulo of x and y. Example, <code> mod(7,2)  = 1 </code> but <code> mod(-0.5,10)  = 9.5 </code>, instead of <code> -0.5 </code>  as given by
     * x%y.
     */
    template<typename T>
    [[nodiscard]] inline T pbc(const T x, const T y) {
        if constexpr(!ndebug)
            if(y == 0) throw("num::pbc(x,y): divisor y == 0");
        if constexpr(std::is_signed_v<T>) {
            if(x >= 0 and x < y) return x;
            if(x < 0 and x >= -2 * y) return x + y;
        } else {
            if(x < y) return x;
        }
        if(x >= y and x < 2 * y) return x - y;
        return num::mod(x, y);
    }

    template<typename T>
    [[nodiscard]] int sign(const T val) noexcept {
        if(val > 0) return +1;
        if(val < 0) return -1;
        return 0;
    }

    template<typename T>
    [[nodiscard]] bool between(const T &value, const T &low, const T &high) noexcept {
        return value >= low and value <= high;
    }

    /*! \brief Python-style range generator, i.e. not-including "last"
     *   \return Range of T's. Example, <code> range(0,8,2) </code> gives a std::vector<int>: <code> [0,2,4,6] </code>
     */
    namespace internal {
        template<typename TA, typename TB>
        using int_or_dbl = typename std::conditional<std::is_floating_point_v<TA> or std::is_floating_point_v<TB>, double, int>::type;
    }

    template<typename T = int, typename T1, typename T2, typename T3 = internal::int_or_dbl<T1, T2>>
    [[nodiscard]] std::vector<T> range(T1 first, T2 last, T3 step = static_cast<T3>(1)) {
        if(step == 0) throw std::runtime_error("Range cannot have step size zero");
        if constexpr(std::is_signed_v<T3>) {
            if(cmp_greater(first, last) and step > 0) return range<T>(first, last, -step);
            if(cmp_less(first, last) and step < 0) return range<T>(first, last, -step);
        } else {
            if(cmp_greater(first, last)) throw std::runtime_error("Range of unsigned step type cannot have first > last");
        }
        if(cmp_equal(first, last)) return {};
        if(cmp_equal(static_cast<T3>(first) + step, last)) return std::vector<T>{static_cast<T>(first)};

        auto num_steps = static_cast<size_t>(
            std::abs<double>((static_cast<double>(last) - static_cast<double>(first) + static_cast<double>(step)) / static_cast<double>(step)));
        if(num_steps > 10000000) throw std::runtime_error("Too many steps for range");

        std::vector<T> vec;
        vec.reserve(num_steps);
        for(T3 current = static_cast<T3>(first); cmp_less(current, last); current += step) vec.emplace_back(current);
        return vec;
    }

    /*! \brief MatLab-style linearly spaced array
     *   \param num number of linearly spaced values
     *   \param a first value in range
     *   \param b last value in range
     *   \return std::vector<T2>. Example,  <code> Linspaced(5,1,5) </code> gives a std::vector<int>: <code> [1,2,3,4,5] </code>
     */
    [[nodiscard]] inline std::vector<double> LinSpaced(std::size_t N, double a, double b) {
        double              h = (b - a) / static_cast<double>(N - 1);
        std::vector<double> xs(N);
        double              val = a;
        for(auto &x : xs) {
            x = val;
            val += h;
        }
        return xs;
    }

    [[nodiscard]] inline std::vector<double> LogSpaced(std::size_t N, double a, double b, double base = 10.0) {
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
    [[nodiscard]] auto prod(const Input &in, const From from, const To to) {
        return std::accumulate(in.data() + from, in.data() + to, 1, std::multiplies<>());
    }

    /*! \brief Checks if multiple values are equal to each other
     *   \param args any number of values
     *   \return bool, true if all args are equal
     */
    template<typename First, typename... T>
    [[nodiscard]] bool all_equal(First &&first, T &&...t) noexcept {
        return ((first == t) && ...);
    }

    template<typename R, typename T>
    [[nodiscard]] R next_power_of_two(T val) {
        return static_cast<R>(std::pow<long>(2, static_cast<long>(std::ceil(std::log2(std::real(val))))));
    }
    template<typename R, typename T>
    [[nodiscard]] R prev_power_of_two(T val) {
        return static_cast<R>(std::pow<long>(2, static_cast<long>(std::floor(std::log2(std::real(val - 1))))));
    }

    template<typename R, typename T>
    [[nodiscard]] inline R next_multiple(const T num, const T mult) {
        if(mult == 0) return num;
        return (num + mult) - mod(num, mult);
    }
    template<typename R, typename T>
    [[nodiscard]] inline R prev_multiple(const T num, const T mult) {
        if(mult == 0) return num;
        auto m = mod(num, mult);
        if(m == 0) return prev_multiple<R>(num - 1, mult);
        return num - m;
    }

    template<typename T>
    [[nodiscard]] inline T round_to_multiple_of(const T number, const T multiple) {
        T result = number + multiple / 2;
        result -= num::mod(result, multiple);
        return result;
    }

}
