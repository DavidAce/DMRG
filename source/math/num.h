#pragma once

#include <cmath>
#include <complex>
#include <iterator>
#include <numeric>
#include <optional>
#include <vector>

/*!
 *  \namespace num
 *  \brief Small convenience-type num functions like modulo
 *  \tableofcontents
 */

namespace num {

    template<typename ContainerType>
    void check_bounds(ContainerType &X, std::optional<long> start_point = std::nullopt, std::optional<long> end_point = std::nullopt) {
        if(start_point.has_value() and (start_point.value() >= static_cast<long>(X.size()) or start_point.value() < 0))
            throw std::range_error("Start point is out of range");
        if(end_point.has_value() and (end_point.value() > static_cast<long>(X.size()) or end_point.value() < start_point))
            throw std::range_error("End point is out of range");
    }

    template<typename ContainerType>
    double mean(ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        try {
            check_bounds(X, start_point, end_point);
        } catch(std::exception &err) { throw std::range_error("mean: " + std::string(err.what())); }
        if(not start_point.has_value()) start_point = 0;
        if(not end_point.has_value()) end_point = X.size();
        if(end_point.value() == start_point.value()) return 0.0;
        if(end_point.value() < start_point.value()) throw std::runtime_error("end_point < start_point");

        auto x_it_start = X.begin();
        auto x_it_end   = X.begin();
        auto n          = static_cast<double>(end_point.value() - start_point.value());
        std::advance(x_it_start, start_point.value());
        std::advance(x_it_end, end_point.value());
        return accumulate(x_it_start, x_it_end, 0.0) / n;
    }

    template<typename ContainerType>
    double stdev(const ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        try {
            check_bounds(X, start_point, end_point);
        } catch(std::exception &err) { throw std::range_error("stdev: " + std::string(err.what())); }
        if(not start_point.has_value()) start_point = 0;
        if(not end_point.has_value()) end_point = X.size();
        if(end_point.value() == start_point.value()) return 0.0;
        if(end_point.value() < start_point.value()) throw std::runtime_error("end_point < start_point");

        double X_mean = num::mean(X, start_point, end_point);
        auto   n      = static_cast<double>(end_point.value() - start_point.value());
        auto   x_it   = X.begin();
        auto   x_en   = X.begin();

        std::advance(x_it, start_point.value());
        std::advance(x_en, end_point.value());

        double sum = std::accumulate(x_it, x_en, 0.0, [&X_mean](auto &x1, auto &x2) { return x1 + (x2 - X_mean) * (x2 - X_mean); });

        return std::sqrt(sum / n);
    }

    template<typename ContainerType>
    double sterr(const ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        return stdev(X, start_point, end_point) / std::sqrt(X.size());
    }

    template<typename ContainerType1, typename ContainerType2>
    double slope(ContainerType1 &X, ContainerType2 &Y, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        if(X.size() != Y.size()) throw std::range_error("slope: Size mismatch in arrays");
        try {
            check_bounds(X, start_point, end_point);
        } catch(std::exception &err) { throw std::range_error("slope: " + std::string(err.what())); }
        if(not start_point.has_value()) start_point = 0;
        if(not end_point.has_value()) end_point = X.size();
        if(end_point.value() == start_point.value()) return 0.0;
        if(end_point.value() < start_point.value()) throw std::runtime_error("end_point < start_point");

        auto x_it = X.begin();
        auto x_en = X.begin();
        auto y_it = Y.begin();
        auto y_en = Y.begin();
        std::advance(x_it, start_point.value());
        std::advance(y_it, start_point.value());
        std::advance(x_en, end_point.value());
        std::advance(y_en, end_point.value());
        auto   n    = static_cast<double>(end_point.value() - start_point.value());
        double avgX = accumulate(x_it, x_en, 0.0) / n;
        double avgY = accumulate(y_it, y_en, 0.0) / n;

        double numerator   = 0.0;
        double denominator = 0.0;
        while(x_it != x_en) {
            numerator += (static_cast<double>(*x_it) - avgX) * (static_cast<double>(*y_it) - avgY);
            denominator += (static_cast<double>(*x_it) - avgX) * (static_cast<double>(*x_it) - avgX);
            y_it++;
            x_it++;
        }
        return std::abs(numerator / denominator);
    }

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
    template<typename T>
    std::vector<T> range(T first, T last, T step = static_cast<T>(1)) {
        if(step == 0) throw std::runtime_error("Range cannot have step size zero");
        if constexpr(std::is_signed_v<T>){
            if(first > last and step > 0) return range(first, last, -step);
            if(first < last and step < 0) return range(first, last, -step);
        }else{
            if(first > last) throw std::runtime_error("Range of unsigned type cannot have first > last");
        }
        if(first == last) return {};
        if(first + step == last) return {first};

        auto num_steps = static_cast<size_t>(
            std::abs<double>((static_cast<double>(last) - static_cast<double>(first) + static_cast<double>(step)) / static_cast<double>(step)));
        if(num_steps > 10000000) throw std::runtime_error("Too many steps for range");

        std::vector<T> vec;
        vec.reserve(num_steps);
        for(T current = first; current < last; current += step) vec.emplace_back(static_cast<T>(current));
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

    //    inline std::vector<double> LogSpaced(std::size_t N, double a, double b, double base = 10.0) {
    //        double expA = std::pow(base, a);
    //        double expB = std::pow(base, b);
    //        auto   vec  = LinSpaced(N, expA, expB);
    //        for(auto &val : vec) val = std::abs(std::log(val) / std::log(base));
    //        return vec;
    //    }

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
