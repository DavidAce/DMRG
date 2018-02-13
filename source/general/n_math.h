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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <numeric>
#include <Eigen/Core>



/*! \brief Prints the content of a vector nicely */
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    if (!v.empty()) {
        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}



/*!
 *  \namespace Math
 *  \brief Small convenience-type math functions like modulo and numerical integration using GSL.
 *  \tableofcontents

 */



namespace Math {
    /*! \fn auto mod(const T1 x, const T2 y)
    *  \brief MatLab-style modulo operator
    *  \param x first number
    *  \param y second number
    *  \return modulo of x and y. Example,  <code> mod(-0.5,10)  = 9.5 </code>, instead of <code> -0.5 </code>  as given by x%y.
    */
    template<typename T1, typename T2>
    inline auto mod(const T1 x, const T2 y) {
        return x >= 0 ? x % y : x % y + y;
    }

    /*! \fn std::vector<T2> LinSpaced(T1 num, T2 min, T2 max )
    *  \brief MatLab-style linearly spaced array
    *  \param num number of linearly spaced values
    *  \param min minimum value in range
    *  \param max maximum value in range
    *  \return std::vector<T2>. Example,  <code> Linspaced(5,1,5) </code> gives a std::vector<int>: <code> [1,2,3,4,5] </code>
    */
    template <typename T1, typename T2>
    inline std::vector<T2> LinSpaced(T1 num, T2 min, T2 max ){
        static_assert(std::is_integral<T1>::value && "Math::LinSpaced -- Given type is not integral!");
        Eigen::Array<T2, Eigen::Dynamic, 1> temp =  Eigen::Array<T2, Eigen::Dynamic, 1> :: LinSpaced(num, min, max);
        return std::vector<T2> (temp.data(), temp.data() + temp.size());
    }


    /*! \fn auto prod(const Input &in, const From from, const To to)
    *   \brief Product operator for containers such as vector
    *   \param in a vector, array or any 1D container with "<code> .data() </code>" method.
    *   \param from first element to multiply
    *   \param to last element to multiply
    *   \return std::vector<T2>. Example, let <code> my_vector = {1,2,3,4}</code>. Then <code> prod(my_vector,0,3) = 24 </code>.
    */
    template<typename Input, typename From, typename To>
    auto prod(const Input &in, const From from, const To to){
        return std::accumulate(in.data() + from, in.data()+to,1,std::multiplies<>());
    };


    /*! \brief "pow", x^p for integers x and p using recursion */
    template<typename IntegerType, typename = typename std::enable_if<std::is_integral<IntegerType>::value>::type>
    IntegerType ipow(IntegerType x, IntegerType p) {
        if (p == 0) return 1;
        if (p == 1) return x;
        int tmp = ipow(x, p / 2);
        if (p % 2 == 0) return tmp * tmp;
        else return x * tmp * tmp;
    }




    /*! \brief class for integration */
    template<typename F>
    class gsl_quad {
        F f;
        int limit;
        std::unique_ptr<gsl_integration_workspace,
                std::function<void(gsl_integration_workspace *)>
        > workspace;

        static double gsl_wrapper(double x, void *p) {
            gsl_quad *t = reinterpret_cast<gsl_quad *>(p);
            return t->f(x);
        }

    public:
        gsl_quad(F f, int limit)
                : f(f), limit(limit), workspace(gsl_integration_workspace_alloc(limit), gsl_integration_workspace_free) {}

        double integrate(double min, double max, double epsabs, double epsrel) {
            gsl_function gsl_f;
            gsl_f.function = &gsl_wrapper;
            gsl_f.params = this;

            double result, error;
            if (!std::isinf(min) && !std::isinf(max)) {
                gsl_integration_qags(&gsl_f, min, max,
                                     epsabs, epsrel, limit,
                                     workspace.get(), &result, &error);
            } else if (std::isinf(min) && !std::isinf(max)) {
                gsl_integration_qagil(&gsl_f, max,
                                      epsabs, epsrel, limit,
                                      workspace.get(), &result, &error);
            } else if (!std::isinf(min) && std::isinf(max)) {
                gsl_integration_qagiu(&gsl_f, min,
                                      epsabs, epsrel, limit,
                                      workspace.get(), &result, &error);
            } else {
                gsl_integration_qagi(&gsl_f,
                                     epsabs, epsrel, limit,
                                     workspace.get(), &result, &error);
            }

            return result;
        }
    };

    /*! \fn compute_integral(F func,
     *                    std::pair<double, double> const &range,
     *                   double epsabs = 1.49e-8, double epsrel = 1.49e-8,
     *                   int limit = 50)
     * \brief Numerical integrations.
     *
     * Examples of integration:
     *
     * \f$ \int_0^1 x^2 dx = \f$
     * ```{cpp}
     *      compute_integral([](double x) { return x*x; }, {0,1})
     * ```
     *
     * \f$ \int_1^\infty x^{-2} dx = \f$
     * ```{cpp}
     *      compute_integral([](double x) { return 1/(x*x); }, {1,INFINITY})
     * ```
     *
     *  \f$ \int_{-\infty}^\infty \exp(-x^2) dx = \f$
     * ```{cpp}
     * compute_integral([](double x) { return std::exp(-x*x); }, {-INFINITY,INFINITY})
     * ```
     *
     */
    template<typename F>
    double compute_integral(F func,
                            std::pair<double, double> const &range,
                            double epsabs = 1.49e-8, double epsrel = 1.49e-8,
                            int limit = 50) {
        return gsl_quad<F>(func, limit).integrate(range.first, range.second, epsabs, epsrel);


    }



}



#endif //TRAINING_FUNCS_H
