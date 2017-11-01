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
#include <n_settings.h>
#include <n_tensor_extra.h>
//#include <GenEigsSolver.h>
//#include <SymGEigsSolver.h>
#include <SymEigsSolver.h>
#include <MatOp/DenseSymMatProd.h>
#include <GenEigsSolver.h>



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


/*! \brief Math functions.*/

/*!
 *  \namespace Math
 *  Contains functions do math on more general types,
 *  print out vectors and do numerical integration on lambdas using GSL.
 *
 * ## Examples of integration:
 *
 * \f$ \int_0^1 x^2 dx = \f$
 * ```{cpp}
 *      compute_integral([](double x) { return x*x; }, {0,1})
 * ```
 * \f$ \int_1^\infty x^{-2} dx = \f$
 * ```{cpp}
 *      compute_integral([](double x) { return 1/(x*x); }, {1,INFINITY})
 * ```
 *  \f$ \int_{-\infty}^\infty \exp(-x^2) dx = \f$
 * ```{cpp}
 * compute_integral([](double x) { return std::exp(-x*x); }, {-INFINITY,INFINITY})
 * ```
 */



namespace Math {
    /*! \brief       MatLab-style modulo operator */
    template<typename T1, typename T2>
    inline auto mod(const T1 x, const T2 y) {
        return x >= 0 ? x % y : x % y + y;
    }

    template <typename T1, typename T2>
    inline std::vector<T2> LinSpaced(T1 num, T2 min, T2 max ){
        static_assert(std::is_integral<T1>::value && "Math::LinSpaced -- Given type is not integral!");
        Eigen::Array<T2, Eigen::Dynamic, 1> temp =  Eigen::Array<T2, Eigen::Dynamic, 1> :: LinSpaced(num, min, max);
        return std::vector<T2> (temp.data(), temp.data() + temp.size());
    }


    /*! \brief  Product operator for containers such as vector */
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



    enum class Handedness {RGHT, LEFT};

template <typename Scalar, Handedness hand = Handedness::RGHT, int Uplo = Eigen::Lower>
class DenseSymMatProd{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef Eigen::Map<const Matrix> MapConstMat;
        typedef Eigen::Map<const Vector> MapConstVec;
        typedef Eigen::Map<Vector> MapVec;

        typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;

        const MapConstMat m_mat;

    public:
        DenseSymMatProd(ConstGenericMatrix& mat_) :
            m_mat(mat_.data(), mat_.rows(), mat_.cols())
        {}

    int rows() const { return m_mat.rows(); }
    int cols() const { return m_mat.cols(); }

    // y_out = A * x_in
    void perform_op(const Scalar* x_in, Scalar* y_out) const
    {
        MapConstVec x(x_in,  m_mat.cols());
        MapVec      y(y_out, m_mat.rows());
        y.noalias() = m_mat.template selfadjointView<Uplo>() * x;
        switch (hand){
            case Handedness::RGHT:
                y.noalias() = m_mat.template selfadjointView<Uplo>() * x;
                break;
            case Handedness::LEFT:
                y.noalias() = x.transpose() * m_mat.template selfadjointView<Uplo>();
                break;
        }

    }
};


    template <typename Scalar, Handedness hand = Handedness::RGHT>
    class DenseGenMatProd
    {
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef Eigen::Map<const Matrix> MapConstMat;
        typedef Eigen::Map<const Vector> MapConstVec;
        typedef Eigen::Map<Vector> MapVec;
        typedef const Eigen::Ref<const Matrix> ConstGenericMatrix;
        const MapConstMat m_mat;

    public:
        DenseGenMatProd(ConstGenericMatrix& mat_) :
                m_mat(mat_.data(), mat_.rows(), mat_.cols())
        {}

        int rows() const { return m_mat.rows(); }
        int cols() const { return m_mat.cols(); }

        // y_out = A * x_in
        void perform_op(const Scalar* x_in, Scalar* y_out) const
        {
            MapConstVec x(x_in,  m_mat.cols());
            MapVec      y(y_out, m_mat.rows());
            switch (hand){
                case Handedness::RGHT:
                    y.noalias() = m_mat * x;
                    break;
                case Handedness::LEFT:
                    y.noalias() = x.transpose() * m_mat;
                    break;
            }
        }
    };


    /*! \brief Returns the largest eigenvector of a matrix */
    template<Handedness hand = Handedness::RGHT>
    Textra::Tensor<1,Textra::cdouble> dominant_eigenvector(const Textra::Tensor<2, double> &tensor){
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
        DenseGenMatProd<double, hand> op(mat);
        int ncv = std::min(settings::precision::eig_max_ncv, 3);
        Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, DenseGenMatProd<double,hand>> eigs(&op, 1, ncv);
        eigs.init();
        eigs.compute(10000, 1e-12, Spectra::LARGEST_REAL);
        if(eigs.info() != Spectra::SUCCESSFUL){
            std::cout << "Eigenvalue solver failed." << '\n';
//            std::cout << tensor << '\n';
//            exit(1);
        }
        std::cout << "Eigenvalue: " << eigs.eigenvalues()<< std::endl;
//        double x = eigs.eigenvectors().real().coeff(0);
//        double sign = (x < 0) ? -1.0/x : 1.0/x;
//        double sign = 1.0/x;
        return Textra::Matrix_to_Tensor1<Textra::cdouble>(eigs.eigenvectors());
    }

    /*! \brief Returns the largest eigenvector of a matrix */
//    template<Handedness hand = Handedness::RGHT>
//    auto dominant_eigenvector_eigenvalue(const Textra::Tensor<2, double> &tensor, const Textra::array2 out_shape){
//        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
//        DenseGenMatProd<double, hand> op(mat);
//        int ncv = std::min(settings::precision::eig_max_ncv, 3);
//        Spectra::GenEigsSolver<double, Spectra::LARGEST_REAL, DenseGenMatProd<double,hand>> eigs(&op, 1, ncv);
//        eigs.init();
//        eigs.compute(10000, 1e-12, Spectra::LARGEST_REAL);
//        if(eigs.info() != Spectra::SUCCESSFUL){
//            std::cout << "Eigenvalue solver failed." << '\n';
////            std::cout << tensor << '\n';
////            exit(1);
//        }
//        std::cout << "Eigenvalue: " << eigs.eigenvalues()<< std::endl;
////        double x = eigs.eigenvectors().real().coeff(0);
////        double sign = (x < 0) ? -1.0/x : 1.0/x;
////        double sign = 1.0/x;
////        return {Textra::Matrix_to_Tensor}
//        return std{Textra::Matrix_to_Tensor1<Textra::cdouble>(eigs.eigenvectors()), eigs.eigenvalues()};
//    }

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

    /*!
    * \fn compute_integral(F func,
                        std::pair<double, double> const &range,
                        double epsabs = 1.49e-8, double epsrel = 1.49e-8,
                        int limit = 50)
    * \brief Function to compute integrals.

    * */
    template<typename F>
    double compute_integral(F func,
                            std::pair<double, double> const &range,
                            double epsabs = 1.49e-8, double epsrel = 1.49e-8,
                            int limit = 50) {
        return gsl_quad<F>(func, limit).integrate(range.first, range.second, epsabs, epsrel);


    }



}



#endif //TRAINING_FUNCS_H
