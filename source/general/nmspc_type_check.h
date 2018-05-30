//
// Created by david on 2018-02-06.
//

#ifndef NMSPC_TYPE_CHECK_H
#define NMSPC_TYPE_CHECK_H
#include <experimental/type_traits>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

namespace TypeCheck{
    template <typename T> using Data_t          = decltype(std::declval<T>().template data());
    template <typename T> using Size_t          = decltype(std::declval<T>().template size());
    template <typename T> using Cstr_t          = decltype(std::declval<T>().template c_str());
    template <typename T> using Imag_t          = decltype(std::declval<T>().template imag());
    template <typename T> using Scal_t          = typename T::Scalar;
    template <typename T> using Valt_t          = typename T::value_type;

    template <typename T> using has_member_data         = std::experimental::is_detected<Data_t, T>;
    template <typename T> using has_member_size         = std::experimental::is_detected<Size_t, T>;
    template <typename T> using has_member_scalar       = std::experimental::is_detected<Scal_t , T>;
    template <typename T> using has_member_value_type   = std::experimental::is_detected<Valt_t , T>;
    template <typename T> using has_member_c_str   = std::experimental::is_detected<Cstr_t , T>;
    template <typename T> using has_member_imag    = std::experimental::is_detected<Imag_t , T>;


    template<typename T> struct is_vector : public std::false_type {};
    template<typename T> struct is_vector<std::vector<T>> : public std::true_type {};


    template<typename T> struct is_eigen_tensor : public std::false_type {};
    template<typename T, size_t N, int StorageOrder>
    struct is_eigen_tensor<Eigen::Tensor<T, N, StorageOrder>> : public std::true_type{};

    template<typename T> struct is_eigen_matrix : public std::false_type {};
    template<typename T, int rows, int cols, int StorageOrder> struct is_eigen_matrix<Eigen::Matrix<T,rows,cols,StorageOrder>> : public std::true_type {};

    template<typename T> struct is_eigen_array : public std::false_type {};
    template<typename T, int rows, int cols, int StorageOrder> struct is_eigen_array<Eigen::Array<T,rows,cols,StorageOrder>> : public std::true_type {};

    template<typename T>
    constexpr bool is_eigen_matrix_or_array(){
        if constexpr(is_eigen_matrix<T>::value or
                     is_eigen_array<T>::value ){
            return true;
        }else{
            return false;
        }
    }
    //This does not work for "non-type" class template parameters.
    //In fact it doesn't seem to work very well at all...
//    template < template <typename...> class Template, typename T >
//    struct is_instance_of : std::false_type {};
//
//    template < template <typename...> class Template, typename... Args >
//    struct is_instance_of< Template, Template<Args...> > : std::true_type {};
//
//    template <typename T> using is_ofEigen              = is_instance_of<Eigen::EigenBase,T>;
}

#endif //PT_NMSPC_TYPE_CHECK_H

