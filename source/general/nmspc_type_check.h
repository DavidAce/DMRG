//
// Created by david on 2018-02-06.
//

#ifndef NMSPC_TYPE_CHECK_H
#define NMSPC_TYPE_CHECK_H

#include <Eigen/Core>

namespace TypeCheck{
    template <typename T> using Data_t          = decltype(std::declval<T>().template data());
    template <typename T> using Size_t          = decltype(std::declval<T>().template size());
    template <typename T> using Dims_t          = decltype(std::declval<T>().template dimensions());
    template <typename T> using Rank_t          = decltype(std::declval<T>().template rank());
    template <typename T> using Setz_t          = decltype(std::declval<T>().template setZero());
    template <typename T> using Matr_t          = decltype(std::declval<T>().template matrix());
    template <typename T> using Scal_t          = typename T::Scalar;
    template <typename T> using Valt_t          = typename T::value_type;

    template <typename T> using has_member_data         = std::experimental::is_detected<Data_t, T>;
    template <typename T> using has_member_size         = std::experimental::is_detected<Size_t, T>;
    template <typename T> using has_member_dimensions   = std::experimental::is_detected<Dims_t , T>;
    template <typename T> using has_member_scalar       = std::experimental::is_detected<Scal_t , T>;
    template <typename T> using has_member_value_type   = std::experimental::is_detected<Valt_t , T>;
    template <typename T> using is_eigen                = std::experimental::is_detected<Setz_t , T>;
    template <typename T> using is_matrix               = std::experimental::is_detected<Matr_t , T>;
    template <typename T> using is_tensor               = std::experimental::is_detected<Rank_t , T>;

    //This does not work for "non-type" class template parameters.
    //In fact it doesn't seem to work very well at all...
    template < template <typename...> class Template, typename T >
    struct is_instance_of : std::false_type {};

    template < template <typename...> class Template, typename... Args >
    struct is_instance_of< Template, Template<Args...> > : std::true_type {};

    template <typename T> using is_ofEigen              = is_instance_of<Eigen::EigenBase,T>;
}

#endif //PT_NMSPC_TYPE_CHECK_H

