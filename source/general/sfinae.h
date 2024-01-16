#pragma once
#include <array>
#include <complex>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

/*!
 * \brief A collection of type-detection and type-analysis utilities using SFINAE
 */
namespace sfinae {

    // SFINAE detection
    // helper constant for static asserts
    template<class>
    inline constexpr bool always_false_v = false;
    template<class>
    inline constexpr bool unrecognized_type_v = false;
    template<class>
    inline constexpr bool invalid_type_v = false;

    template<typename...>
    struct print_type_and_exit_compile_time;

    template<typename T>
    constexpr auto type_name() {
        std::string_view name, prefix, suffix;
#ifdef __clang__
        name   = __PRETTY_FUNCTION__;
        prefix = "auto sfinae::type_name() [T = ";
        suffix = "]";
#elif defined(__GNUC__)
        name   = __PRETTY_FUNCTION__;
        prefix = "constexpr auto sfinae::type_name() [with T = ";
        suffix = "]";
#elif defined(_MSC_VER)
        name   = __FUNCSIG__;
        prefix = "auto __cdecl sfinae::type_name<";
        suffix = ">(void)";
#endif
        name.remove_prefix(prefix.size());
        name.remove_suffix(suffix.size());
        return name;
    }

#if __cplusplus >= 202002L
    template<typename T1, typename T2>
    concept type_is = std::same_as<std::remove_cvref_t<T1>, T2>;

    template<typename T, typename... Ts>
    concept is_any_v = (type_is<T, Ts> || ...);

    template<typename T, typename... Ts>
    concept are_same_v = (type_is<T, Ts> && ...);

    template<typename T>
    concept is_pointer_type = std::is_pointer_v<T>;

    template<typename T>
    concept is_arithmetic_type = std::is_arithmetic_v<T>;

    template<typename T>
    concept has_data_v = requires(T m) {
        { m.data() } -> is_pointer_type;
    };

    template<typename T>
    concept has_size_v = requires(T m) {
        { m.size() } -> std::integral;
    };

    template<typename T>
    concept has_resize_v = requires(T m) {
        { m.resize() } -> std::same_as<void>;
    };

    template<typename T>
    concept has_resize2_v = requires(T m) {
        { m.resize(0, 0) } -> std::same_as<void>;
    };

    template<typename T>
    concept has_value_type_v = requires(T m) {
        { T::value_type };
    };

    template<typename T>
    concept has_c_str_v = requires(T m) {
        { m.c_str() };
    };

    template<typename T>
    concept has_imag_v = requires(T m) {
        { m.imag() } -> is_arithmetic_type;
    };

    template<typename T>
    concept has_Scalar_v = requires(T m) {
        { T::Scalar };
    };

    template<typename T>
    concept has_NumIndices_v = requires(T m) {
        { T::NumIndices };
    };

    template<typename T>
    concept has_dimensions_v = requires(T m) {
        { m.dimensions };
    };

    template<typename T>
    concept has_x_v = requires(T m) {
        { m.x };
    };
    template<typename T>
    concept has_y_v = requires(T m) {
        { m.y };
    };
    template<typename T>
    concept has_z_v = requires(T m) {
        { m.z };
    };

    template<typename T>
    concept is_std_optional_v = type_is<T, std::optional<typename T::value_type>>;

    template<typename T>
    concept is_vector_v = type_is<T, std::vector<typename T::value_type>>;

    template<typename T>
    concept is_array_v = type_is<T, typename std::array<typename T::value_type, std::tuple_size<T>::value>>;

    template<typename T>
    concept is_std_complex_v = type_is<T, std::complex<typename T::value_type>>;

    template<typename T, typename T2>
    concept is_pair_v = type_is<T, std::pair<typename T::first_type, typename T::second_type>>;

    template<typename Test, template<typename...> class Ref>
    struct is_specialization : std::false_type {};
    template<template<typename...> class Ref, typename... Args>
    struct is_specialization<Ref<Args...>, Ref> : std::true_type {};
    template<typename T, template<typename...> class Ref>
    inline constexpr bool is_specialization_v = is_specialization<T, Ref>::value;

    template<typename T>
    concept is_iterable_v = requires(T m) {
        { m.begin() };
        { m.end() };
    };

    template<typename T>
    concept is_text_v = has_c_str_v<T> || std::convertible_to<T, std::string_view> || std::is_same_v<T, char>;
    //    static_assert(is_text_v<std::string>);
    //    static_assert(is_text_v<std::string_view>);
    //    static_assert(is_text_v<char *>);
    //    static_assert(is_text_v<char []>);
    //    static_assert(is_text_v<char>);

    template<typename T>
    concept has_text_v = !is_text_v<T> && (is_text_v<typename std::remove_all_extents_t<T>> || is_text_v<typename std::remove_pointer_t<T>> ||
                                           is_text_v<typename T::value_type>);
    //    static_assert(has_text_v<std::vector<std::string>>);
    //    static_assert(has_text_v<std::array<std::string_view, 0>>);
    //    static_assert(has_text_v<std::string_view[]>);

    //    template<typename T>
    //    concept is_Scalar2_v = has_x_v<T> && has_y_v<T> && (sizeof(T) == sizeof(T::x) + sizeof(T::y));

    template<typename T>
    concept is_Scalar2_v = requires {
        { T::x };
        { T::y };
        { std::is_same_v<typename T::x, typename T::y> };
        { sizeof(T) == sizeof(T::x) + sizeof(T::y) };
    };
    template<typename T>
    concept is_Scalar3_v = requires {
        { T::x };
        { T::y };
        { T::z };
        { std::is_same_v<typename T::x, typename T::y> };
        { std::is_same_v<typename T::x, typename T::z> };
        { sizeof(T) == sizeof(T::x) + sizeof(T::y) + sizeof(T::z) };
    };

    template<typename O, typename I>
    concept is_Scalar2_of_type = is_Scalar2_v<O> && std::is_same_v<I, decltype(O::x)>;

    template<typename O, typename I>
    concept is_Scalar3_of_type = is_Scalar3_v<O> && std::is_same_v<I, decltype(O::x)>;

    template<typename T>
    concept is_ScalarN = is_Scalar2_v<T> || is_Scalar3_v<T>;

#else
    template<class T, class... Ts>
    struct is_any : std::disjunction<std::is_same<T, Ts>...> {};
    template<class T, class... Ts>
    inline constexpr bool is_any_v = is_any<T, Ts...>::value;

    template<class T, class... Ts>
    struct are_same : std::conjunction<std::is_same<T, Ts>...> {};
    template<class T, class... Ts>
    inline constexpr bool are_same_v = are_same<T, Ts...>::value;

    template<typename T, typename = std::void_t<>>
    struct has_data : public std::false_type {};
    template<typename T>
    struct has_data<T, std::void_t<decltype(std::declval<T>().data())>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_data_v = has_data<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_size : public std::false_type {};
    template<typename T>
    struct has_size<T, std::void_t<decltype(std::declval<T>().size())>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_size_v = has_size<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_resize : public std::false_type {};
    template<typename T>
    struct has_resize<T, std::void_t<decltype(std::declval<T>().resize(0))>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_resize_v = has_resize<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_resize2 : public std::false_type {};
    template<typename T>
    struct has_resize2<T, std::void_t<decltype(std::declval<T>().resize(0, 0))>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_resize2_v = has_resize2<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_value_type : public std::false_type {};
    template<typename T>
    struct has_value_type<T, std::void_t<typename T::value_type>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_value_type_v = has_value_type<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_c_str : public std::false_type {};
    template<typename T>
    struct has_c_str<T, std::void_t<decltype(std::declval<T>().c_str())>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_c_str_v = has_c_str<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_imag : public std::false_type {};
    template<typename T>
    struct has_imag<T, std::void_t<decltype(std::declval<T>().imag())>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_imag_v = has_imag<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_Scalar : public std::false_type {};
    template<typename T>
    struct has_Scalar<T, std::void_t<typename T::Scalar>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_Scalar_v = has_Scalar<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_NumIndices : public std::false_type {};
    template<typename T>
    struct has_NumIndices<T, std::void_t<decltype(std::declval<T>().NumIndices)>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_NumIndices_v = has_NumIndices<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_dimensions : public std::false_type {};
    template<typename T>
    struct has_dimensions<T, std::void_t<decltype(std::declval<T>().dimensions())>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_dimensions_v = has_dimensions<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_x : public std::false_type {};
    template<typename T>
    struct has_x<T, std::void_t<decltype(std::declval<T>().x)>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_x_v = has_x<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_y : public std::false_type {};
    template<typename T>
    struct has_y<T, std::void_t<decltype(std::declval<T>().y)>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_y_v = has_y<T>::value;

    template<typename T, typename = std::void_t<>>
    struct has_z : public std::false_type {};
    template<typename T>
    struct has_z<T, std::void_t<decltype(std::declval<T>().z)>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_z_v = has_z<T>::value;

    template<typename T>
    struct is_std_optional : public std::false_type {};
    template<typename T>
    struct is_std_optional<std::optional<T>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_std_optional_v = is_std_optional<T>::value;

    template<typename T>
    struct is_std_vector : public std::false_type {};
    template<typename T>
    struct is_std_vector<std::vector<T>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_std_vector_v = is_std_vector<T>::value;

    template<typename T>
    struct is_std_array : public std::false_type {};
    template<typename T, auto N>
    struct is_std_array<std::array<T, N>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_std_array_v = is_std_array<T>::value;

    template<typename T>
    struct is_std_complex : public std::false_type {};
    template<typename T>
    struct is_std_complex<std::complex<T>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_std_complex_v = is_std_complex<T>::value;

    template<typename T>
    struct is_pair : std::false_type {};
    template<typename T, typename U>
    struct is_pair<std::pair<T, U>> : std::true_type {};
    template<typename T>
    inline constexpr bool is_pair_v = is_pair<T>::value;

    template<typename Test, template<typename...> class Ref>
    struct is_specialization : std::false_type {};
    template<template<typename...> class Ref, typename... Args>
    struct is_specialization<Ref<Args...>, Ref> : std::true_type {};
    template<typename T, template<typename...> class Ref>
    inline constexpr bool is_specialization_v = is_specialization<T, Ref>::value;

    template<typename T, typename = std::void_t<>>
    struct is_iterable : public std::false_type {};
    template<typename T>
    struct is_iterable<T, std::void_t<decltype(std::begin(std::declval<T>())), decltype(std::end(std::declval<T>()))>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_iterable_v = is_iterable<T>::value;

    template<typename T>
    struct is_text {
        private:
        template<typename U>
        static constexpr bool test() {
            using DecayType = typename std::decay<U>::type;
            // No support for wchar_t, char16_t and char32_t
            if constexpr(has_c_str_v<DecayType>) return true;
            if constexpr(std::is_convertible_v<DecayType, std::string_view>) return true;
            if constexpr(std::is_same_v<DecayType, char>)
                return true;
            else
                return false;
        }

        public:
        static constexpr bool value = test<T>();
    };
    template<typename T>
    inline constexpr bool is_text_v = is_text<T>::value;

    template<typename T>
    struct has_text {
        private:
        template<typename U>
        static constexpr bool test() {
            using DecayType = typename std::decay<U>::type;
            if constexpr(is_text_v<U>) return false;
            if constexpr(std::is_array_v<DecayType>) return is_text_v<typename std::remove_all_extents_t<DecayType>>;
            if constexpr(std::is_pointer_v<DecayType>) return is_text_v<typename std::remove_pointer_t<DecayType>>;
            if constexpr(has_value_type_v<DecayType>) return is_text_v<typename DecayType::value_type>;
            return false;
        }

        public:
        static constexpr bool value = test<T>();
    };
    template<typename T>
    inline constexpr bool has_text_v = has_text<T>::value;

    template<typename T>
    struct is_Scalar2 {
        private:
        static constexpr bool test() {
            if constexpr(has_x_v<T> and has_y_v<T> and not has_z_v<T>) {
                constexpr size_t t_size = sizeof(T);
                constexpr size_t x_size = sizeof(T::x);
                constexpr size_t y_size = sizeof(T::y);
                return t_size == x_size + y_size;
            } else {
                return false;
            }
        }

        public:
        static constexpr bool value = test();
    };
    template<typename T>
    inline constexpr bool is_Scalar2_v = is_Scalar2<T>::value;

    template<typename T>
    struct is_Scalar3 {
        private:
        static constexpr bool test() {
            if constexpr(has_x_v<T> and has_y_v<T> and has_z_v<T>) {
                constexpr size_t t_size = sizeof(T);
                constexpr size_t x_size = sizeof(T::x);
                constexpr size_t y_size = sizeof(T::y);
                constexpr size_t z_size = sizeof(T::z);
                return t_size == x_size + y_size + z_size;
            } else {
                return false;
            }
        }

        public:
        static constexpr bool value = test();
    };
    template<typename T>
    inline constexpr bool is_Scalar3_v = is_Scalar3<T>::value;

    template<typename T1, typename T2>
    constexpr bool is_Scalar2_of_type() {
        if constexpr(is_Scalar2_v<T1>) return std::is_same<decltype(T1::x), T2>::value;
        return false;
    }

    template<typename T1, typename T2>
    constexpr bool is_Scalar3_of_type() {
        if constexpr(is_Scalar3_v<T1>)
            return std::is_same<decltype(T1::x), T2>::value;
        else
            return false;
    }

    template<typename T>
    constexpr bool is_ScalarN() {
        return is_Scalar2_v<T> or is_Scalar3_v<T>;
    }
#endif

}
