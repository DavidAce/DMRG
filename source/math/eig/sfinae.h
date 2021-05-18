#pragma once

namespace eig::sfinae {
    // A bit of sfinae deduction
    template<typename T, typename = std::void_t<>>
    struct has_RawEigenvaluesImag : public std::false_type {};
    template<typename T>
    struct has_RawEigenvaluesImag<T, std::void_t<decltype(std::declval<T>().size())>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_RawEigenvaluesImag_v = has_RawEigenvaluesImag<T>::value;

    template<typename T>
    constexpr auto type_name() {
        std::string_view name, prefix, suffix;
#ifdef __clang__
        name   = __PRETTY_FUNCTION__;
        prefix = "auto eig::sfinae::type_name() [T = ";
        suffix = "]";
#elif defined(__GNUC__)
        name   = __PRETTY_FUNCTION__;
        prefix = "constexpr auto eig::sfinae::type_name() [with T = ";
        suffix = "]";
#elif defined(_MSC_VER)
        name   = __FUNCSIG__;
        prefix = "auto __cdecl eig::sfinae::type_name<";
        suffix = ">(void)";
#endif
        name.remove_prefix(prefix.size());
        name.remove_suffix(suffix.size());
        return name;
    }

}