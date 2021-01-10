#pragma once

namespace eig::sfinae {
    // A bit of sfinae deduction
    template<typename T, typename = std::void_t<>>
    struct has_RawEigenvaluesImag : public std::false_type {};
    template<typename T>
    struct has_RawEigenvaluesImag<T, std::void_t<decltype(std::declval<T>().size())>> : public std::true_type {};
    template<typename T>
    inline constexpr bool has_RawEigenvaluesImag_v = has_RawEigenvaluesImag<T>::value;
}