#pragma once
#ifndef F128_H
    #define F128_H
    #include <stdexcept>
    #include <string>
    #include <string_view>
    #include <type_traits>
    #if defined(USE_QUADMATH)
        #include <complex>
        #include <quadmath.h>
// Wrapper around the __float128 type
struct f128_t {
    private:
    __float128 val;

    public:
    using format_type = long double;
    using value_type  = __float128;
    f128_t()          = default;
    f128_t(std::string_view s) : val(strtoflt128(s.data(), nullptr)) {}
    f128_t(const char *c) : val(strtoflt128(c, nullptr)) {}
    template<typename T>
    f128_t(const T &v) : val(v) {}
    operator __float128() const { return val; }
    [[nodiscard]] __float128 value() const { return val; }
    template<typename T>
    bool operator==(const T &v) const {
        static_assert(std::is_arithmetic_v<T> or std::is_same_v<T, __float128>);
        if constexpr(std::is_same_v<T, f128_t>)
            return val == v.val;
        else
            return val == v;
    }
    template<typename T>
    f128_t operator=(T v) {
        val = v;
        return val;
    }

    template<typename T>
    f128_t operator*(const T &v) const {
        return val * v;
    }
    template<typename T>
    f128_t operator/(const T &v) const {
        return val / v;
    }
    template<typename T>
    f128_t operator+(const T &v) const {
        return val + v;
    }
    template<typename T>
    f128_t operator-(const T &v) const {
        return val - v;
    }
    [[nodiscard]] std::string string(int prec = 36, int width = 0, char pres = 'f', std::string_view align = "") const {
        if(prec < 0) prec = 36;
        if(width < 0) width = 0;
        std::string buf;
        char        fstr[16] = {0};
        auto        res = std::snprintf(fstr, sizeof fstr, "%%%s%d.%uQ%c", align == "<" ? "-" : align.data(), width, prec, pres);
        if(res < 0) throw std::runtime_error("f128_t.string(): snprintf() returned < 0");
        auto size = quadmath_snprintf(nullptr, 0, fstr, val);
        if(size >= 0) {
            buf.resize(static_cast<size_t>(size));
        } else {
            throw std::runtime_error("f128_t.string(): quadmath_snprintf() returned < 0");
        }
        quadmath_snprintf(buf.data(), static_cast<size_t>(size), fstr, val);
        return buf;
    }
};

struct c128_t : std::complex<f128_t> {
    using complex<f128_t>::complex;
    template<typename T>
    c128_t(const T &v) {
        this->real() = std::real(v);
        this->imag() = std::imag(v);
    }
};

    #else
// Dummy wrapper around the long double type
struct f128_t {
    private:
    long double val;

    public:
    using format_type = long double;
    using value_type  = long double;
    f128_t(std::string_view s) : val(strtold(s.data(), nullptr)) {}
    f128_t(const char *c) : val(strtold(c, nullptr)) {}
    template<typename T>
    f128_t(T v) : val(v) {}
    operator long double() const { return val; }
    template<typename T>
    bool operator==(T v) const {
        static_assert(std::is_arithmetic_v<T> or std::is_same_v<T, __float128>);
        if constexpr(std::is_same_v<T, __float128>)
            return val == v.val;
        else
            return val == v;
    }
    template<typename T>
    f128_t operator=(T v) {
        val = v;
        return val;
    }

    template<typename T>
    f128_t operator*(T v) const {
        return val * v;
    }
    template<typename T>
    f128_t operator/(T v) const {
        return val / v;
    }
    template<typename T>
    f128_t operator+(T v) const {
        return val + v;
    }
    template<typename T>
    f128_t operator-(T v) const {
        return val - v;
    }
    [[nodiscard]] std::string string(int prec = 36, int width = 0, char pres = 'f', std::string_view align = "") const {
        std::string buf;
        char        fstr[16];
        std::snprintf(fstr, sizeof fstr, "%%%s%d.%dL%c", align == "<" ? "-" : align.data(), width, prec, pres);
        auto size = std::snprintf(nullptr, 0, fstr, val);
        if(size >= 0) {
            buf.resize(static_cast<size_t>(size) + 1);
        } else {
            throw std::runtime_error("f128_t.string(): quadmath_snprintf() returned size < 0");
        }
        std::snprintf(buf.data(), buf.size(), fstr, val);
        return buf;
    }
    [[nodiscard]] long double value() const { return val; }
    bool                      operator==(f128_t other) const { return val == other.val; }
};
    #endif

#endif