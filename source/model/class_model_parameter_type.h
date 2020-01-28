#pragma once
#include <string>
template<typename T>
class ParameterType {
    protected:
    const std::string name_;
    T                 value_;

    public:
    ParameterType() = delete;
    ParameterType(const ParameterType &other) : name_(other.name()), value_(other.value()) {}
    ParameterType(std::string_view name, const T &value) : name_(name), value_(value) {}
    explicit ParameterType(const std::pair<std::string_view, T> &pair) : name_(pair.first), value_(pair.second) {}

    using value_type = std::decay<T>;
    [[nodiscard]] const std::string &name() const { return name_; }
    [[nodiscard]] const T &          value() const { return value_; }

    template<typename RHS, typename = std::enable_if_t<std::is_arithmetic_v<RHS>>>
    ParameterType &operator=(const RHS &rhs) {
        value_ = rhs;
        return *this;
    }

    ParameterType &operator=(const ParameterType &rhs) {
        if(this == &rhs) return *this;
        value_ = rhs.value();
        return *this;
    }

    template<typename RHS>
    bool operator==(const RHS &rhs) const {
        return value_ == rhs;
    }
    template<typename RHS>
    bool operator!=(const RHS &rhs) const {
        return value_ != rhs;
    }
    template<typename RHS>
    bool operator<=(const RHS &rhs) const {
        return value_ <= rhs;
    }
    template<typename RHS>
    bool operator>=(const RHS &rhs) const {
        return value_ >= rhs;
    }
    template<typename RHS>
    bool operator<(const RHS &rhs) const {
        return value_ < rhs;
    }
    template<typename RHS>
    bool operator>(const RHS &rhs) const {
        return value_ > rhs;
    }
    operator T() const { return value_; }
    //    virtual             operator std::string() const = 0;
    explicit             operator std::string_view() const { return std::string_view(value_); }
    friend std::ostream &operator<<(std::ostream &os, const ParameterType &p) { return os << p.value_; }
};

class ParamDouble : public ParameterType<double> {
    public:
    explicit operator std::string() const { return std::to_string(value_); }
};

class ParamString : public ParameterType<std::string> {};

class ParamInteger : public ParameterType<int> {
    public:
    explicit operator std::string() const { return std::to_string(value_); }
};
class ParamBool : public ParameterType<bool> {
    public:
    explicit operator std::string() const { return std::to_string(value_); }
};
class ParamSize_t : public ParameterType<size_t> {
    public:
    explicit operator std::string() const { return std::to_string(value_); }
};