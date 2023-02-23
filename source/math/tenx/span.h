#pragma once
#include "sfinae.h"
#include <iterator>

namespace tenx {

    template<typename T>
    class span {
        T          *ptr_;
        std::size_t len_;

        public:
        using value_type = T;
        template<typename size_type>
        span(T *ptr, size_type len) noexcept : ptr_{ptr}, len_{static_cast<std::size_t>(len)} {
            static_assert(std::is_integral_v<size_type>);
        }
        span(T *bgn, T *end) noexcept : ptr_{bgn}, len_{static_cast<std::size_t>(std::distance(bgn, end))} {}
        template<auto rank, int options>
        span(Eigen::Tensor<T, rank, options> &t) noexcept : ptr_{t.data()}, len_{static_cast<std::size_t>(t.size())} {}
        span(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &m) noexcept : ptr_{m.data()}, len_{static_cast<std::size_t>(m.size())} {}
        span(Eigen::Matrix<T, Eigen::Dynamic, 1> &m) noexcept : ptr_{m.data()}, len_{static_cast<std::size_t>(m.size())} {}
        span(Eigen::Matrix<T, 1, Eigen::Dynamic> &m) noexcept : ptr_{m.data()}, len_{static_cast<std::size_t>(m.size())} {}
        span(Eigen::Array<T, Eigen::Dynamic, 1> &a) noexcept : ptr_{a.data()}, len_{static_cast<std::size_t>(a.size())} {}
        span(Eigen::Array<T, 1, Eigen::Dynamic> &a) noexcept : ptr_{a.data()}, len_{static_cast<std::size_t>(a.size())} {}
        span(std::vector<T> &v) noexcept : ptr_{v.data()}, len_{v.size()} {}

        T &operator[](int i) noexcept { return *ptr_[i]; }

        T const &operator[](int i) const noexcept { return *ptr_[i]; }

        [[nodiscard]] std::size_t size() const noexcept { return len_; }

        T       *data() noexcept { return ptr_; }
        const T *data() const noexcept { return ptr_; }
        T       *begin() noexcept { return ptr_; }
        T       *end() noexcept { return ptr_ + len_; }
        const T *begin() const noexcept { return ptr_; }
        const T *end() const noexcept { return ptr_ + len_; }
    };
}
