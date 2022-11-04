#pragma once

#include "general/iter.h"
#include "num.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <numeric>
#include <optional>
#include <vector>

/*!
 *  \namespace stat
 *  \brief Small convenience-type statistical functions, mean and slope
 *  \tableofcontents
 */

namespace stat {
    template<typename T, typename = std::void_t<>>
    struct is_iterable : public std::false_type {};
    template<typename T>
    struct is_iterable<T, std::void_t<decltype(std::begin(std::declval<T>())), decltype(std::end(std::declval<T>()))>> : public std::true_type {};
    template<typename T>
    inline constexpr bool is_iterable_v = is_iterable<T>::value;

    template<typename ContainerType>
    void check_bounds(ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        if(start_point.has_value() and num::cmp_greater_equal(start_point.value(), std::size(X))) throw std::range_error("Start point is out of range");
        if(end_point.has_value() and (num::cmp_greater(end_point.value(), std::size(X)) or end_point.value() < start_point))
            throw std::range_error("End point is out of range");
    }

    template<typename ContainerType>
    [[nodiscard]] auto get_start_end_iterators(ContainerType &X, std::optional<size_t> start_point = std::nullopt,
                                               std::optional<size_t> end_point = std::nullopt) {
        if(end_point.has_value() and end_point.value() == -1ul) end_point = X.size();
        try {
            check_bounds(X, start_point, end_point);
        } catch(std::exception &err) { throw std::range_error("check_bounds failed: " + std::string(err.what())); }
        if(not start_point.has_value()) start_point = 0;
        if(not end_point.has_value()) end_point = X.size();
        return std::make_pair(X.begin() + start_point.value(), X.begin() + end_point.value());
    }

    template<typename ContainerType>
    [[nodiscard]] auto min(ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        auto [x_it, x_en] = get_start_end_iterators(X, start_point, end_point);
        return *std::min_element(x_it, x_en);
    }
    template<typename ContainerType>
    [[nodiscard]] auto max(ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        auto [x_it, x_en] = get_start_end_iterators(X, start_point, end_point);
        return *std::max_element(x_it, x_en);
    }

    template<typename ContainerType>
    [[nodiscard]] typename ContainerType::value_type mean(ContainerType &X, std::optional<size_t> start_point = std::nullopt,
                                                          std::optional<size_t> end_point = std::nullopt) {
        auto [x_it, x_en] = get_start_end_iterators(X, start_point, end_point);
        auto n            = static_cast<double>(std::distance(x_it, x_en));
        return std::accumulate(x_it, x_en, static_cast<typename ContainerType::value_type>(0.0)) / n;
    }

    template<typename ContainerType>
    [[nodiscard]] typename ContainerType::value_type median(ContainerType &X, std::optional<size_t> start_point = std::nullopt,
                                                            std::optional<size_t> end_point = std::nullopt) {
        auto [x_it, x_en] = get_start_end_iterators(X, start_point, end_point);
        auto n            = static_cast<double>(std::distance(x_it, x_en));
        if(n == 0) return 0;
        // Find the starting point
        auto x_it_mid = X.begin();
        std::advance(x_it_mid, start_point.value());
        // ... and then the mid-point between start to end
        std::advance(x_it_mid, static_cast<long>(n / 2));

        if(n > 0 and num::mod<size_t>(static_cast<size_t>(n), 2) == 0) {
            // Even number of elements, take mean between middle elements
            auto a = *x_it_mid;
            std::advance(x_it_mid, -1);
            auto b = *x_it_mid;
            return 0.5 * (a + b);
        } else {
            // odd number of elements
            return *x_it_mid;
        }
    }

    template<typename ContainerType>
    [[nodiscard]] auto typical(const ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        ContainerType Xlog = X;
        for(auto &x : Xlog) x = std::log(std::abs(x));
        return std::exp(mean(Xlog, start_point, end_point));
    }
    template<typename ContainerType>
    [[nodiscard]] double variance(const ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        auto [x_it, x_en] = get_start_end_iterators(X, start_point, end_point);
        auto n            = static_cast<double>(std::distance(x_it, x_en));
        if(n == 0) return 0.0;
        double X_mean = stat::mean(X, start_point, end_point);
        double sum    = std::accumulate(x_it, x_en, 0.0, [&X_mean](auto &x1, auto &x2) { return x1 + (x2 - X_mean) * (x2 - X_mean); });
        return sum / n;
    }
    template<typename ContainerType>
    [[nodiscard]] double stdev(const ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        return std::sqrt(variance(X, start_point, end_point));
    }

    template<typename ContainerType>
    [[nodiscard]] double sterr(const ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        auto [x_it, x_en] = get_start_end_iterators(X, start_point, end_point);
        if(x_it == x_en) return 0.0;
        auto n = static_cast<double>(std::distance(x_it, x_en));
        return stdev(X, start_point, end_point) / std::sqrt(n);
    }

    template<typename ContainerType>
    [[nodiscard]] std::vector<double> sterr_moving(const ContainerType &X, double width = 0.1, std::optional<size_t> start_point = std::nullopt,
                                                   std::optional<size_t> end_point = std::nullopt) {
        std::vector<double> res;
        width             = std::clamp<double>(width, 0.0, 1.0);
        auto [x_it, x_en] = get_start_end_iterators(X, start_point, end_point);
        long idx0         = std::distance(X.begin(), x_it);
        long idxN         = std::distance(X.begin(), x_en);
        long xlen         = std::distance(x_it, x_en);                                               // Number of points to consider in X
        long wlen         = std::max<long>(2, static_cast<long>(width * static_cast<double>(xlen))); // Number of points in the moving window
        while(x_it != x_en) {
            long idx_beg = std::distance(X.begin(), x_it); // Center around x_it
            long idx_end = std::distance(X.begin(), x_it); // Center around x_it
            while(idx_end - idx_beg < wlen) {
                idx_end = std::clamp(idx_end + 1, idx_end, idxN); // Move forward 1 step if possible.
                if(idx_end - idx_beg < wlen) {
                    idx_beg = std::clamp(idx_beg - 1, idx0, idx_beg); // Move back 1 step if possible.
                }
                if(idx_beg == idx0 and idx_end == idxN) break;
            }
            if(idx_end < 0 or idx_end > X.size()) break;
            res.emplace_back(sterr(X, idx_beg, idx_end));
            x_it++;
        }
        return res;
    }

    template<typename ContainerType1, typename ContainerType2, typename = std::enable_if_t<is_iterable_v<ContainerType1> and is_iterable_v<ContainerType2>>>
    [[nodiscard]] std::vector<double> cumtrapz(const ContainerType1 &X, const ContainerType2 &Y, std::optional<size_t> start_point = std::nullopt,
                                               std::optional<size_t> end_point = std::nullopt) {
        auto [x_it, x_en]       = get_start_end_iterators(X, start_point, end_point);
        auto [y_it, y_en]       = get_start_end_iterators(Y, start_point, end_point);
        auto                num = std::min(std::distance(x_it, x_en), std::distance(y_it, y_en));
        std::vector<double> res;
        res.reserve(num);
        for(long i = 0; i < num; i++) { res.emplace_back(num::trapz(X, Y, i, num - i)); }
        return res;
    }
    template<typename ContainerType1, typename ContainerType2, typename = std::enable_if_t<is_iterable_v<ContainerType1> and is_iterable_v<ContainerType2>>>
    [[nodiscard]] std::vector<double> cumtrapz_avg(const ContainerType1 &X, const ContainerType2 &Y, std::optional<size_t> start_point = std::nullopt,
                                                   std::optional<size_t> end_point = std::nullopt) {
        auto [x_it, x_en]           = get_start_end_iterators(X, start_point, end_point);
        auto [y_it, y_en]           = get_start_end_iterators(Y, start_point, end_point);
        auto                num     = std::min(std::distance(x_it, x_en), std::distance(y_it, y_en));
        auto                idx_end = std::min(std::distance(X.begin(), x_en), std::distance(Y.begin(), y_en));
        std::vector<double> res;
        res.reserve(num);
        while(x_it != x_en and y_it != y_en) {
            auto idx_start = std::min(std::distance(X.begin(), x_it), std::distance(Y.begin(), y_it));
            auto div       = static_cast<double>(X[idx_end] - X[idx_start]);
            num            = std::min(std::distance(x_it, x_en), std::distance(y_it, y_en));
            res.emplace_back(num::trapz(X, Y, idx_start, num) / div);
            x_it++;
            y_it++;
        }
        return res;
    }

    template<typename ContainerType1, typename ContainerType2, typename = std::enable_if_t<is_iterable_v<ContainerType1> and is_iterable_v<ContainerType2>>>
    [[nodiscard]] std::pair<double, double> slope(const ContainerType1 &X, const ContainerType2 &Y, std::optional<size_t> start_point = std::nullopt,
                                                  std::optional<size_t> end_point = std::nullopt) {
        if(X.size() != Y.size())
            throw std::range_error("slope: size mismatch in arrays: X.size() == " + std::to_string(X.size()) + " | Y.size() == " + std::to_string(Y.size()));
        auto [x_it, x_en] = get_start_end_iterators(X, start_point, end_point);
        auto [y_it, y_en] = get_start_end_iterators(Y, start_point, end_point);
        auto n            = static_cast<double>(std::distance(x_it, x_en));
        if(n <= 1) return std::make_pair(std::numeric_limits<double>::infinity(), 0.0); // Need at least 2 points
        double avgX = std::accumulate(x_it, x_en, 0.0) / n;
        double avgY = std::accumulate(y_it, y_en, 0.0) / n;
        double sxx  = 0.0;
        double syy  = 0.0;
        double sxy  = 0.0;
        while(x_it != x_en) {
            sxx += std::pow((static_cast<double>(*x_it) - avgX), 2);
            syy += std::pow((static_cast<double>(*y_it) - avgY), 2);
            sxy += (static_cast<double>(*x_it) - avgX) * (static_cast<double>(*y_it) - avgY);
            y_it++;
            x_it++;
        }
        double slope    = sxy / sxx;
        double residual = syy * (1 - sxy * sxy / sxx / syy);
        return std::make_pair(slope, residual);
    }

    template<typename ContainerType>
    [[nodiscard]] std::pair<double, double> slope(const ContainerType &Y, std::optional<size_t> start_point = std::nullopt,
                                                  std::optional<size_t> end_point = std::nullopt) {
        ContainerType X(Y.size());
        std::iota(X.begin(), X.end(), 0);
        return slope(X, Y, start_point, end_point);
    }

    template<typename ContainerType1, typename ContainerType2>
    [[nodiscard]] double slope_at(const ContainerType1 &X, const ContainerType2 &Y, size_t at, size_t width = 1) {
        auto min_idx = static_cast<size_t>(std::max(static_cast<long>(at) - static_cast<long>(width), 0l));
        auto max_idx = static_cast<size_t>(std::min(at + width, Y.size()));
        return slope(X, Y, min_idx, max_idx).first;
    }

    template<typename ContainerType>
    [[nodiscard]] ContainerType smooth(const ContainerType &X, long width = 2) {
        if(X.size() <= 2) return X;
        width = std::min<long>(4, static_cast<long>(X.size()) / 2);
        ContainerType S;
        S.reserve(X.size());
        for(auto &&[i, x] : iter::enumerate(X)) {
            auto min_idx = std::clamp<long>(static_cast<long>(i) - width, 0l, static_cast<long>(i));
            auto max_idx = std::clamp<long>(static_cast<long>(i) + width, static_cast<long>(i), static_cast<long>(X.size()) - 1l);
            S.push_back(stat::mean(X, min_idx, max_idx));
        }
        return S;
    }

    template<typename ContainerType>
    [[nodiscard]] size_t find_last_valid_point(const ContainerType &Y, std::optional<size_t> start_point = std::nullopt,
                                               std::optional<size_t> end_point = std::nullopt) {
        auto [y_it, y_en] = get_start_end_iterators(Y, start_point, end_point);
        auto n            = static_cast<double>(std::distance(y_it, y_en));
        if(n == 0) return start_point.has_value() ? start_point.value() : 0;
        auto y_invalid_it = std::find_if(y_it, y_en, [](auto &val) { return std::isinf(val) or std::isnan(val); });
        if(y_invalid_it != y_en and y_invalid_it != y_it) {
            // Found an invalid point! Back-track once
            std::advance(y_invalid_it, -1);
        }
        end_point = std::distance(Y.begin(), y_invalid_it);
        return end_point.value();
    }

    template<typename ContainerType>
    [[nodiscard]] size_t find_saturation_point(const ContainerType &Y, double slope_tolerance = 1, double std_tolerance = 1,
                                               std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        auto [y_it, y_en] = get_start_end_iterators(Y, start_point, end_point);
        auto n            = static_cast<double>(std::distance(y_it, y_en));
        auto idx          = start_point.has_value() ? start_point.value() : 0;
        if(n == 0) return idx;

        // Sometimes we check saturation using logarithms. It is important that the start-end point range doesn't have nan or infs
        end_point = find_last_valid_point(Y, start_point, end_point);

        // Consider Y vs X: a noisy signal decaying in the shape of a hockey-club, say.
        // We want to identify the point at which the signal stabilizes. We use the fact that the
        // standard deviation is high if it includes parts of the non-stable signal, and low if
        // it includes only the stable part.
        // Here we monitor the standard deviation of the signal between [start_point, end_point],
        // and move "start_point" towards the end. If the standard deviation goes below a certain
        // threshold, i.e. threshold < max_std, then we have found the stabilization point.

        auto X = num::range<size_t>(0, Y.size());
        while(idx < end_point.value()) {
            auto std        = stat::stdev(Y, idx, end_point.value());
            auto [slp, res] = stat::slope(X, Y, idx, end_point.value());
            printf("std: %g | slp: %g\n", std, slp);
            if(std::abs(slp) < slope_tolerance or std < std_tolerance) {
                if(idx > 0) {
                    printf("found idx %ld\n", idx - 1);
                    return idx - 1; // Backtrack
                } else {
                    printf("found idx %ld\n", idx);
                    return idx;
                }
            }
            idx++;
        }
        return end_point.value() - 1;
    }
}