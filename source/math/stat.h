#pragma once

#include <cmath>
#include <complex>
#include <iterator>
#include <numeric>
#include <optional>
#include <vector>
#include "num.h"
/*!
 *  \namespace stat
 *  \brief Small convenience-type statistical functions, mean and slope
 *  \tableofcontents
 */


namespace stat{

    template<typename ContainerType>
    void check_bounds(ContainerType &X, std::optional<long> start_point = std::nullopt, std::optional<long> end_point = std::nullopt) {
        if(start_point.has_value() and (start_point.value() >= static_cast<long>(X.size()) or start_point.value() < 0))
            throw std::range_error("Start point is out of range");
        if(end_point.has_value() and (end_point.value() > static_cast<long>(X.size()) or end_point.value() < start_point))
            throw std::range_error("End point is out of range");
    }

    template<typename ContainerType>
    double mean(ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        try {
            check_bounds(X, start_point, end_point);
        } catch(std::exception &err) { throw std::range_error("mean: " + std::string(err.what())); }
        if(not start_point.has_value()) start_point = 0;
        if(not end_point.has_value()) end_point = X.size();
        if(end_point.value() == start_point.value()) return 0.0;
        if(end_point.value() < start_point.value()) throw std::runtime_error("end_point < start_point");

        auto x_it_start = X.begin();
        auto x_it_end   = X.begin();
        auto n          = static_cast<double>(end_point.value() - start_point.value());
        std::advance(x_it_start, start_point.value());
        std::advance(x_it_end, end_point.value());
        return accumulate(x_it_start, x_it_end, 0.0) / n;
    }

    template<typename ContainerType>
    double stdev(const ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        try {
            check_bounds(X, start_point, end_point);
        } catch(std::exception &err) { throw std::range_error("stdev: " + std::string(err.what())); }
        if(not start_point.has_value()) start_point = 0;
        if(not end_point.has_value()) end_point = X.size();
        if(end_point.value() == start_point.value()) return 0.0;
        if(end_point.value() < start_point.value()) throw std::runtime_error("end_point < start_point");

        double X_mean = stat::mean(X, start_point, end_point);
        auto   n      = static_cast<double>(end_point.value() - start_point.value());
        auto   x_it   = X.begin();
        auto   x_en   = X.begin();

        std::advance(x_it, start_point.value());
        std::advance(x_en, end_point.value());

        double sum = std::accumulate(x_it, x_en, 0.0, [&X_mean](auto &x1, auto &x2) { return x1 + (x2 - X_mean) * (x2 - X_mean); });

        return std::sqrt(sum / n);
    }

    template<typename ContainerType>
    double sterr(const ContainerType &X, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        return stdev(X, start_point, end_point) / std::sqrt(X.size());
    }

    template<typename ContainerType1, typename ContainerType2>
    std::pair<double,double> slope(const ContainerType1 &X, const ContainerType2 &Y, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        if(X.size() != Y.size()) throw std::range_error("slope: size mismatch in arrays: X.size() == " + std::to_string(X.size()) + " | Y.size() == "  + std::to_string(Y.size()));
        try {
            check_bounds(X, start_point, end_point);
        } catch(std::exception &err) { throw std::range_error("slope: " + std::string(err.what())); }
        if(not start_point.has_value()) start_point = 0;
        if(not end_point.has_value()) end_point = X.size();
        if(end_point.value() < start_point.value()) throw std::runtime_error("end_point < start_point");
        if(end_point.value() - start_point.value() <= 1) return std::make_pair(0.0,0.0); // Need at least 2 points

        auto x_it = X.begin();
        auto x_en = X.begin();
        auto y_it = Y.begin();
        auto y_en = Y.begin();
        std::advance(x_it, start_point.value());
        std::advance(y_it, start_point.value());
        std::advance(x_en, end_point.value());
        std::advance(y_en, end_point.value());
        auto   n    = static_cast<double>(end_point.value() - start_point.value());
        double avgX = std::accumulate(x_it, x_en, 0.0) / n;
        double avgY = std::accumulate(y_it, y_en, 0.0) / n;
        double sxx = 0.0;
        double syy = 0.0;
        double sxy = 0.0;
        while(x_it != x_en) {
            sxx += std::pow((static_cast<double>(*x_it) - avgX),2);
            syy += std::pow((static_cast<double>(*y_it) - avgY),2);
            sxy += (static_cast<double>(*x_it) - avgX) * (static_cast<double>(*y_it) - avgY);
            y_it++;
            x_it++;
        }
        double slope = sxy / sxx;
        double residual = syy * (1 - sxy*sxy / sxx / syy);
        return std::make_pair(slope,residual);
    }

    template<typename ContainerType>
    size_t find_saturation_point(const ContainerType &Y, double slope_tolerance = 1, double std_tolerance = 1, std::optional<size_t> start_point = std::nullopt, std::optional<size_t> end_point = std::nullopt) {
        try {
            check_bounds(Y, start_point, end_point);
        } catch(std::exception &err) { throw std::range_error("find_saturation_point: " + std::string(err.what())); }
        if(not start_point.has_value()) start_point = 0;
        if(not end_point.has_value()) end_point = Y.size();
        if(end_point.value() == start_point.value()) return start_point.value();
        if(end_point.value() < start_point.value()) throw std::runtime_error("find_saturation_point: end_point < start_point");
        printf("checking saturation from: %ld to slp: %ld\n", start_point.value(),end_point.value());

        // Sometimes we check saturation using logarithms. It is important that the start-end point range doesn't have nan or infs

        auto y_it = Y.begin();
        auto y_en = Y.begin();
        std::advance(y_it, start_point.value());
        std::advance(y_en, end_point.value());
        auto y_invalid_it = std::find_if(y_it,y_en, [](auto & val){return std::isinf(val) or std::isnan(val);});
        if(y_invalid_it <= y_en){
            // Found an invalid point! Back-track once
            std::advance(y_invalid_it, -1);
            if(y_it == y_invalid_it) return start_point.value();
            else end_point = std::distance(Y.begin(), y_invalid_it);
            printf("checking saturation from: %ld to slp: %ld instead\n", start_point.value(),end_point.value());
        }


        // Consider Y vs X: a noisy signal decaying in the shape of a hockey-club, say.
        // We want to identify the point at which the signal stabilizes. We use the fact that the
        // standard deviation is high if it includes parts of the non-stable signal, and low if
        // it includes only the stable part.
        // Here we monitor the standard deviation of the signal between [start_point, end_point],
        // and move "start_point" towards the end. If the standard deviation goes below a certain
        // threshold, i.e. threshold < max_std, then we have found the stabilization point.

        size_t idx = start_point.value();
        auto X = num::range<size_t>(0, Y.size());
        while(idx < end_point.value()){
            auto std = stat::stdev(Y,idx,end_point.value());
            auto [slp,res] = stat::slope(X,Y,idx,end_point.value());
            printf("std: %g | slp: %g\n", std,slp);
            if(std::abs(slp) < slope_tolerance or std < std_tolerance){
                if(idx > 0) return idx-1; // Backtrack
                else return idx;
            }
            idx++;
        }
        return end_point.value()-1;
    }
}