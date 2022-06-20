#include "AlgorithmBase.h"
#include "config/settings.h"
#include "general/iter.h"
#include "math/num.h"
#include "math/stat.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include <complex>
#include <debug/exceptions.h>
#include <h5pp/h5pp.h>
using Scalar = AlgorithmBase::Scalar;

AlgorithmBase::AlgorithmBase(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type_) : h5file(std::move(h5ppFile_)) {
    status.algo_type = algo_type_;
    tools::log->set_error_handler([](const std::string &msg) { throw except::runtime_error(msg); });
    tools::log = tools::Logger::setLogger(status.algo_type_str(), settings::console::loglevel, settings::console::timestamp);
    tools::log->trace("Constructing class_algorithm_base");
}

void AlgorithmBase::copy_from_tmp(StorageReason storage_reason, std::optional<CopyPolicy> copy_policy) {
    if(h5file) tools::common::h5::tmp::copy_from_tmp(status, *h5file, storage_reason, copy_policy);
}

void AlgorithmBase::init_bond_dimension_limits() {
    status.bond_init = settings::get_bond_init(status.algo_type);
    status.bond_max  = settings::get_bond_max(status.algo_type);
    if(settings::strategy::bond_grow_mode != BondGrow::OFF)
        status.bond_lim = settings::get_bond_init(status.algo_type);
    else
        status.bond_lim = settings::get_bond_max(status.algo_type);

    // Sanity check
    if(status.bond_lim == 0) throw except::runtime_error("Bond dimension limit invalid: {}", status.bond_lim);
}

void AlgorithmBase::write_enable() { write_enabled = true; }
void AlgorithmBase::write_disable() { write_enabled = false; }

size_t AlgorithmBase::count_convergence(const std::vector<double> &Y_vec, double threshold, size_t start_idx) {
    size_t scount = 0; // Counts how many converged points there have been since saturation (start_idx)
    for(const auto &[i, y] : iter::enumerate(Y_vec)) {
        if(i < start_idx) continue;
        if(y <= threshold) scount++;
    }
    size_t rcount = 0; // Counts in reverse how many converged points there have been in a row. Useful with noisy signals that can't saturate.
    for(const auto &y : iter::reverse(Y_vec)) {
        if(y <= threshold)
            rcount++;
        else
            break;
    }
    return std::max(scount, rcount);
}

AlgorithmBase::SaturationReport AlgorithmBase::check_saturation(const std::vector<double> &Y_vec, double sensitivity) {
    SaturationReport report;
    constexpr size_t min_data_points = 2;
    if(Y_vec.size() < min_data_points) { return report; }

    // Running average [i:end]
    std::vector<double> Y_avg;
    Y_avg.reserve(Y_vec.size());
    for(const auto &[i, y] : iter::enumerate(Y_vec)) Y_avg.push_back(stat::mean(Y_vec, i));

    // Get the standard deviations from i to end
    std::vector<double> Y_std;
    Y_std.reserve(Y_avg.size());
    for(const auto &[i, y] : iter::enumerate(Y_avg)) Y_std.push_back(stat::stdev(Y_avg, i));

    // "Normalize" the standard deviations so that this becomes scale invariant w.r.t Y_vec
    std::vector<double> Y_stn;
    Y_stn.reserve(Y_std.size());
    for(const auto &[i, y] : iter::enumerate(Y_std)) {
        double divisor = Y_avg[i] == 0.0 ? 1.0 : Y_avg[i];
        Y_stn.push_back(Y_std[i] / divisor);
    }

    // Get also the slopes of the smoothened data
    std::vector<double> Y_slp, Y_log;
    Y_slp.reserve(Y_vec.size());
    Y_log.reserve(Y_vec.size());
    auto Y_smt = stat::smooth(Y_vec, 2);
    for(auto &y : Y_smt) Y_log.push_back(-std::log10(std::abs(y)));
    // Normalize so the last element is 1
    double yback = Y_log.back();
    for(auto &y : Y_log) y /= yback;
    for(const auto &[i, y] : iter::enumerate(Y_vec)) {
        auto [slp, res] = stat::slope(Y_log, i);
        Y_slp.push_back(std::abs(slp));
    }

    size_t saturated_from_idx = 0;
    for(const auto &[i, a] : iter::enumerate(Y_avg)) {
        saturated_from_idx = i;
        auto median        = stat::median(Y_avg, i);
        auto bwidth        = 10 * Y_std[i]; // Band width
                                            //        tools::log->info("Y_vec[{:3}] = {:7.4e} | band = {:7.4e} +- {:7.4e}",i, Y_vec[i], median,bwidth);
        bool rel_cond = Y_stn[i] < sensitivity;
        bool abs_cond = Y_std[i] < 1e-10;
        bool win_cond = Y_vec[i] == std::clamp(Y_vec[i], median - bwidth, median + bwidth);
        bool slp_cond = Y_slp[i] < 0.1 * sensitivity;
        if((rel_cond or abs_cond or slp_cond) and win_cond) break;
    }

    report.has_computed    = true;
    report.saturated_point = saturated_from_idx;
    report.saturated_count = Y_vec.size() - saturated_from_idx - 1;
    report.has_saturated   = report.saturated_count > 0;
    report.Y_avg           = Y_avg;
    report.Y_vec           = Y_vec;
    report.Y_std           = Y_std;
    report.Y_stn           = Y_stn;
    report.Y_slp           = Y_slp;
    return report;
}
