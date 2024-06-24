#include "AlgorithmBase.h"
#include "config/settings.h"
#include "debug/exceptions.h"
#include "general/iter.h"
#include "math/cast.h"
#include "math/num.h"
#include "math/stat.h"
#include "tools/common/h5.h"
#include "tools/common/log.h"
#include <complex>
#include <h5pp/h5pp.h>
AlgorithmBase::AlgorithmBase(OptRitz opt_ritz_, AlgorithmType algo_type_) {
    status.opt_ritz  = opt_ritz_;
    status.algo_type = algo_type_;
    tools::log->set_error_handler([](const std::string &msg) { throw except::runtime_error(msg); });
    tools::log = tools::Logger::setLogger(fmt::format("{}", status.algo_type_sv()), settings::console::loglevel, settings::console::timestamp);
    tools::log->trace("Constructing class_algorithm_base");
    if(settings::test_unwind) throw std::runtime_error("Testing stack unwinding");
}
AlgorithmBase::AlgorithmBase(std::shared_ptr<h5pp::File> h5file_, OptRitz opt_ritz_, AlgorithmType algo_type_) : AlgorithmBase(opt_ritz_, algo_type_) {
    h5file = h5file_;
}

void AlgorithmBase::copy_from_tmp(StorageEvent storage_event, CopyPolicy copy_policy) {
    tools::common::h5::tmp::copy_from_tmp(*h5file, status.iter, status.step, storage_event, copy_policy);
}

void AlgorithmBase::init_bond_dimension_limits() {
    status.bond_max = settings::get_bond_max(status.algo_type); // Bond max from the loaded config file
    if(has_any_flags(status.algo_type, AlgorithmType::fDMRG, AlgorithmType::xDMRG, AlgorithmType::fLBIT)) {
        // Finite systems have limited bond dimension!
        // Be careful not to overflow long when the systems are very large (e.g. > 128 sites)
        double long_max = static_cast<double>(std::numeric_limits<long>::max());
        double bond_max = std::min(long_max, std::pow(2.0, settings::model::model_size / 2));
        status.bond_max      = std::min(status.bond_max, safe_cast<long>(bond_max));
    }
    status.bond_lim = std::min(status.bond_max, settings::get_bond_min(status.algo_type));
    status.bond_min = std::min(status.bond_max, settings::get_bond_min(status.algo_type));
    // Sanity check
    if(status.bond_lim == 0) throw except::runtime_error("Bond dimension limit invalid: {}", status.bond_lim);
}

void AlgorithmBase::init_truncation_error_limits() {
    if(settings::strategy::trnc_decrease_when == UpdatePolicy::NEVER) {
        status.trnc_min = settings::precision::svd_truncation_lim;
        status.trnc_lim = settings::precision::svd_truncation_lim;
        status.trnc_max = settings::precision::svd_truncation_lim;
    } else {
        status.trnc_min = settings::precision::svd_truncation_lim;
        status.trnc_lim = settings::precision::svd_truncation_init;
        status.trnc_max = settings::precision::svd_truncation_init;
    }
    // Sanity check
    if(status.trnc_lim == 0.0) throw except::runtime_error("Truncation error limit invalid: {}", status.trnc_lim);
    tools::log->info("Initialized truncation error limits: init {:8.2e} lim {:8.2e} min {:8.2e}", status.trnc_max, status.trnc_lim, status.trnc_min);
}

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

AlgorithmBase::SaturationReport AlgorithmBase::check_saturation(const std::vector<double> &Y_vec, double sensitivity, SaturationScale scale) {
    SaturationReport report;
    constexpr size_t min_data_points = 2;
    if(Y_vec.size() < min_data_points) { return report; }

    report.Y_vec = Y_vec;

    // Rescale
    if(scale == SaturationScale::log)
        for(auto &y : report.Y_vec) y = std::log10(y);

    // Get the cumulative max/min of the rescaled signal
    report.Y_max = num::cummax(report.Y_vec);
    report.Y_min = num::cummin(report.Y_vec);
    report.Y_sat.resize(report.Y_vec.size());

    if(not num::all_equal(report.Y_vec.size(), report.Y_min.size(), report.Y_max.size(), report.Y_sat.size()))
        throw except::logic_error("Report vectors are not equal size:\n Y_vec {}\n Y_min {}\n Y_max {}\n Y_sat {}", report.Y_vec.size(), report.Y_min.size(),
                                  report.Y_max.size(), report.Y_sat.size());
    auto Y_avg_full    = stat::mean(report.Y_vec);
    auto Y_avg_left    = stat::mean(report.Y_vec, 0, report.Y_vec.size() / 2);
    auto Y_avg_rite    = stat::mean(report.Y_vec, report.Y_vec.size() / 2);
    bool is_decreasing = Y_avg_left > Y_avg_full and Y_avg_full > Y_avg_rite;
    bool is_increasing = Y_avg_left < Y_avg_full and Y_avg_full < Y_avg_rite;
    bool is_uncreasing = not is_increasing and not is_decreasing; // May be more complicated
    // Calculate the running standard deviation from i to end
    report.Y_vec_std.resize(report.Y_vec.size());
    report.Y_min_std.resize(report.Y_min.size());
    report.Y_max_std.resize(report.Y_max.size());
    for(size_t i = 0; i < report.Y_vec.size() - 1; ++i) report.Y_vec_std[i] = stat::stdev(report.Y_vec, i);
    for(size_t i = 0; i < report.Y_min.size() - 1; ++i) report.Y_min_std[i] = stat::stdev(report.Y_min, i);
    for(size_t i = 0; i < report.Y_max.size() - 1; ++i) report.Y_max_std[i] = stat::stdev(report.Y_max, i);
    report.Y_vec_std.back() = report.Y_vec_std.rbegin()[1]; // The last entry is always zero, so just reuse the penultimate entry
    report.Y_min_std.back() = report.Y_min_std.rbegin()[1]; // The last entry is always zero, so just reuse the penultimate entry
    report.Y_max_std.back() = report.Y_max_std.rbegin()[1]; // The last entry is always zero, so just reuse the penultimate entry
    report.Y_mov_ste        = stat::sterr_moving(report.Y_vec, 0.25);

    for(size_t i = 0; i < report.Y_sat.size(); ++i) {
        // Below sensitivity --> saturated
        bool vec_sat = report.Y_vec_std[i] < sensitivity;
        bool min_sat = report.Y_min_std[i] < sensitivity;
        bool max_sat = report.Y_max_std[i] < sensitivity;
        bool mov_sat = report.Y_mov_ste[i] < sensitivity;
        if(is_decreasing) report.Y_sat[i] = vec_sat or mov_sat or min_sat;
        if(is_increasing) report.Y_sat[i] = vec_sat or mov_sat or max_sat;
        if(is_uncreasing) report.Y_sat[i] = vec_sat or mov_sat;
    }

    // Since the last element in Y_vec is always zero, Y_sat is always zero. we just copy the saturation state of the second to last element.
    // if(report.Y_sat.size() > 1) { report.Y_sat.back() = report.Y_sat.rbegin()[1]; }

    // From the end, count how many Y_sat[i] are 1,  before finding a 0.
    report.saturated_point = report.Y_sat.size() - 1;
    for(const auto &[i, y] : iter::enumerate_reverse(report.Y_sat)) {
        if(y == 0) break;
        report.saturated_point = i;
        report.saturated_count++;
    }
    report.has_computed    = true;
    report.saturated_scale = scale;
    report.has_saturated   = report.saturated_count > 0;
    return report;
}
