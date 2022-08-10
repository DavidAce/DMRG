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
    if(settings::strategy::bond_increase_when == UpdateWhen::NEVER) {
        status.bond_max  = settings::get_bond_max(status.algo_type);
        status.bond_lim  = settings::get_bond_max(status.algo_type);
        status.bond_init = settings::get_bond_max(status.algo_type);
    } else {
        status.bond_max  = settings::get_bond_max(status.algo_type);
        status.bond_lim  = settings::get_bond_init(status.algo_type);
        status.bond_init = settings::get_bond_init(status.algo_type);
    }
    // Sanity check
    if(status.bond_lim == 0) throw except::runtime_error("Bond dimension limit invalid: {}", status.bond_lim);
}

void AlgorithmBase::init_truncation_error_limits() {
    if(settings::strategy::trnc_decrease_when == UpdateWhen::NEVER) {
        status.trnc_min  = settings::precision::svd_truncation_lim;
        status.trnc_lim  = settings::precision::svd_truncation_lim;
        status.trnc_init = settings::precision::svd_truncation_lim;
    } else {
        status.trnc_min  = settings::precision::svd_truncation_lim;
        status.trnc_lim  = settings::precision::svd_truncation_init;
        status.trnc_init = settings::precision::svd_truncation_init;
    }
    // Sanity check
    if(status.trnc_lim == 0.0) throw except::runtime_error("Truncation error limit invalid: {}", status.trnc_lim);
    tools::log->info("Initialized truncation error limits: init {:8.2e} lim {:8.2e} min {:8.2e}", status.trnc_init, status.trnc_lim, status.trnc_min);
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

AlgorithmBase::SaturationReport AlgorithmBase::check_saturation(const std::vector<double> &Y_vec, double sensitivity, SaturationScale scale) {
    SaturationReport report;
    constexpr size_t min_data_points = 2;
    if(Y_vec.size() < min_data_points) { return report; }

    report.Y_vec = Y_vec;
    if(scale == SaturationScale::log)
        for(auto &y : report.Y_vec) y = std::log10(y);

    // Running average and standard deviation from i to end
    report.Y_avg.resize(report.Y_vec.size());
    report.Y_std.resize(report.Y_vec.size());
    for(const auto &[i, y] : iter::enumerate(report.Y_vec)) {
        report.Y_avg[i] = stat::mean(report.Y_vec, i);
        report.Y_std[i] = stat::stdev(report.Y_vec, i);
    }
    report.Y_sat.resize(report.Y_vec.size());
    for(const auto &[i, s] : iter::enumerate(report.Y_std)) report.Y_sat[i] = static_cast<int>(s < sensitivity); // Below sensitivity --> saturated

    // Since the last element is always zero, we just copy the saturation state of the second to last element.
    if(report.Y_sat.size() > 1) report.Y_sat.back() = report.Y_sat.rbegin()[1];

    // From the end, count how many Y_sat[i] are 1,  before finding a 0.
    for(const auto &[i, y] : iter::enumerate_reverse(report.Y_sat)) {
        report.saturated_point = static_cast<size_t>(i);
        if(y == 0) break;
        report.saturated_count++;
    }
    report.has_computed    = true;
    report.saturated_scale = scale;
    report.has_saturated   = report.saturated_count > 0;
    return report;
}
