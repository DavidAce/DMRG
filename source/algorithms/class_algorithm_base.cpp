//
// Created by david on 2018-01-18.
//

#include "class_algorithm_base.h"
#include <complex>
#include <config/nmspc_settings.h>
#include <h5pp/h5pp.h>
#include <math/num.h>
#include <math/stat.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

using Scalar = class_algorithm_base::Scalar;

class_algorithm_base::class_algorithm_base(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type_)
    : h5pp_file(std::move(h5ppFile_)), algo_type(algo_type_) {
    algo_name = enum2str(algo_type_);
    tools::common::profile::set_default_prof(algo_type);
    tools::common::profile::init_profiling();
    tools::log->set_error_handler([](const std::string &msg) { throw std::runtime_error(msg); });
    tools::log = tools::Logger::setLogger(std::string(enum2str(algo_type)), settings::console::verbosity, settings::console::timestamp);
    tools::log->trace("Constructing class_algorithm_base");
}

void class_algorithm_base::copy_from_tmp(StorageReason storage_reason, std::optional<CopyPolicy> copy_policy) {
    if(not h5pp_file) return;
    if(not settings::output::use_temp_dir) return;
    if(not copy_policy) return copy_from_tmp(storage_reason, CopyPolicy::TRY);
    if(copy_policy == CopyPolicy::OFF) return;

    // Check if we already copied the file this iteration and step
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    auto                                                                  save_point = std::make_pair(status.iter, status.step);

    if(copy_policy == CopyPolicy::TRY) {
        if(save_log[h5pp_file->getFilePath()] == save_point) return;
        switch(storage_reason) {
            case StorageReason::SAVEPOINT:
            case StorageReason::CHECKPOINT:
                if(num::mod(status.iter, settings::output::copy_from_temp_freq) != 0) return; // Check that we write according to the frequency given
            case StorageReason::FINISHED:
            case StorageReason::CHI_UPDATE:
            case StorageReason::PROJ_STATE:
            case StorageReason::INIT_STATE:
            case StorageReason::EMIN_STATE:
            case StorageReason::EMAX_STATE:
            case StorageReason::MODEL: break;
        }
        tools::common::io::h5tmp::copy_from_tmp(h5pp_file->getFilePath());
    } else if(copy_policy == CopyPolicy::FORCE)
        tools::common::io::h5tmp::copy_from_tmp(h5pp_file->getFilePath());

    save_log[h5pp_file->getFilePath()] = save_point;
}

void class_algorithm_base::init_bond_dimension_limits() {
    status.chi_lim_init = cfg_chi_lim_init();
    status.chi_lim_max  = cfg_chi_lim_max();
    if(cfg_chi_lim_grow())
        status.chi_lim = cfg_chi_lim_init();
    else
        status.chi_lim = cfg_chi_lim_max();

    // Sanity check
    if(status.chi_lim == 0) throw std::runtime_error(fmt::format("Bond dimension limit invalid: {}", status.chi_lim));
}

/*! \brief Checks convergence based on slope.
 * We want to check once every "rate" steps. First, check the sim_state.iteration number when you last measured.
 * If the measurement happened less than rate iterations ago, return.
 * Otherwise, compute the slope of the last 25% of the measurements that have been made.
 * The slope here is defined as the relative slope, i.e. \f$ \frac{1}{ \langle y\rangle} * \frac{dy}{dx} \f$.
 */
class_algorithm_base::SaturationReport class_algorithm_base::check_saturation_using_slope(std::vector<double> &Y_vec, std::vector<size_t> &X_vec,
                                                                                          double new_data, size_t iter, size_t rate, double tolerance) {
    SaturationReport report;
    size_t           last_measurement = X_vec.empty() ? 0 : X_vec.back();
    if(iter < rate + last_measurement) { return report; }

    // It's time to check. Insert current numbers
    Y_vec.push_back(new_data);
    X_vec.push_back(iter);
    size_t min_data_points = 2;
    if(Y_vec.size() < min_data_points) { return report; }
    size_t start_point = 0;
    double band_size   = 2.0 + 2.0 * tolerance; // Between 2 and  4 standard deviations away

    // Consider Y_vec vs X_vec: a noisy signal decaying in the shape of a hockey-club, say.
    // We want to identify the point at which the signal stabilizes. We use the fact that the
    // standard deviation is high if it includes parts of the non-stable signal, and low if
    // it includes only the stable part.
    // Here we monitor the standard deviation of the signal between [some_point, X_vec.end()],
    // and move "some_point" towards the end. If the standard deviation goes below a certain
    // threshold, we've found the stabilization point.
    auto recent_point       = static_cast<size_t>(std::floor(0.75 * static_cast<double>(Y_vec.size())));
    recent_point            = std::min(Y_vec.size() - min_data_points, recent_point);
    double recent_point_std = stat::stdev(Y_vec, recent_point); // Computes the standard dev of Y_vec from recent_point to end
    for(size_t some_point = 0; some_point < Y_vec.size(); some_point++) {
        double some_point_std = stat::stdev(Y_vec, some_point); // Computes the standard dev of Y_vec from some_point to end
        if(some_point_std < band_size * recent_point_std and start_point == 0) {
            start_point = some_point;
            break;
        }
    }
    // Scale the slope so that it can be interpreted as change in percent, just as the tolerance.
    double avgY  = stat::mean(Y_vec, start_point);
    auto [slope,res] = stat::slope(X_vec, Y_vec, start_point);
    slope = std::abs(slope) / avgY * 100 / std::sqrt(Y_vec.size() - start_point); // TODO: Is dividing by sqrt(elems) reasonable?
    slope        = std::isnan(slope) ? 0.0 : slope;
    report.slope = slope;
    report.check_from   = start_point;
    report.avgY         = avgY;
    report.has_computed = true;
    return report;
}

void class_algorithm_base::print_profiling_lap() {
    if(not settings::profiling::extra) return;
    tools::common::profile::print_profiling_laps();
}