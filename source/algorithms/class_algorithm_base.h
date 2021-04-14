//
// Created by david on 2018-01-18.
//

#pragma once

#include <algorithms/class_algorithm_status.h>
#include <complex>
#include <config/enums.h>
#include <memory>
#include <optional>
#include <vector>

namespace h5pp {
    class File;
}
namespace spdlog {
    class logger;
}

class class_algorithm_base {
    protected:
    public:
    using Scalar           = std::complex<double>;
    class_algorithm_base() = default;
    class_algorithm_base(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type_);
    std::shared_ptr<h5pp::File> h5pp_file;
    class_algorithm_status      status;
    StopReason                  stop_reason = StopReason::NONE;
    AlgorithmType               algo_type;
    std::string                 algo_name;
    static constexpr double     quietNaN = std::numeric_limits<double>::quiet_NaN();

    // Virtual Functions
    virtual void   run()                                                                                                                         = 0;
    virtual void   check_convergence()                                                                                                           = 0;
    virtual void   write_to_file(StorageReason storage_reason = StorageReason::CHECKPOINT, std::optional<CopyPolicy> copy_policy = std::nullopt) = 0;
    virtual bool   cfg_algorithm_is_on()                                                                                                         = 0;
    virtual size_t cfg_print_freq()                                                                                                              = 0;
    virtual long   cfg_chi_lim_max()                                                                                                             = 0;
    virtual bool   cfg_chi_lim_grow()                                                                                                            = 0;
    virtual long   cfg_chi_lim_init()                                                                                                            = 0;
    virtual void   print_status_update()                                                                                                         = 0;
    virtual void   print_status_full()                                                                                                           = 0;
    virtual void   clear_convergence_status()                                                                                                    = 0;
    virtual void   update_variance_max_digits(std::optional<double> energy = std::nullopt)                                                       = 0;
    virtual void   update_bond_dimension_limit(std::optional<long> max_bond_dim = std::nullopt)                                                  = 0;

    // common functions
    void copy_from_tmp(StorageReason storage_reason = StorageReason::SAVEPOINT, std::optional<CopyPolicy> copy_policy = std::nullopt);
    void init_bond_dimension_limits();
    void update_variance_convergence_threshold();
    void print_profiling_lap();

    protected:
    //    using SaturationReport = std::tuple<bool,bool,double,double,int>; //slopes computed, has saturated, rel slope, avgY, check from
    struct SaturationReport {
        bool                has_computed = false;
        bool                has_saturated = false;
        size_t              saturated_count = 0;
        size_t              saturated_point = 0;
        double              Y_avg; // Average from the saturation point onward
        std::vector<double> Y_vec; // The values used to gauge saturation
        std::vector<double> Y_log; // Normalized values to check saturation. Let y = -log10(Y_vec). Then Y_log = y/y.back()
        std::vector<double> Y_std; // The "moving" standard deviation of Y_log. (std from x -> end, moving x towards end)
        std::vector<double> Y_ste; // The "moving" standard error of Y_log. (std from x -> end, moving x towards end)
    };
    size_t count_convergence(const std::vector<double> & Y_vec, double threshold);
    SaturationReport check_saturation(const std::vector<double> &Y_vec, double sensitivity);
};
