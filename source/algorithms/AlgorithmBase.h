#pragma once

#include <algorithms/AlgorithmStatus.h>
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

class AlgorithmBase {
    protected:
    public:
    using Scalar    = std::complex<double>;
    AlgorithmBase() = default;
    AlgorithmBase(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type_);
    std::shared_ptr<h5pp::File> h5pp_file;
    AlgorithmStatus             status;
    static constexpr double     quietNaN = std::numeric_limits<double>::quiet_NaN();

    // Virtual Functions
    virtual void run()                                                                                                                         = 0;
    virtual void check_convergence()                                                                                                           = 0;
    virtual void write_to_file(StorageReason storage_reason = StorageReason::CHECKPOINT, std::optional<CopyPolicy> copy_policy = std::nullopt) = 0;
    virtual void print_status_update()                                                                                                         = 0;
    virtual void print_status_full()                                                                                                           = 0;
    virtual void clear_convergence_status()                                                                                                    = 0;
    virtual void update_variance_max_digits(std::optional<double> energy = std::nullopt)                                                       = 0;
    virtual void update_bond_dimension_limit()                                                  = 0;

    // common functions
    void copy_from_tmp(StorageReason storage_reason = StorageReason::SAVEPOINT, std::optional<CopyPolicy> copy_policy = std::nullopt);
    void init_bond_dimension_limits();

    protected:
    //    using SaturationReport = std::tuple<bool,bool,double,double,int>; //slopes computed, has saturated, rel slope, avgY, check from
    struct SaturationReport {
        bool                has_computed    = false;
        bool                has_saturated   = false;
        size_t              saturated_count = 0;
        size_t              saturated_point = 0;
        std::vector<double> Y_vec; // The values used to gauge saturation
        std::vector<double> Y_avg; // Running average from [i:end]
        std::vector<double> Y_std; // The "moving" standard deviation of Y_avg from [i:end]
        std::vector<double> Y_stn; // The "normalized" scale invariant version of Y_std
        //        std::vector<double> Y_log; // Normalized values to check saturation. Let y = -log10(Y_vec). Then Y_log = y/y.back()
        //        std::vector<double> Y_ste; // The "moving" standard error of Y_log. (std from x -> end, moving x towards end)
        std::vector<double> Y_slp; // The "moving" slope of Y_log. (slope from x -> end, moving x towards end)
    };
    size_t           count_convergence(const std::vector<double> &Y_vec, double threshold, size_t start_idx = 0);
    SaturationReport check_saturation(const std::vector<double> &Y_vec, double sensitivity);
};
