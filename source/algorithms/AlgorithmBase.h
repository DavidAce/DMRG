#pragma once

#include "algorithms/AlgorithmStatus.h"
#include "config/enums.h"
#include <complex>
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
    public:
    using cplx      = std::complex<double>;
    AlgorithmBase() = default;
    AlgorithmBase(std::shared_ptr<h5pp::File> h5ppFile_, AlgorithmType algo_type_);
    std::shared_ptr<h5pp::File> h5file;
    AlgorithmStatus             status;
    static constexpr double     quietNaN = std::numeric_limits<double>::quiet_NaN();

    // Virtual Functions
    virtual void run()                                                                                                                         = 0;
    virtual void run_algorithm()                                                                                                               = 0;
    virtual void run_preprocessing()                                                                                                           = 0;
    virtual void run_postprocessing()                                                                                                          = 0;
    virtual void update_state()                                                                                                                = 0;
    virtual void check_convergence()                                                                                                           = 0;
    virtual void write_to_file(StorageEvent storage_event = StorageEvent::ITERATION, CopyPolicy copy_policy = CopyPolicy::TRY)                = 0;
    virtual void print_status()                                                                                                                = 0;
    virtual void print_status_full()                                                                                                           = 0;
    virtual void clear_convergence_status()                                                                                                    = 0;
    virtual void update_variance_max_digits(std::optional<double> energy = std::nullopt)                                                       = 0;
    virtual void update_bond_dimension_limit()                                                                                                 = 0;
    virtual void update_truncation_error_limit()                                                                                               = 0;

    // common functions
    void copy_from_tmp(StorageEvent storage_event = StorageEvent::ITERATION, CopyPolicy copy_policy = CopyPolicy::TRY);
    void init_bond_dimension_limits();
    void init_truncation_error_limits();

    protected:
    enum class SaturationScale { lin, log };
    struct SaturationReport {
        bool                has_computed    = false;
        bool                has_saturated   = false;
        size_t              saturated_count = 0;
        size_t              saturated_point = 0;
        SaturationScale     saturated_scale = SaturationScale::lin;
        std::vector<double> Y_vec;     // The values used to measure saturation. This is the log of given data when log is on
        std::vector<double> Y_min;     // The cumulative minimum of Y_vec
        std::vector<double> Y_max;     // The cumulative maximum of Y_vec
        std::vector<double> Y_vec_std; // The "moving start" standard deviation of Y_vec from [i:end]
        std::vector<double> Y_min_std; // The "moving start" standard deviation of Y_min from [i:end]
        std::vector<double> Y_max_std; // The "moving start" standard deviation of Y_max from [i:end]
        std::vector<double> Y_mov_ste; // The "moving window" standard error of 25% width.
        std::vector<int>    Y_sat;     // Flags that tell if Y_vec has saturated at that index

        [[nodiscard]] constexpr std::string_view get_saturated_scale() noexcept {
            switch(saturated_scale) {
                case SaturationScale::lin: return "lin";
                case SaturationScale::log: return "log";
            }
            return "error";
        }
    };
    size_t           count_convergence(const std::vector<double> &Y_vec, double threshold, size_t start_idx = 0);
    SaturationReport check_saturation(const std::vector<double> &Y_vec, double sensitivity, SaturationScale scale);
};
