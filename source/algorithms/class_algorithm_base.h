//
// Created by david on 2018-01-18.
//

#pragma once

#include <algorithms/class_algorithm_status.h>
#include <complex>
#include <config/enums.h>
#include <general/class_tic_toc.h>
#include <list>
#include <map>
#include <memory>
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
    std::string                 state_name;
    static constexpr double     quietNaN = std::numeric_limits<double>::quiet_NaN();

    // Virtual Functions
    virtual void   run()                                                                                                                           = 0;
    virtual void   check_convergence()                                                                                                             = 0;
    virtual void   write_to_file(StorageReason storage_reason = StorageReason::CHECKPOINT)                                                         = 0;
    virtual void   copy_from_tmp(StorageReason storage_reason = StorageReason::CHECKPOINT)                                                         = 0;
    virtual bool   cfg_algorithm_is_on()                                                                                                                       = 0;
    virtual size_t cfg_print_freq()                                                                                                                    = 0;
    virtual long   cfg_chi_lim_max()                                                                                                                   = 0;
    virtual bool   cfg_chi_lim_grow()                                                                                                                  = 0;
    virtual long   cfg_chi_lim_init()                                                                                                                  = 0;
    virtual void   print_status_update()                                                                                                           = 0;
    virtual void   print_status_full()                                                                                                             = 0;
    virtual void randomize_into_product_state(ResetReason reason, std::optional<std::string> sector = std::nullopt,
                                                 std::optional<long> bitfield = std::nullopt, std::optional<bool> use_eigenspinors = std::nullopt) = 0;
    virtual void   randomize_from_current_state(std::optional<std::vector<std::string>> pauli_strings = std::nullopt, std::optional<std::string> sector = std::nullopt,
                                         std::optional<long> chi_lim = std::nullopt, std::optional<double> svd_threshold = std::nullopt)           = 0;
    virtual void   clear_convergence_status()                                                                                                         = 0;
    virtual void update_truncation_limit()                                                                                                         = 0;
    virtual void update_bond_dimension_limit(std::optional<long> max_bond_dim = std::nullopt)                                                      = 0;

    // common functions
    void   print_profiling();
    void   init_bond_dimension_limits();
    double process_memory_in_mb(std::string name);

    protected:
    //    using SaturationReport = std::tuple<bool,bool,double,double,int>; //slopes computed, has saturated, rel slope, avgY, check from
    struct SaturationReport {
        bool   has_computed = false;
        size_t check_from   = 0;
        double slope;
        double avgY;
    };

    SaturationReport check_saturation_using_slope(std::list<double> &Y_vec, std::list<size_t> &X_vec, double new_data, size_t iter, size_t rate,
                                                  double tolerance);
    //    SaturationReport2 check_saturation_using_slope2(std::list<double> &Y_vec, std::list<int> &X_vec, double new_data, int iter, int rate, double
    //    tolerance);

};
