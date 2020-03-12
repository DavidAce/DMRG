//
// Created by david on 2018-01-18.
//

#pragma once

#include <complex>
#include <general/class_tic_toc.h>
#include <list>
#include <map>
#include <memory>
#include <simulation/class_simulation_status.h>
#include <simulation/enums.h>
#include <vector>

enum class StopReason { SUCCEEDED, SATURATED, MAX_ITERS, MAX_RESET, RANDOMIZE,NONE };

namespace h5pp {
    class File;
}
namespace spdlog {
    class logger;
}

class class_algorithm_base {
    protected:
    std::shared_ptr<spdlog::logger> log;

    public:
    using Scalar           = std::complex<double>;
    class_algorithm_base() = default;
    class_algorithm_base(std::shared_ptr<h5pp::File> h5ppFile_, std::string sim_name_, SimulationType sim_type_);
    std::shared_ptr<h5pp::File> h5pp_file;
    StopReason                  stop_reason = StopReason::NONE;
    std::string                 sim_tag;
    std::string                 sim_name;
    SimulationType              sim_type;
    class_simulation_status     sim_status;

    static constexpr double quietNaN = std::numeric_limits<double>::quiet_NaN();

    // Virtual Functions
    virtual void   run()                                                                                         = 0;
    virtual void   check_convergence()                                                                           = 0;
    virtual void   write_state(bool result = false)                                                              = 0;
    virtual void   write_measurements(bool result = false)                                                       = 0;
    virtual void   write_sim_status(bool result = false)                                                         = 0;
    virtual void   write_profiling(bool result = false)                                                          = 0;
    virtual void   copy_from_tmp(bool result = false)                                                            = 0;
    virtual bool   sim_on()                                                                                      = 0;
    virtual long   chi_max()                                                                                     = 0;
    virtual size_t num_sites()                                                                                   = 0;
    virtual size_t write_freq()                                                                                  = 0;
    virtual size_t print_freq()                                                                                  = 0;
    virtual bool   chi_grow()                                                                                    = 0;
    virtual long   chi_init()                                                                                    = 0;
    virtual void   print_status_update()                                                                         = 0;
    virtual void   print_status_full()                                                                           = 0;
    virtual void   reset_to_random_product_state(const std::string &parity = "random")                           = 0;
    virtual void   reset_to_random_current_state(std::optional<double> chi_lim = std::nullopt)                   = 0;
    virtual void   reset_to_initial_state()                                                                      = 0;
    virtual void   clear_saturation_status()                                                                     = 0;
    virtual void   update_truncation_limit()                                                                     = 0;
    virtual void   update_bond_dimension_limit(std::optional<long> max_bond_dim = std::nullopt)                  = 0;

    // common functions
    void   print_profiling();
    double process_memory_in_mb(std::string name);

    protected:
    //    using SaturationReport = std::tuple<bool,bool,double,double,int>; //slopes computed, has saturated, rel slope, avgY, check from
    struct SaturationReport {
        bool   has_computed = false;
        size_t check_from   = 0;
        double slope;
        double avgY;
    };
//    struct SaturationReport2 {
//        bool                has_computed  = false;
//        bool                has_saturated = false;
//        size_t              saturated_for = 0;
//        std::vector<double> slopes;
//        std::vector<double> avgY;
//    };

    SaturationReport  check_saturation_using_slope(std::list<double> &Y_vec, std::list<int> &X_vec, double new_data, int iter, int rate, double tolerance);
//    SaturationReport2 check_saturation_using_slope2(std::list<double> &Y_vec, std::list<int> &X_vec, double new_data, int iter, int rate, double tolerance);
};
