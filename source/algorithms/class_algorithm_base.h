//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_DMRG_BASE_H
#define DMRG_CLASS_DMRG_BASE_H

#include <memory>
#include <map>
#include <vector>
#include <complex>
#include <list>
#include <math/nmspc_eigutils.h>
#include <general/class_tic_toc.h>
#include <simulation/nmspc_settings.h>
#include <simulation/class_simulation_status.h>


template <typename log_type> class class_hdf5_log;
class class_log_profiling;
class class_log_simulation_status;
namespace h5pp{class File;}
namespace spdlog{class logger;}


class class_algorithm_base {
protected:
    std::shared_ptr<spdlog::logger> log;
public:
    using Scalar = std::complex<double>;
    class_algorithm_base() = default;
    class_algorithm_base(std::shared_ptr<h5pp::File> h5ppFile_,
                         std::string sim_name_,
                         SimulationType sim_type_);
    enum class StopReason {CONVERGED, SATURATED, MAX_STEPS} stop_reason;
    void set_profiling_labels ();

    std::shared_ptr<h5pp::File>                                   h5pp_file;
    std::shared_ptr<class_hdf5_log<class_log_profiling>>          log_profiling;
    std::shared_ptr<class_hdf5_log<class_log_simulation_status>>  log_sim_status;

    std::string              sim_name;
    SimulationType           sim_type;
    class_simulation_status  sim_status;




    static constexpr double quietNaN = std::numeric_limits<double>::quiet_NaN();

    //Virtual Functions
    virtual void   run()                                                                                      = 0;
    virtual void   compute_observables()                                                                      = 0;
    virtual void   check_convergence()                                                                        = 0;
    virtual void   write_measurements(bool force = false)                                                     = 0;
    virtual void   write_state       (bool force = false)                                                     = 0;
    virtual void   write_status      (bool force = false)                                                     = 0;
    virtual void   write_logs        (bool force = false)                                                     = 0;
    virtual bool   sim_on()                                                                                   = 0;
    virtual long   chi_max()                                                                                  = 0;
    virtual size_t num_sites()                                                                                = 0;
    virtual size_t write_freq()                                                                               = 0;
    virtual size_t print_freq()                                                                               = 0;
    virtual bool   chi_grow()                                                                                 = 0;
    virtual void   print_status_update()                                                                      = 0;
    virtual void   print_status_full()                                                                        = 0;
    virtual void   reset_to_random_state(const std::string parity = "none", int seed_state = -1)              = 0;
    virtual void   clear_saturation_status()                                                                  = 0;


    //common functions
    void update_bond_dimension();
    void print_profiling();
    double process_memory_in_mb(std::string name);

    // Profiling
    class_tic_toc t_tot;    /*!< Total time */
    class_tic_toc t_pre;    /*!< Preprocessing time */
    class_tic_toc t_run;    /*!< Simulation run time */
    class_tic_toc t_pos;    /*!< Postprocessing time*/
    class_tic_toc t_prt;    /*!< Printing time */
    class_tic_toc t_con;    /*!< Convergence checks time */

protected:
//    using SaturationReport = std::tuple<bool,bool,double,double,int>; //slopes computed, has saturated, rel slope, avgY, check from
    struct SaturationReport {
        bool   has_computed  = false;
        bool   has_saturated = false;
        double slope         = quietNaN;
        double avgY          = quietNaN;
        int    check_from    = -1;
    };

    SaturationReport
    check_saturation_using_slope(std::list<bool> &B_vec,
                                 std::list<double> &Y_vec,
                                 std::list<int> &X_vec,
                                 double new_data,
                                 int iter,
                                 int rate,
                                 double tolerance);

};











#endif //DMRG_CLASS_DMRG_BASE_H
