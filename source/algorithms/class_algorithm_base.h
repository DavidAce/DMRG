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
class class_hdf5_file;
template <typename table_type> class class_hdf5_table;
class class_table_profiling;
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

    std::shared_ptr<h5pp::File>                               h5pp_file;
    std::unique_ptr<class_hdf5_table<class_table_profiling>>  table_profiling;

    std::string             sim_name;
    SimulationType          sim_type;
    class_simulation_status  sim_status;




    static constexpr double quietNaN = std::numeric_limits<double>::quiet_NaN();

    //Virtual Functions
    virtual void   run()                                                                                      = 0;
    virtual void   compute_observables()                                                                      = 0;
    virtual void   check_convergence()                                                                        = 0;
    virtual void   store_state_and_measurements_to_file(bool force = false)                                   = 0;
//    virtual void   store_table_entry_progress(bool force = false)                                             = 0;
    virtual bool   sim_on()                                                                                   = 0;
    virtual long   chi_max()                                                                                  = 0;
    virtual size_t num_sites()                                                                                = 0;
    virtual size_t store_freq()                                                                               = 0;
    virtual size_t print_freq()                                                                               = 0;
    virtual bool   chi_grow()                                                                                 = 0;
    virtual void   print_status_update()                                                                      = 0;
    virtual void   print_status_full()                                                                        = 0;
    virtual void   print_profiling()                                                                          = 0;
    virtual void   print_profiling_sim(class_tic_toc &t_parent)                                               = 0;
    virtual void   reset_to_random_state(const std::string parity)                                            = 0;
    virtual void   clear_saturation_status()                                                                  = 0;


    //common functions
    void store_algorithm_state_to_file();

    void update_bond_dimension();



    //Functions for finite_chains
    void insert_superblock_to_chain();
    void copy_superblock_mps_to_chain();
    void copy_superblock_mpo_to_chain();
    void copy_superblock_env_to_chain();
    void copy_superblock_to_chain();
    void move_center_point();

    double process_memory_in_mb(std::string name);

    // Profiling
    void store_profiling_deltas(bool force = false);
    void store_profiling_totals(bool force = false);

    class_tic_toc t_tot;
    class_tic_toc t_opt;
    class_tic_toc t_eig;
    class_tic_toc t_ham;
    class_tic_toc t_sim;
    class_tic_toc t_svd;
    class_tic_toc t_env;
    class_tic_toc t_evo;
    class_tic_toc t_udt;
    class_tic_toc t_sto;
    class_tic_toc t_ste;
    class_tic_toc t_prt;
    class_tic_toc t_obs;
    class_tic_toc t_mps;
    class_tic_toc t_con;

protected:
    void check_saturation_using_slope(std::list<bool> &B_vec,
                                      std::list<double> &Y_vec,
                                      std::list<int> &X_vec,
                                      double new_data,
                                      int iter,
                                      int rate,
                                      double tolerance,
                                      double &slope,
                                      bool &has_saturated);




};











#endif //DMRG_CLASS_DMRG_BASE_H
