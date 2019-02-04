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
#include <IO/class_custom_cout.h>
#include <general/class_tic_toc.h>
#include <general/nmspc_eigsolver_props.h>
#include <sim_parameters/nmspc_sim_settings.h>

class class_superblock;
class class_finite_chain_state;
class class_measurement;
class class_hdf5_file;
class class_table_profiling;
template <typename table_type> class class_hdf5_table;

class class_algorithm_base {

public:
    using Scalar = std::complex<double>;
    class_algorithm_base() = default;
    class_algorithm_base(std::shared_ptr<class_hdf5_file> hdf5_,
                         std::string sim_name_,
                         SimulationType sim_type_);
    void set_profiling_labels ();

    //Storage
    std::shared_ptr<class_hdf5_file>         hdf5;
    std::unique_ptr<class_hdf5_table<class_table_profiling>>  table_profiling;

    std::string sim_name;
    SimulationType sim_type;

    //MPS
    std::shared_ptr<class_superblock>            superblock;
    std::shared_ptr<class_measurement>           measurement;
    std::shared_ptr<class_finite_chain_state>    state;

    //Console
    class_custom_cout ccout;
    //Settings.

    // Common variables
    int    seed       = 1;
    size_t iteration = 0; //In idmrg and itebd: iterations, in fdmrg and xdmrg: full sweeps along the chain.
    size_t step      = 0; //In fdmrg and xdmrg: how many individual moves along the chain.
    size_t num_sites    ;
    long   chi_max      ;
    bool   chi_grow     ;
    size_t print_freq   ;
    size_t store_freq   ;
    long   chi_temp   = 16;
    bool   simulation_has_converged = false;
    bool   simulation_has_to_stop   = false;
    bool   bond_dimension_has_reached_max = false;
    bool   entanglement_has_converged = false;
    bool   entanglement_has_saturated = false;
    size_t variance_mpo_saturated_for = 0;
    bool   variance_mpo_has_converged = false;
    bool   variance_mpo_has_saturated = false;
    size_t variance_ham_saturated_for = 0;
    bool   variance_ham_has_converged = false;
    bool   variance_ham_has_saturated = false;
    size_t variance_mom_saturated_for = 0;
    bool   variance_mom_has_converged = false;
    bool   variance_mom_has_saturated = false;



    //Virtual Functions
    virtual void run()                                          = 0;
    virtual void initialize_constants()                         = 0;
    virtual void print_profiling()                              = 0;
    virtual void print_profiling_sim(class_tic_toc &t_parent)   = 0;
    virtual void store_table_entry_to_file(bool force = false)  = 0;

    //Common functions
    void print_status_update();
    void print_status_full();
    void single_DMRG_step(eigsolver_properties::Ritz ritz = eigsolver_properties::Ritz::SR);
    void store_state_to_file(bool force = false);

    virtual void check_convergence();
    static constexpr double quietNaN = std::numeric_limits<double>::quiet_NaN();
    void check_convergence_variance_mpo(double threshold = quietNaN, double slope_threshold = quietNaN);
    void check_convergence_variance_ham(double threshold = quietNaN, double slope_threshold = quietNaN);
    void check_convergence_variance_mom(double threshold = quietNaN, double slope_threshold = quietNaN);
    void check_convergence_entanglement(double slope_threshold = quietNaN);
    void update_bond_dimension(size_t min_saturation_length = 1);
    void clear_saturation_status();

    void initialize_state(std::string initial_state);
    void reset_chain_mps_to_random_product_state(std::string parity = "none");
    void set_random_fields_in_chain_mpo();
    void compute_observables();
    void enlarge_environment();
    void enlarge_environment(int direction);
    void swap();

    //Functions for finite_chains
    void insert_superblock_to_chain();
    void copy_superblock_mps_to_chain();
    void copy_superblock_mpo_to_chain();
    void copy_superblock_env_to_chain();
    void copy_superblock_to_chain();
    void move_center_point();


    // Profiling
    void store_profiling_to_file_delta(bool force = false);
    void store_profiling_to_file_total(bool force = false);

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

private:
    void check_saturation_using_slope(std::list<bool> &B_vec,
                                      std::list<double> &Y_vec,
                                      std::list<int> &X_vec,
                                      double new_data,
                                      int iter,
                                      int rate,
                                      double tolerance,
                                      double &slope,
                                      bool &has_saturated);

    std::list<bool>   B_mpo_vec; //History of saturation true/false
    std::list<double> V_mpo_vec; //History of variances
    std::list<int>    X_mpo_vec; //History of step numbers
    double V_mpo_slope = 1;

    std::list<bool>   B_ham_vec; //History of saturation true/false
    std::list<double> V_ham_vec;
    std::list<int>    X_ham_vec;
    double V_ham_slope = 1;

    std::list<bool>   B_mom_vec; //History of saturation true/false
    std::list<double> V_mom_vec;
    std::list<int>    X_mom_vec;
    double V_mom_slope = 1;

    std::list<bool>   BS_vec; //History of saturation true/false
    std::list<double> S_vec;
    std::list<int>    XS_vec;
    double S_slope = 1;


};











#endif //DMRG_CLASS_DMRG_BASE_H
