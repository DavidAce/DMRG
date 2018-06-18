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
class class_finite_chain_sweeper;
class class_measurement;
class class_hdf5_file;
class class_table_profiling;
template <typename table_type> class class_hdf5_table;

class class_base_algorithm {
private:
    void check_convergence_using_slope(std::list<double> &Y_vec,
                                       std::list<int> &X_vec,
                                       double new_data,
                                       int    rate,
                                       double tolerance,
                                       double &slope,
                                       bool &has_converged);

public:
    using Scalar = std::complex<double>;
    class_base_algorithm() = default;
    class_base_algorithm(std::shared_ptr<class_hdf5_file> hdf5_,
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
    std::shared_ptr<class_finite_chain_sweeper>  env_storage;

    //Console
    class_custom_cout ccout;
    //Settings.

    // Common variables
    int    iteration = 0; //In idmrg and itebd: steps, in fdmrg and xdmrg: iteration.
    long   chi_max      ;
    bool   chi_grow     ;
    int    print_freq   ;
    int    store_freq   ;
    int    seed       = 1;
    long   chi_temp   = 4;
    bool   simulation_has_converged = false;
    bool   bond_dimension_has_converged = false;
    bool   entanglement_has_converged = false;
    bool   variance_mpo_has_converged = false;
    bool   variance_ham_has_converged = false;
    bool   variance_mom_has_converged = false;



    //Virtual Functions
    virtual void run()                                          = 0;
    virtual void initialize_constants()                         = 0;
    virtual void print_profiling()                              = 0;
    virtual void print_profiling_sim(class_tic_toc &t_parent)   = 0;
    virtual void store_table_entry_to_file()                    = 0;

    //Common functions
    void print_status_update();
    void print_status_full();
    void single_DMRG_step(long chi_max, Ritz ritz = Ritz::SR);
    void single_TEBD_step(long chi_max);

    virtual void check_convergence_overall();
    void check_convergence_variance_mpo();
    void check_convergence_variance_ham();
    void check_convergence_variance_mom();
    void check_convergence_entanglement();
    void check_convergence_bond_dimension();
    void clear_convergence_checks();

    void initialize_state(std::string initial_state);

    void compute_observables();
    void enlarge_environment();
    void enlarge_environment(int direction);
    void swap();

    //Functions for finite_chains
    void env_storage_insert();
    void env_storage_overwrite_local_MPS();
    void env_storage_overwrite_local_MPO();
    void env_storage_overwrite_local_ENV();
    void env_storage_overwrite_local_ALL();
    void env_storage_move();


    // Profiling
    void store_profiling_to_file();

    class_tic_toc t_tot;
    class_tic_toc t_opt;
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
//    template<typename T>
//    class class_limited_vector : public std::vector<T>{
//    public:
//        void push_back_limited(T val){
//            if (this->size() > 10){
//                this->erase(this->begin());
//            }
//            this->push_back(val);
//        }
//    };
//    class_limited_vector<double> V_mpo_vec;
//    class_limited_vector<double> S_vec;
//    class_limited_vector<int>    X_vec;

    std::list<double> V_mpo_vec;
    std::list<int>    X_mpo_vec;
    double V_mpo_slope;

    std::list<double> V_ham_vec;
    std::list<int>    X_ham_vec;
    double V_ham_slope;

    std::list<double> V_mom_vec;
    std::list<int>    X_mom_vec;
    double V_mom_slope;

    std::list<double> S_vec;
    std::list<int>    X2_vec;
    double S_slope;


};











#endif //DMRG_CLASS_DMRG_BASE_H
