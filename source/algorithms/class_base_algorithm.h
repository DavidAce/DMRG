//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_DMRG_BASE_H
#define DMRG_CLASS_DMRG_BASE_H

#include <memory>
#include <map>
#include <list>
#include <IO/class_custom_cout.h>
#include <general/class_tic_toc.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <vector>
#include <unsupported/Eigen/CXX11/Tensor>
#include <mps_routines/class_finite_chain_storage.h>

class class_superblock;
class class_hdf5_file;
class class_hdf5_table_buffer;
class class_measurement;




class class_base_algorithm {
public:
    using Scalar = std::complex<double>;
    class_base_algorithm() = default;
    class_base_algorithm(std::shared_ptr<class_hdf5_file> hdf5_,
                         std::string sim_name_,
                         std::string table_name_,
                         SimulationType sim_type_);
    void set_profiling_labels ();

    //Storage
    std::shared_ptr<class_hdf5_file>         hdf5;
    std::shared_ptr<class_hdf5_table_buffer> table_buffer;

    std::string sim_name;
    std::string table_name;
    SimulationType sim_type;

    //MPS
    std::shared_ptr<class_superblock>            superblock;
    std::shared_ptr<class_measurement>           measurement;
    std::shared_ptr<class_finite_chain_storage>  env_storage;

    //Console
    class_custom_cout ccout;
    //Settings.

    // Static variables
    long   chi_max      ;
    long   chi_min      ;
    long   chi_num      ;
    bool   chi_grow     ;
    int    max_length   ;
    int    max_sweeps   ;
    int    print_freq   ;
    int    store_freq   ;
    int    max_steps    ;
    double delta_t0     ;
    double delta_tmin   ;
    int    suzuki_order ;

    // Variables for excited state DMRG
    int    seed         ;
    double r_strength   ;  //Randomness strength for the random field.


    // Variables that monitor the simulation
    int    sweeps    = 0;
    int    iteration = 0;
    int    position  = 0;
    int    direction = 1;
    double phys_time = 0;
    double delta_t   = 0; //Make sure this one gets initialized to delta_t0!
    long   chi_temp  = 2;
    bool   simulation_has_converged = false;




    //Virtual Functions
    virtual void run() = 0;
    virtual void print_profiling()   = 0;
    virtual void print_profiling_sim(class_tic_toc &t_parent) = 0;

    //Common functions
    void store_table_entry();
    void print_status_update();
    void print_status_full();
    void single_DMRG_step(long chi_max);
    void single_TEBD_step(long chi_max);

    void update_chi();
    bool entropy_has_converged();
    void initialize_constants();
    void initialize_state(std::string initial_state);

    void compute_observables();
    int  enlarge_environment();
    int  enlarge_environment(int direction);
    void swap();

    //Functions for finite_chains
    int  env_storage_insert();
    int  env_storage_insert_edges();
    void env_storage_overwrite_MPS();
    int  env_storage_move();


    // Profiling
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
    class_tic_toc t_chi;
private:
    template<typename T>
    class class_limited_vector : public std::vector<T>{
    public:
        void push_back_limited(T val){
            if (this->size() > 10){
                this->erase(this->begin());
            }
            this->push_back(val);
        }
    };
    class_limited_vector<double> S_vec;
    class_limited_vector<int>    X_vec;

};


#endif //DMRG_CLASS_DMRG_BASE_H
