//
// Created by david on 2018-01-18.
//

#ifndef DMRG_CLASS_DMRG_BASE_H
#define DMRG_CLASS_DMRG_BASE_H

#include <memory>
#include <map>
#include <IO/class_custom_cout.h>
#include <general/class_tic_toc.h>
#include <sim_parameters/n_sim_settings.h>

class class_superblock;
class class_hdf5_file;
class class_hdf5_table_buffer;
class class_measurement;

enum class SimulationType{iDMRG,fDMRG, FES_iDMRG, iTEBD, FES_iTEBD, NONE};
class class_base {
public:

    class_base(std::shared_ptr<class_hdf5_file>         hdf5_,
               std::shared_ptr<class_hdf5_table_buffer> table_buffer_,
               std::shared_ptr<class_superblock>        superblock_,
               std::shared_ptr<class_measurement>       observables_);

    class_base();
    SimulationType simtype;
    std::string simulation_name;

    void set_profiling_labels ();
    //Storage
    std::shared_ptr<class_hdf5_file>         hdf5;
    std::shared_ptr<class_hdf5_table_buffer> table_buffer;

    //MPS
    std::shared_ptr<class_superblock>         superblock;
    std::shared_ptr<class_measurement>        observables;

    //Console
    class_custom_cout ccout;

    //Functions
    virtual void run() = 0;
    virtual void store_table_entry() = 0;
    virtual void print_profiling()   = 0;

    virtual void print_status_full()   = 0;  /*!< Print out status of all quantities.*/
    virtual void print_status_update() = 0;  /*!< Print out status of all quantities.*/


    void single_DMRG_step(long chi_max);
    void single_TEBD_step(long chi_max);

    void enlarge_environment();
    void enlarge_environment(int direction);
    void swap();

    // Profiling
    class_tic_toc t_tot;
    class_tic_toc t_eig;
    class_tic_toc t_sim;
    class_tic_toc t_svd;
    class_tic_toc t_env;
    class_tic_toc t_evo;
    class_tic_toc t_udt;
    class_tic_toc t_sto;
    class_tic_toc t_ste;
    class_tic_toc t_obs;
    class_tic_toc t_mps;

};


#endif //DMRG_CLASS_DMRG_BASE_H
