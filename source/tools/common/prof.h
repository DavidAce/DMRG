#pragma once

#include <general/class_tic_toc.h>
#include <memory>
#include <optional>
class class_tic_toc;
enum class AlgorithmType;
namespace tools::common::profile{
    // Profiling
    inline std::unique_ptr<class_tic_toc> t_tot;             // + Total Time
    inline std::unique_ptr<class_tic_toc> t_pre;             // |- Preprocessing
    inline std::unique_ptr<class_tic_toc> t_pos;             // |- Postprocessing
    inline std::unique_ptr<class_tic_toc> t_sim;             // |+ Simulation
    inline std::unique_ptr<class_tic_toc> t_con;             //  |- Convergence checks
    inline std::unique_ptr<class_tic_toc> t_eig;             //  |- Eig. decomp.
    inline std::unique_ptr<class_tic_toc> t_svd;             //  |- Svd. decomp.
    inline std::unique_ptr<class_tic_toc> t_evo;             //  |- Time evolution
    inline std::unique_ptr<class_tic_toc> t_env;             //  |- Environment upd.
    inline std::unique_ptr<class_tic_toc> t_ent;             //  |- Entanglement entropy
    inline std::unique_ptr<class_tic_toc> t_ene;             //  |- Energy
    inline std::unique_ptr<class_tic_toc> t_var;             //  |- Variance
    inline std::unique_ptr<class_tic_toc> t_prj;             //  |- Projections
    inline std::unique_ptr<class_tic_toc> t_chk;             //  |- Checks
    inline std::unique_ptr<class_tic_toc> t_hdf;             //  |- h5pp storage
    inline std::unique_ptr<class_tic_toc> t_ene_ham;         //  |- Energy (HAM)
    inline std::unique_ptr<class_tic_toc> t_ene_mom;         //  |- Energy (MOM)
    inline std::unique_ptr<class_tic_toc> t_var_ham;         //  |- Variance (HAM)
    inline std::unique_ptr<class_tic_toc> t_var_mom;         //  |- Variance (MOM)
    inline std::unique_ptr<class_tic_toc> t_mps;             //  |- Multisite-MPS
    inline std::unique_ptr<class_tic_toc> t_mpo;             //  |- Multisite-MPO
    inline std::unique_ptr<class_tic_toc> t_opt;             //  |+ Optimization (xdmrg)
    inline std::unique_ptr<class_tic_toc> t_opt_dir;         //  ||+ Direct
    inline std::unique_ptr<class_tic_toc> t_opt_dir_bfgs;    //  |||+ L-BFGS steps
    inline std::unique_ptr<class_tic_toc> t_opt_dir_vH2;     //  || |- vH2
    inline std::unique_ptr<class_tic_toc> t_opt_dir_vH2v;    //  || |- vH2v
    inline std::unique_ptr<class_tic_toc> t_opt_dir_vH;      //  || |- vH
    inline std::unique_ptr<class_tic_toc> t_opt_dir_vHv;     //  || |- vHv
    inline std::unique_ptr<class_tic_toc> t_opt_sub;         //  ||+ Subspace
    inline std::unique_ptr<class_tic_toc> t_opt_sub_ham;     //  | |- Hamiltonian Matrix
    inline std::unique_ptr<class_tic_toc> t_opt_sub_hsq;     //  | |- Hamiltonian Matrix²
    inline std::unique_ptr<class_tic_toc> t_opt_sub_lu;      //  | |- LU decomposition
    inline std::unique_ptr<class_tic_toc> t_opt_sub_eig;     //  | |- Eigenvalue decomp
    inline std::unique_ptr<class_tic_toc> t_opt_sub_bfgs;    //  | |+ L-BFGS steps
    inline std::unique_ptr<class_tic_toc> t_opt_sub_vH2;     //  |  |- vH2
    inline std::unique_ptr<class_tic_toc> t_opt_sub_vH2v;    //  |  |- vH2v
    inline std::unique_ptr<class_tic_toc> t_opt_sub_vH;      //  |  |- vH
    inline std::unique_ptr<class_tic_toc> t_opt_sub_vHv;     //  |  |- vHv


    extern void print_profiling();
    extern void print_profiling(std::optional<AlgorithmType> algo_type);
    extern void init_profiling();
    extern void reset_profiling();

    extern void print_mem_usage();

    extern double mem_usage_in_mb(std::string_view name);
    extern double mem_rss_in_mb();
    extern double mem_hwm_in_mb();
    extern double mem_vm_in_mb();

}