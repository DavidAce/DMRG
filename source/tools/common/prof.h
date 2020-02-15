#pragma once

#include <general/class_tic_toc.h>
#include <memory>
class class_tic_toc;
namespace tools::common::profile{
    // Profiling
    inline std::unique_ptr<class_tic_toc> t_tot;        /*!< Total  */
    inline std::unique_ptr<class_tic_toc> t_pre;        /*!< Preprocessing  */
    inline std::unique_ptr<class_tic_toc> t_pos;        /*!< Postprocessing */
    inline std::unique_ptr<class_tic_toc> t_sim;        /*!< Simulation runtime  */
    inline std::unique_ptr<class_tic_toc> t_con;        /*!< Convergence checks */
    inline std::unique_ptr<class_tic_toc> t_eig;        /*!< Eigenvalue decomposition */
    inline std::unique_ptr<class_tic_toc> t_svd;        /*!< SVD decomposition */
    inline std::unique_ptr<class_tic_toc> t_opt;        /*!< Optimization, i.e. Ceres */
    inline std::unique_ptr<class_tic_toc> t_evo;        /*!< Time evolution */
    inline std::unique_ptr<class_tic_toc> t_env;        /*!< Update environments */
    inline std::unique_ptr<class_tic_toc> t_ent;        /*!< Entanglement entropy */
    inline std::unique_ptr<class_tic_toc> t_ene;        /*!< Energy */
    inline std::unique_ptr<class_tic_toc> t_var;        /*!< Variance */
    inline std::unique_ptr<class_tic_toc> t_prj;        /*!< Projections */
    inline std::unique_ptr<class_tic_toc> t_chk;        /*!< Integrity checks/debugging */
    inline std::unique_ptr<class_tic_toc> t_hdf;        /*!< hdf5 -- writing to file */
    inline std::unique_ptr<class_tic_toc> t_ene_ham;
    inline std::unique_ptr<class_tic_toc> t_ene_mom;
    inline std::unique_ptr<class_tic_toc> t_var_ham;
    inline std::unique_ptr<class_tic_toc> t_var_mom;

    inline std::unique_ptr<class_tic_toc> t_ham ;
    inline std::unique_ptr<class_tic_toc> t_vH2v;
    inline std::unique_ptr<class_tic_toc> t_vHv ;
    inline std::unique_ptr<class_tic_toc> t_vH2 ;
    inline std::unique_ptr<class_tic_toc> t_vH  ;
    inline std::unique_ptr<class_tic_toc> t_op  ;


    extern void print_profiling();
    extern void init_profiling();
}