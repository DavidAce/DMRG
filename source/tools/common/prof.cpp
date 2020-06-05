//
// Created by david on 2019-06-08.
//

#include <config/enums.h>
#include <config/nmspc_settings.h>
#include <fstream>
#include <sstream>
#include <tools/common/log.h>
#include <tools/common/prof.h>

void tools::common::profile::print_profiling() { print_profiling(std::nullopt); }

void tools::common::profile::print_profiling(std::optional<AlgorithmType> algo_type) {
    if(t_tot == nullptr) return;
    if(settings::profiling::on) {
        auto t_tot_percent = 100.0 / std::max(1.0, t_tot->get_measured_time());
        auto t_sim_percent = 100.0 / std::max(1.0, t_sim->get_measured_time());

        /* clang-format off */
        tools::log->info("{:<30}{:>10.3f} s"                       , t_tot->get_name(), t_tot->get_measured_time());
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_pre->get_name(), t_pre->get_measured_time()     , t_pre->get_measured_time()     * t_tot_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_pos->get_name(), t_pos->get_measured_time()     , t_pos->get_measured_time()     * t_tot_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_sim->get_name(), t_sim->get_measured_time()     , t_sim->get_measured_time()     * t_tot_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_con->get_name(), t_con->get_measured_time()     , t_con->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_eig->get_name(), t_eig->get_measured_time()     , t_eig->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_svd->get_name(), t_svd->get_measured_time()     , t_svd->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_evo->get_name(), t_evo->get_measured_time()     , t_evo->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_env->get_name(), t_env->get_measured_time()     , t_env->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_ent->get_name(), t_ent->get_measured_time()     , t_ent->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_ene->get_name(), t_ene->get_measured_time()     , t_ene->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_var->get_name(), t_var->get_measured_time()     , t_var->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_prj->get_name(), t_prj->get_measured_time()     , t_prj->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_chk->get_name(), t_chk->get_measured_time()     , t_chk->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_hdf->get_name(), t_hdf->get_measured_time()     , t_hdf->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_ene->get_name(), t_ene_ham->get_measured_time() , t_ene_ham->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_ene->get_name(), t_ene_mom->get_measured_time() , t_ene_mom->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_var->get_name(), t_var_ham->get_measured_time() , t_var_ham->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_var->get_name(), t_var_mom->get_measured_time() , t_var_mom->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_mps->get_name(), t_mps->get_measured_time()     , t_mps->get_measured_time()     * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_mpo->get_name(), t_mpo->get_measured_time()     , t_mpo->get_measured_time()     * t_sim_percent);
        if(algo_type and algo_type.value() != AlgorithmType::xDMRG) return;
        auto t_opt_percent = 100.0/std::max(1.0,t_opt->get_measured_time());
        auto t_opt_dir_percent = 100.0/std::max(1.0,t_opt_dir->get_measured_time());
        auto t_opt_sub_percent = 100.0/std::max(1.0,t_opt_sub->get_measured_time());
        auto t_opt_dir_bfgs_percent = 100.0/std::max(1.0,t_opt_dir_bfgs->get_measured_time());
        auto t_opt_sub_bfgs_percent = 100.0/std::max(1.0,t_opt_sub_bfgs->get_measured_time());
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt->get_name()          ,t_opt->get_measured_time()            , t_opt->get_measured_time()              * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_dir->get_name()      ,t_opt_dir->get_measured_time()        , t_opt_dir->get_measured_time()          * t_opt_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_dir_bfgs->get_name() ,t_opt_dir_bfgs->get_measured_time()   , t_opt_dir_bfgs->get_measured_time()     * t_opt_dir_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_dir_vH2->get_name()  ,t_opt_dir_vH2->get_measured_time()    , t_opt_dir_vH2->get_measured_time()      * t_opt_dir_bfgs_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_dir_vH2v->get_name() ,t_opt_dir_vH2v->get_measured_time()   , t_opt_dir_vH2v->get_measured_time()     * t_opt_dir_bfgs_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_dir_vH->get_name()   ,t_opt_dir_vH->get_measured_time()     , t_opt_dir_vH->get_measured_time()       * t_opt_dir_bfgs_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_dir_vHv->get_name()  ,t_opt_dir_vHv->get_measured_time()    , t_opt_dir_vHv->get_measured_time()      * t_opt_dir_bfgs_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub->get_name()      ,t_opt_sub->get_measured_time()        , t_opt_sub->get_measured_time()          * t_opt_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub_ham->get_name()  ,t_opt_sub_ham->get_measured_time()    , t_opt_sub_ham->get_measured_time()      * t_opt_sub_percent);
        tools::log->info("{:<31}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub_hsq->get_name()  ,t_opt_sub_hsq->get_measured_time()    , t_opt_sub_hsq->get_measured_time()      * t_opt_sub_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub_lu->get_name()   ,t_opt_sub_lu->get_measured_time()     , t_opt_sub_lu->get_measured_time()       * t_opt_sub_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub_eig->get_name()  ,t_opt_sub_eig->get_measured_time()    , t_opt_sub_eig->get_measured_time()      * t_opt_sub_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub_bfgs->get_name() ,t_opt_sub_bfgs->get_measured_time()   , t_opt_sub_bfgs->get_measured_time()     * t_opt_sub_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub_vH2->get_name()  ,t_opt_sub_vH2->get_measured_time()    , t_opt_sub_vH2->get_measured_time()      * t_opt_sub_bfgs_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub_vH2v->get_name() ,t_opt_sub_vH2v->get_measured_time()   , t_opt_sub_vH2v->get_measured_time()     * t_opt_sub_bfgs_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub_vH->get_name()   ,t_opt_sub_vH->get_measured_time()     , t_opt_sub_vH->get_measured_time()       * t_opt_sub_bfgs_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", t_opt_sub_vHv->get_name()  ,t_opt_sub_vHv->get_measured_time()    , t_opt_sub_vHv->get_measured_time()      * t_opt_sub_bfgs_percent);
        /* clang-format on */
    }
}

void tools::common::profile::init_profiling() {
    if(t_tot != nullptr) return;
    /* clang-format off */
    t_tot             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "+ Total Time              ");
    t_pre             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocessing          ");
    t_pos             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocessing         ");
    t_sim             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             ");
    t_con             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    ");
    t_eig             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Eig. decomp.          ");
    t_svd             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          ");
    t_evo             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Time evolution        ");
    t_env             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Environment upd.      ");
    t_ent             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  ");
    t_ene             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy                ");
    t_var             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance              ");
    t_prj             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Projections           ");
    t_chk             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Checks                ");
    t_hdf             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          ");
    t_ene_ham         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (HAM)          ");
    t_ene_mom         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (MOM)          ");
    t_var_ham         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (HAM)        ");
    t_var_mom         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (MOM)        ");
    t_mps             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPS         ");
    t_mpo             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPO         ");
    t_opt             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |+ Optimization (xdmrg)  ");
    t_opt_dir         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " ||+ Direct               ");
    t_opt_dir_bfgs    = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |||+ L-BFGS              ");
    t_opt_dir_vH2     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vH2                ");
    t_opt_dir_vH2v    = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vH2v               ");
    t_opt_dir_vH      = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vH                 ");
    t_opt_dir_vHv     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vHv                ");
    t_opt_sub         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " ||+ Subspace             ");
    t_opt_sub_ham     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- Hamiltonian Matrix  ");
    t_opt_sub_hsq     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- Hamiltonian Matrix² ");
    t_opt_sub_lu      = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- LU decomposition    ");
    t_opt_sub_eig     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- Eigenvalue decomp   ");
    t_opt_sub_bfgs    = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |+ L-BFGS              ");
    t_opt_sub_vH2     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vH2                ");
    t_opt_sub_vH2v    = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vH2v               ");
    t_opt_sub_vH      = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vH                 ");
    t_opt_sub_vHv     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vHv                ");
    /* clang-format on */
}

void tools::common::profile::reset_profiling() {
    if(t_tot == nullptr) throw std::runtime_error("Profiling timers have not been initialized");
    // Do not reset total time!
    /* clang-format off */
    t_pre             ->reset(); // |- Preprocessing
    t_pos             ->reset(); // |- Postprocessing
    t_sim             ->reset(); // |+ Simulation
    t_con             ->reset(); //  |- Convergence checks
    t_eig             ->reset(); //  |- Eig. decomp.
    t_svd             ->reset(); //  |- Svd. decomp.
    t_evo             ->reset(); //  |- Time evolution
    t_env             ->reset(); //  |- Environment upd.
    t_ent             ->reset(); //  |- Entanglement entropy
    t_ene             ->reset(); //  |- Energy
    t_var             ->reset(); //  |- Variance
    t_prj             ->reset(); //  |- Projections
    t_chk             ->reset(); //  |- Checks
    t_hdf             ->reset(); //  |- h5pp storage
    t_ene_ham         ->reset(); //  |- Energy (HAM)
    t_ene_mom         ->reset(); //  |- Energy (MOM)
    t_var_ham         ->reset(); //  |- Variance (HAM)
    t_var_mom         ->reset(); //  |- Variance (MOM)
    t_mps             ->reset(); //  |- Multisite-MPS
    t_mpo             ->reset(); //  |- Multisite-MPO
    t_opt             ->reset(); //  |+ Optimization (xdmrg)
    t_opt_dir         ->reset(); //  ||+ Direct
    t_opt_dir_bfgs    ->reset(); //  |||+ L-BFGS steps
    t_opt_dir_vH2     ->reset(); //  || |- vH2
    t_opt_dir_vH2v    ->reset(); //  || |- vH2v
    t_opt_dir_vH      ->reset(); //  || |- vH
    t_opt_dir_vHv     ->reset(); //  || |- vHv
    t_opt_sub         ->reset(); //  ||+ Subspace
    t_opt_sub_ham     ->reset(); //  | |- Hamiltonian Matrix
    t_opt_sub_hsq     ->reset(); //  | |- Hamiltonian Matrix²
    t_opt_sub_lu      ->reset(); //  | |- LU decomposition
    t_opt_sub_eig     ->reset(); //  | |- Eigenvalue decomp
    t_opt_sub_bfgs    ->reset(); //  | |+ L-BFGS steps
    t_opt_sub_vH2     ->reset(); //  |  |- vH2
    t_opt_sub_vH2v    ->reset(); //  |  |- vH2v
    t_opt_sub_vH      ->reset(); //  |  |- vH
    t_opt_sub_vHv     ->reset(); //  |  |- vHv
    /* clang-format on */
}

double tools::common::profile::mem_usage_in_mb(std::string_view name) {
    std::ifstream filestream("/proc/self/status");
    std::string   line;
    while(std::getline(filestream, line)) {
        std::istringstream is_line(line);
        std::string        key;
        if(std::getline(is_line, key, ':')) {
            if(key == name) {
                std::string value_str;
                if(std::getline(is_line, value_str)) {
                    // Extract the number
                    std::string::size_type sz; // alias of size_t
                    int                    value = std::stoi(value_str, &sz);
                    // Now we have the value in kb
                    return value / 1024.0;
                }
            }
        }
    }
    return -1.0;
}

double tools::common::profile::mem_rss_in_mb() { return mem_usage_in_mb("VmRSS"); }
double tools::common::profile::mem_hwm_in_mb() { return mem_usage_in_mb("VmHWM"); }
double tools::common::profile::mem_vm_in_mb() { return mem_usage_in_mb("VmPeak"); }

void tools::common::profile::print_mem_usage() {
    tools::log->info("{:<30}{:>10.1f} MB", "Memory RSS", mem_rss_in_mb());
    tools::log->info("{:<30}{:>10.1f} MB", "Memory Peak", mem_hwm_in_mb());
    tools::log->info("{:<30}{:>10.1f} MB", "Memory Vm", mem_vm_in_mb());
}