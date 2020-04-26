//
// Created by david on 2019-06-08.
//

#include <simulation/nmspc_settings.h>
#include <tools/common/prof.h>
#include <tools/common/log.h>
#include <fstream>
#include <sstream>

void tools::common::profile::print_profiling() {
    if(settings::profiling::on) {
        auto t_tot_percent = 100.0/std::max(1.0,t_tot->get_age());
        auto t_sim_percent = 100.0/std::max(1.0,t_sim->get_measured_time());
        /* clang-format off */
        tools::log->info("{:<30}{:>10.3f} s"                      , "+ Total Time              ",t_tot->get_age());
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", "|- Preprocessing          ",t_pre->get_measured_time(), t_pre->get_measured_time() * t_tot_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", "|- Postprocessing         ",t_pos->get_measured_time(), t_pos->get_measured_time() * t_tot_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", "|+ Simulation             ",t_sim->get_measured_time(), t_sim->get_measured_time() * t_tot_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Convergence checks    ",t_con->get_measured_time(), t_con->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Eig. decomp.          ",t_eig->get_measured_time(), t_eig->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Svd. decomp.          ",t_svd->get_measured_time(), t_svd->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Hamiltonian Matrix    ",t_ham->get_measured_time(), t_ham->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Hamiltonian Matrix Sq ",t_hsq->get_measured_time(), t_hsq->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Multisite-MPO         ",t_mpo->get_measured_time(), t_mpo->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Optimization (Ceres)  ",t_opt->get_measured_time(), t_opt->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Time evolution        ",t_evo->get_measured_time(), t_evo->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Environment upd.      ",t_env->get_measured_time(), t_env->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Entanglement entropy  ",t_ent->get_measured_time(), t_ent->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Energy                ",t_ene->get_measured_time(), t_ene->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Variance              ",t_var->get_measured_time(), t_var->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Projections           ",t_prj->get_measured_time(), t_prj->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Checks                ",t_chk->get_measured_time(), t_chk->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- h5pp storage          ",t_hdf->get_measured_time(), t_hdf->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Energy (HAM)          ",t_ene_ham->get_measured_time(), t_ene_ham->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Energy (MOM)          ",t_ene_mom->get_measured_time(), t_ene_mom->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Variance (HAM)        ",t_var_ham->get_measured_time(), t_var_ham->get_measured_time() * t_sim_percent);
        tools::log->info("{:<30}{:>10.3f} s ({:<3.2f} % of parent)", " |- Variance (MOM)        ",t_var_mom->get_measured_time(), t_var_mom->get_measured_time() * t_sim_percent);
        /* clang-format on */
    }
}

void tools::common::profile::init_profiling() {
    t_tot     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "+ Total Time              ");
    t_pre     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocessing          ");
    t_pos     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocessing         ");
    t_sim     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             ");
    t_con     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    ");
    t_eig     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Eig. decomp.          ");
    t_svd     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          ");
    t_ham     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Hamiltonian Matrix    ");
    t_hsq     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Hamiltonian Matrix Sq ");
    t_mpo     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPO         ");
    t_opt     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Optimization (Ceres)  ");
    t_evo     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Time evolution        ");
    t_env     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Environment upd.      ");
    t_ent     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  ");
    t_ene     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy                ");
    t_var     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance              ");
    t_prj     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Projections           ");
    t_chk     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Checks                ");
    t_hdf     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          ");
    t_ene_ham = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (HAM)          ");
    t_ene_mom = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (MOM)          ");
    t_var_ham = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (HAM)        ");
    t_var_mom = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (MOM)        ");

    t_vH2v = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "t_vH2v");
    t_vHv  = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "t_vHv ");
    t_vH2  = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "t_vH2 ");
    t_vH   = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "t_vH  ");
    t_op   = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "t_op  ");
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
double tools::common::profile::mem_vm_in_mb()  { return mem_usage_in_mb("VmPeak"); }

void tools::common::profile::print_mem_usage() {
    tools::log->info("{:<30}{:>10.1f} MB", "Memory RSS", mem_rss_in_mb());
    tools::log->info("{:<30}{:>10.1f} MB", "Memory Peak", mem_hwm_in_mb());
    tools::log->info("{:<30}{:>10.1f} MB", "Memory Vm", mem_vm_in_mb());
}