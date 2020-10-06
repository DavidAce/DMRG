//
// Created by david on 2019-06-08.
//

#include <fstream>
#include <sstream>

#include <config/enums.h>
#include <config/nmspc_settings.h>

#include <tools/common/log.h>
#include <tools/common/prof.h>

tools::common::profile::internal::MapTicTocUnique &tools::common::profile::get_default_prof() {
    if(not internal::default_algo_type) throw std::runtime_error("No default algorithm type has been set for profiling");
    return prof[internal::default_algo_type.value()];
}
void          tools::common::profile::set_default_prof(AlgorithmType algo_type) { internal::default_algo_type = algo_type; }
AlgorithmType tools::common::profile::get_current_algo_type() {
    if(not internal::default_algo_type) throw std::runtime_error("No default algorithm type has been set");
    return internal::default_algo_type.value();
}

void tools::common::profile::print_profiling_all() { tools::common::profile::print_profiling(std::nullopt); }

void tools::common::profile::print_profiling(std::optional<AlgorithmType> algo_type) {
    if(t_tot == nullptr) return;
    if(settings::profiling::on) {
        auto t_tot_percent = 100.0 / std::max(1.0, t_tot->get_measured_time());
        //        auto t_sim_percent = 100.0 / std::max(1.0, t_sim->get_measured_time());

        /* clang-format off */
        tools::log->info("{:<6}{:<30}{:>10.3f} s"                       ,"DMRG ", t_tot->get_name(), t_tot->get_measured_time());

        if(algo_type)
            for(const auto &[t_key, t_val] : prof[algo_type.value()]) {
                if(t_val->get_measured_time() == 0) continue;
                tools::log->info("{:<8}  {:<30}{:>10.3f} s ({:<3.2f} % of total)", enum2str(algo_type.value()), t_val->get_name(), t_val->get_measured_time(),
                                 t_val->get_measured_time() * t_tot_percent);
            }
        else
            for(const auto &[a_key, a_val] : prof)
                for(const auto &[t_key, t_val] : a_val) {
                    if(t_val->get_measured_time() == 0) continue;
                    tools::log->info("{:<8}  {:<30}{:>10.3f} s ({:<3.2f} % of total)", enum2str(a_key), t_val->get_name(), t_val->get_measured_time(),
                                     t_val->get_measured_time() * t_tot_percent);
                }
    }
}

void tools::common::profile::print_profiling_laps(std::optional<AlgorithmType> algo_type) {
    if(t_tot == nullptr) return;
    if(settings::profiling::on) {
        auto t_tot_percent = 100.0 / std::max(1.0, t_tot->get_lap());
        //        auto t_sim_percent = 100.0 / std::max(1.0, t_sim->get_measured_time());

        /* clang-format off */
        tools::log->info("{:<6}{:<30}{:>10.3f} s"                       ,"DMRG ", t_tot->get_name(), t_tot->get_lap());
        t_tot->start_lap();
        if(algo_type)
            for(const auto &[t_key, t_val] : prof[algo_type.value()]) {
                if(t_val->get_lap() == 0) continue;
                tools::log->info("{:<6}{:<30}{:>10.3f} s ({:<3.2f} % of total)", enum2str(algo_type.value()), t_val->get_name(), t_val->get_lap(),
                                 t_val->get_lap() * t_tot_percent);
                t_val->start_lap();
            }
        else
            for(const auto &[a_key, a_val] : prof)
                for(const auto &[t_key, t_val] : a_val) {
                    if(t_val->get_lap() == 0) continue;
                    tools::log->info("{:<6}{:<30}{:>10.3f} s ({:<3.2f} % of total)", enum2str(a_key), t_val->get_name(), t_val->get_lap(),
                                     t_val->get_lap() * t_tot_percent);
                    t_val->start_lap();
                }
    }
}


void tools::common::profile::init_profiling() {
    if(t_tot != nullptr) return;
    /* clang-format off */
    t_tot             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "+ Total Time              ");

    prof[AlgorithmType::xDMRG] = internal::MapTicTocUnique();
    prof[AlgorithmType::fDMRG] = internal::MapTicTocUnique();
    prof[AlgorithmType::iDMRG] = internal::MapTicTocUnique();
    prof[AlgorithmType::iTEBD] = internal::MapTicTocUnique();

    prof[AlgorithmType::xDMRG]["t_pre"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocess             ");
    prof[AlgorithmType::xDMRG]["t_rnd"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Randomize state        ");
    prof[AlgorithmType::xDMRG]["t_pos"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocess            ");
    prof[AlgorithmType::xDMRG]["t_sim"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             ");
    prof[AlgorithmType::xDMRG]["t_con"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    ");
    prof[AlgorithmType::xDMRG]["t_eig"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Eig. decomp.          ");
    prof[AlgorithmType::xDMRG]["t_svd"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          ");
    prof[AlgorithmType::xDMRG]["t_env"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Environment upd.      ");
    prof[AlgorithmType::xDMRG]["t_ent"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  ");
    prof[AlgorithmType::xDMRG]["t_ene"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy                ");
    prof[AlgorithmType::xDMRG]["t_var"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance              ");
    prof[AlgorithmType::xDMRG]["t_prj"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Projections           ");
    prof[AlgorithmType::xDMRG]["t_chk"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Checks                ");
    prof[AlgorithmType::xDMRG]["t_hdf"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          ");
    prof[AlgorithmType::xDMRG]["t_mps"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPS         ");
    prof[AlgorithmType::xDMRG]["t_mpo"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPO         ");
    prof[AlgorithmType::xDMRG]["t_opt"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |+ Optimization          ");
    prof[AlgorithmType::xDMRG]["t_opt_dir"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " ||+ Direct               ");
    prof[AlgorithmType::xDMRG]["t_opt_dir_bfgs"]    = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |||+ L-BFGS              ");
    prof[AlgorithmType::xDMRG]["t_opt_dir_vH2"]     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vH2                ");
    prof[AlgorithmType::xDMRG]["t_opt_dir_vH2v"]    = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vH2v               ");
    prof[AlgorithmType::xDMRG]["t_opt_dir_vH"]      = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vH                 ");
    prof[AlgorithmType::xDMRG]["t_opt_dir_vHv"]     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vHv                ");
    prof[AlgorithmType::xDMRG]["t_opt_sub"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " ||+ Subspace             ");
    prof[AlgorithmType::xDMRG]["t_opt_sub_ham"]     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- Hamiltonian Matrix  ");
    prof[AlgorithmType::xDMRG]["t_opt_sub_hsq"]     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- Hamiltonian Matrix² ");
    prof[AlgorithmType::xDMRG]["t_opt_sub_lu"]      = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- LU decomposition    ");
    prof[AlgorithmType::xDMRG]["t_opt_sub_eig"]     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- Eigenvalue decomp   ");
    prof[AlgorithmType::xDMRG]["t_opt_sub_bfgs"]    = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |+ L-BFGS              ");
    prof[AlgorithmType::xDMRG]["t_opt_sub_vH2"]     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vH2                ");
    prof[AlgorithmType::xDMRG]["t_opt_sub_vH2v"]    = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vH2v               ");
    prof[AlgorithmType::xDMRG]["t_opt_sub_vH"]      = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vH                 ");
    prof[AlgorithmType::xDMRG]["t_opt_sub_vHv"]     = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vHv                ");


    prof[AlgorithmType::fDMRG]["t_pre"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocess             ");
    prof[AlgorithmType::fDMRG]["t_rnd"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Randomize state        ");
    prof[AlgorithmType::fDMRG]["t_pos"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocess            ");
    prof[AlgorithmType::fDMRG]["t_sim"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             ");
    prof[AlgorithmType::fDMRG]["t_con"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    ");
    prof[AlgorithmType::fDMRG]["t_eig"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Eig. decomp.          ");
    prof[AlgorithmType::fDMRG]["t_svd"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          ");
    prof[AlgorithmType::fDMRG]["t_env"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Environment upd.      ");
    prof[AlgorithmType::fDMRG]["t_ent"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  ");
    prof[AlgorithmType::fDMRG]["t_ene"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy                ");
    prof[AlgorithmType::fDMRG]["t_var"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance              ");
    prof[AlgorithmType::fDMRG]["t_chk"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Checks                ");
    prof[AlgorithmType::fDMRG]["t_hdf"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          ");
    prof[AlgorithmType::fDMRG]["t_mps"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPS         ");
    prof[AlgorithmType::fDMRG]["t_mpo"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPO         ");

    prof[AlgorithmType::iDMRG]["t_pre"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocess             ");
    prof[AlgorithmType::iDMRG]["t_rnd"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Randomize state        ");
    prof[AlgorithmType::iDMRG]["t_pos"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocess            ");
    prof[AlgorithmType::iDMRG]["t_sim"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             ");
    prof[AlgorithmType::iDMRG]["t_con"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    ");
    prof[AlgorithmType::iDMRG]["t_eig"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Eig. decomp.          ");
    prof[AlgorithmType::iDMRG]["t_svd"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          ");
    prof[AlgorithmType::iDMRG]["t_env"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Environment upd.      ");
    prof[AlgorithmType::iDMRG]["t_ent"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  ");
    prof[AlgorithmType::iDMRG]["t_ene"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy                ");
    prof[AlgorithmType::iDMRG]["t_var"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance              ");
    prof[AlgorithmType::iDMRG]["t_chk"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Checks                ");
    prof[AlgorithmType::iDMRG]["t_hdf"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          ");
    prof[AlgorithmType::iDMRG]["t_ene_ham"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (HAM)          ");
    prof[AlgorithmType::iDMRG]["t_ene_mom"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (MOM)          ");
    prof[AlgorithmType::iDMRG]["t_var_ham"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (HAM)        ");
    prof[AlgorithmType::iDMRG]["t_var_mom"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (MOM)        ");

    prof[AlgorithmType::iTEBD]["t_pre"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocessing          ");
    prof[AlgorithmType::iTEBD]["t_pos"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocessing         ");
    prof[AlgorithmType::iTEBD]["t_sim"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             ");
    prof[AlgorithmType::iTEBD]["t_con"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    ");
    prof[AlgorithmType::iTEBD]["t_svd"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          ");
    prof[AlgorithmType::iTEBD]["t_evo"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Time evolution        ");
    prof[AlgorithmType::iTEBD]["t_ent"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  ");
    prof[AlgorithmType::iTEBD]["t_chk"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Checks                ");
    prof[AlgorithmType::iTEBD]["t_hdf"]             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          ");
    prof[AlgorithmType::iTEBD]["t_ene_ham"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (HAM)          ");
    prof[AlgorithmType::iTEBD]["t_ene_mom"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (MOM)          ");
    prof[AlgorithmType::iTEBD]["t_var_ham"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (HAM)        ");
    prof[AlgorithmType::iTEBD]["t_var_mom"]         = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (MOM)        ");
    /* clang-format on */
}

void tools::common::profile::reset_profiling(std::optional<AlgorithmType> algo_type, const std::vector<std::string> &excl) {
    if(t_tot == nullptr) throw std::runtime_error("Profiling timers have not been initialized");
    // Do not reset total time!
    tools::log->trace("Resetting profiling timers");
    if(algo_type)
        for(auto &[t_key, t_val] : prof[algo_type.value()]) {
            for(const auto &e : excl)
                if(t_key == e) continue;
            t_val->reset();
        }
    else
        for(auto &[a_key, a_val] : prof)
            for(auto &[t_key, t_val] : a_val) {
                for(const auto &e : excl)
                    if(t_key == e) continue;
                t_val->reset();
            }
}

//void tools::common::profile::reset_for_run_algorithm(std::optional<AlgorithmType> algo_type, const std::vector<std::string> &excl) {
//    if(t_tot == nullptr) throw std::runtime_error("Profiling timers have not been initialized");
//    if(not algo_type) algo_type = tools::common::profile::get_current_algo_type();
//    // Do not reset total time!
//    tools::log->trace("Resetting profiling timers for new {} run", enum2str(algo_type.value()));
//    if(algo_type)
//        for(auto &[t_key, t_val] : prof[algo_type.value()]) {
//            for(const auto &e : excl)
//                if(t_key == e) continue;
//            t_val->reset();
//        }
//    else
//        for(auto &[a_key, a_val] : prof)
//            for(auto &[t_key, t_val] : a_val) {
//                for(const auto &e : excl)
//                    if(t_key == e) continue;
//                t_val->reset();
//            }
//}

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
                    // Filter non-digit characters
                    value_str.erase(std::remove_if(value_str.begin(), value_str.end(), [](auto const &c) -> bool { return not std::isdigit(c); }),
                                    value_str.end());
                    // Extract the number
                    long long value = 0;
                    try {
                        std::string::size_type sz; // alias of size_t
                        value = std::stoll(value_str, &sz);
                    } catch(const std::exception &ex) {
                        tools::log->error("Could not read mem usage from /proc/self/status: Failed to parse string [{}]: {}", value_str, ex.what());
                    }
                    // Now we have the value in kb
                    return static_cast<double>(value) / 1024.0;
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
