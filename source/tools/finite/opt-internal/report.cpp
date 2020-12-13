
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <tools/finite/opt-internal/report.h>
#include <tools/finite/opt_state.h>

/* clang-format off */
void tools::finite::opt::internal::reports::print_bfgs_report(){
    if (tools::log->level() > spdlog::level::debug) return;
    if (bfgs_log.empty()) return;
    std::string format_hdr = "{:<52} {:<7} {:<7} {:<20} {:<12} {:<18} {:<18} {:<5} {:<7} {:<18} {:<18}";
    std::string format_num = "- {:<50} {:<7} {:<7} {:<20.15f} {:<12.8f} {:<18.15f} {:<18.15f} {:<5} {:<7} {:<18.3f} {:<18.3f}";
    tools::log->debug(format_hdr.c_str(),
                      "Optimization report",
                      "size",
                      "space",
                      "energy/L",
                      "log₁₀ var", // Special characters are counted properly in fmt 1.7.0
                      "overlap",
                      "norm",
                      "iter",
                      "counter",
                      "Elapsed time [ms]",
                      "Time per count [ms]");

    for(auto &entry : bfgs_log){
        tools::log->debug(format_num.c_str(),
        entry.description,
        entry.size,entry.space,entry.energy,std::log10(entry.variance),entry.overlap,entry.norm,entry.iter,entry.counter,entry.time*1000,entry.time*1000.0/std::max(1.0,static_cast<double>(entry.counter)));
    }
    bfgs_log.clear();
}

void tools::finite::opt::internal::reports::print_eigs_report(){
    if (tools::log->level() > spdlog::level::debug) return;
    if (eigs_log.empty()) return;
    std::string format_hdr = "- {:<5} {:<20} {:<20} {:<20} {:<12} {:<12} {:<12} {:<6}";
    std::string format_num = "- {:<5} {:<20.15f} {:<20.15f} {:<20.8f} {:<12.3f} {:<12.3f} {:<12.3f} {:<6}";
    tools::log->debug(format_hdr.c_str(),
                       "nev",
                       "max <θ_i|θ>",
                       "min <θ_i|θ>",
                       "log₁₀(1-Σ|<θ_i|θ>|²)",  // Special characters are counted properly in fmt 1.7.0
                       "Eig Time[ms]",
                       "Ham Time[ms]",
                       "LU Time[ms]",
                       "Steps" );

    for(auto &entry : eigs_log){
        tools::log->debug(format_num.c_str(),
                          entry.nev,
                          entry.max_olap,
                          entry.min_olap,
                          entry.eps ,
                          entry.eig_time * 1000,
                          entry.ham_time * 1000,
                          entry.lu_time  * 1000,
                          entry.steps
                            );
    }
    eigs_log.clear();
}




void tools::finite::opt::internal::reports::print_time_report(){
    if (tools::log->level() > spdlog::level::trace) return;
    if(time_log.empty()) return;
    std::string format_hdr = "LBFGS Time report [ms] {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}";
    std::string format_num = "                       {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f}";
    tools::log->trace(format_hdr.c_str(),
                      "vH2v",
                      "vHv",
                      "vH2",
                      "vH",
                      "sum",
                      "step",
                      "l-bfgs");
    for(auto &entry : time_log){
    tools::log->trace(format_num.c_str(),
                     1000 * entry.vH2v,
                     1000 * entry.vHv,
                     1000 * entry.vH2,
                     1000 * entry.vH,
                     1000 *(entry.vH2v
                          + entry.vHv
                          + entry.vH2
                          + entry.vH),
                     1000 * entry.step,
                     1000 * entry.bfgs);
    }
    time_log.clear();
}




void tools::finite::opt::internal::reports::bfgs_add_entry(const std::string & description, long size,long space, double energy, double variance,
                    double overlap, double norm, size_t iter, size_t counter, double time){
    if (tools::log->level() > spdlog::level::debug) return;
    bfgs_log.push_back({description, size,space, energy, variance, overlap, norm, iter, counter, time});
}
void tools::finite::opt::internal::reports::bfgs_add_entry(const std::string & mode,const std::string & tag, const opt_state & tensor, std::optional<long> space){
    if (tools::log->level() > spdlog::level::debug) return;
    if(not space) space = tensor.get_tensor().size();
    std::string description = fmt::format("{:<8} {:<16} {}",mode,tensor.get_name(),tag);
    bfgs_log.push_back({description, tensor.get_tensor().size(),space.value(),tensor.get_energy_per_site(), tensor.get_variance(), tensor.get_overlap(), tensor.get_norm(),
                        tensor.get_iter(), tensor.get_counter(), tensor.get_time()});
}

void tools::finite::opt::internal::reports::time_add_dir_entry(){
    if (tools::log->level() > spdlog::level::trace) return;
    time_log.push_back({  tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vH2v"]->get_last_interval(),
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vHv"]->get_last_interval(),
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vH2"]->get_last_interval() ,
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_vH"]->get_last_interval(),
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_bfgs"]->get_last_interval(),
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_dir_step"]->get_last_interval()});
}
void tools::finite::opt::internal::reports::time_add_sub_entry(){
    if (tools::log->level() > spdlog::level::trace) return;
    time_log.push_back ({ tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2v"]->get_last_interval(),
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vHv"]->get_last_interval(),
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH2"]->get_last_interval(),
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_vH"]->get_last_interval(),
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_bfgs"]->get_last_interval(),
                          tools::common::profile::prof[AlgorithmType::xDMRG]["t_opt_sub_step"]->get_last_interval()});
}
void tools::finite::opt::internal::reports::eigs_add_entry(long nev, double max_olap, double min_olap, double eps, double eig_time,double ham_time, double lu_time, size_t steps){
    if (tools::log->level() > spdlog::level::debug) return;
    eigs_log.push_back({nev, max_olap, min_olap, eps, eig_time,ham_time,lu_time,steps});
}
