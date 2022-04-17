
#include "report.h"
#include "../opt_mps.h"
#include <general/iter.h>
#include <tid/tid.h>
#include <tools/common/log.h>

/* clang-format off */
void tools::finite::opt::reports::print_bfgs_report(){
    if (tools::log->level() > spdlog::level::debug) return;
    if (bfgs_log.empty()) return;
    tools::log->debug(FMT_STRING("{:<52} {:<7} {:<7} {:<22} {:<8} {:<18} {:<18} {:<8} {:<5} {:<7} {:<8} {:<8} {:<10} {:<10}"),
                      "Optimization report",
                      "size",
                      "space",
                      "E/L",
                      "σ²H", // Special characters are counted properly in fmt 1.7.0
                      "overlap",
                      "norm",
                      "rnorm",
                      "iter",
                      "ops",
                      "|Δf|",
                      "|∇Ψ|ᵐᵃˣ",
                      "time [s]",
                      "avg [op/s]");
    for(auto &entry : bfgs_log){
        tools::log->debug(FMT_STRING("- {:<50} {:<7} {:<7} {:<22.16f} {:<8.2e} {:<18.15f} {:<18.15f} {:<8.2e} {:<5} {:<7} {:<8.2e} {:<8.2e} {:<10.2e} {:<10.2e}"),
        entry.description, entry.size,
        entry.space, entry.energy,
        entry.variance,
        entry.overlap,entry.norm, entry.rnorm,
        entry.iter, entry.counter,
        entry.delta_f, entry.max_grad_norm,
        entry.time,
        static_cast<double>(entry.counter)/entry.time);
    }
    bfgs_log.clear();
}

void tools::finite::opt::reports::print_subs_report(){
    if (tools::log->level() > spdlog::level::debug) return;
    if (subs_log.empty()) return;
    tools::log->debug(FMT_STRING("- {:<5} {:<18} {:<18} {:<18} {:<11} {:<11} {:<11} {:<6} {:<6} {:<6}"),
                       "nev",
                       "max <φ_i|ψ>",
                       "min <φ_i|ψ>",
                       "ε:(1-Σ|<φ_i|ψ>|²)",  // Special characters are counted properly in spdlog 1.7.0
                       "eig time[s]",
                       "ham time[s]",
                       "lu Time[s]",
                       "iter",
                       "mv",
                       "pc");

    for(auto &entry : subs_log){
        tools::log->debug(FMT_STRING("- {:<5} {:<18.16f} {:<18.16f} {:<18.2e} {:<11.2e} {:<11.2e} {:<11.2e} {:<6} {:<6} {:<6}"),
                          entry.nev,
                          entry.max_olap,
                          entry.min_olap,
                          entry.eps ,
                          entry.eig_time,
                          entry.ham_time,
                          entry.lu_time ,
                          entry.iter,
                          entry.mv,
                          entry.pc
                          );
    }
    subs_log.clear();
}




void tools::finite::opt::reports::print_time_report(){
    if (tools::log->level() > spdlog::level::trace) return;
    if(time_log.empty()) return;
    tools::log->trace(FMT_STRING("BFGS Time report [s] {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}"),
                      "<ψ|H²|ψ>",
                      "<ψ|H|ψ>",
                      "H²|ψ>",
                      "H|ψ>",
                      "sum",
                      "step",
                      "l-bfgs");
    for(auto &entry : time_log){
    tools::log->trace(FMT_STRING("{:<21} {:<10.2e} {:<10.2e} {:<10.2e} {:<10.2e} {:<10.2e} {:<10.2e} {:<10.2e}"),
                     " ",
                      entry.vH2v,
                      entry.vHv,
                      entry.vH2,
                      entry.vH,
                      (
                        entry.vH2v +
                        entry.vHv +
                        entry.vH2 +
                        entry.vH
                      ),
                     entry.step,
                     entry.bfgs);
    }
    time_log.clear();
}

void tools::finite::opt::reports::print_eigs_report(std::optional<size_t> max_entries){
    if (eigs_log.empty()) return;
    auto level = eigs_log.front().level;
    if (level < tools::log->level()) {
        eigs_log.clear();
        return;
    }
    tools::log->log(level, FMT_STRING("{:<52} {:<7} {:<4} {:<4} {:<4} {:<4} {:<8} {:<8} {:<22} {:<22} {:<8} {:<18} {:<18} {:<8} {:<5} {:<7} {:<7} {:<10} {:<10}"),
                      "Optimization report",
                      "size",
                      "ritz",
                      "idx",
                      "nev",
                      "ncv",
                      "tol",
                      "|∇Ψ|ᵐᵃˣ",
                      "E/L",
                      "λ",
                      "σ²H", // Special characters are counted properly in fmt 1.7.0
                      "overlap",
                      "norm",
                      "rnorm",
                      "iter",
                      "mv",
                      "pc",
                      "time [s]",
                      "avg [mv/s]");

    for(const auto &[idx,entry] : iter::enumerate(eigs_log)){
        if(max_entries and max_entries.value() <= idx) break;
        tools::log->log(level, FMT_STRING("- {:<50} {:<7} {:<4} {:<4} {:<4} {:<4} {:<8.2e} {:<8.2e} {:<+22.15f} {:<+22.15f} {:<8.2e} {:<18.15f} {:<18.15f} {:<8.2e} {:<5} {:<7} {:<7} {:<10.2e} {:<10.2e}"),
                          entry.description,
                          entry.size, entry.ritz,entry.idx, entry.nev, entry.ncv, entry.tol, entry.grad,
                          entry.energy,entry.eigval,
                          entry.variance,
                          entry.overlap,entry.norm, entry.rnorm,
                          entry.iter, entry.mv, entry.pc,
                          entry.time,
                          static_cast<double>(entry.mv)/entry.time);
    }
    eigs_log.clear();
}

/* clang-format on */
void tools::finite::opt::reports::bfgs_add_entry(const std::string &description, long size, long space, double energy, double variance, double overlap,
                                                 double norm, double rnorm, double delta_f, double grad_norm, size_t iter, size_t counter, double time) {
    if(tools::log->level() > spdlog::level::debug) return;
    bfgs_log.push_back({description, size, space, energy, variance, overlap, norm, rnorm, delta_f, grad_norm, iter, counter, time});
}
void tools::finite::opt::reports::bfgs_add_entry(std::string_view mode, std::string_view tag, const opt_mps &mps, std::optional<long> space) {
    if(tools::log->level() > spdlog::level::debug) return;
    if(not space) space = mps.get_tensor().size();
    std::string description = fmt::format("{:<8} {:<16} {}", mode, mps.get_name(), tag);
    bfgs_log.push_back(bfgs_entry{description, mps.get_tensor().size(), space.value(), mps.get_energy_per_site(), mps.get_variance(), mps.get_overlap(),
                                  mps.get_norm(), mps.get_eigs_rnorm(), mps.get_delta_f(), mps.get_grad_max(), mps.get_iter(), mps.get_mv(), mps.get_time()});
}

void tools::finite::opt::reports::time_add_entry() {
    if(tools::log->level() > spdlog::level::trace) return;
    time_log.push_back({tid::get("vH2v").get_last_interval(), tid::get("vHv").get_last_interval(), tid::get("vH2").get_last_interval(),
                        tid::get("vH").get_last_interval(), tid::get("bfgs").get_last_interval(), tid::get("step").get_last_interval()});
}

void tools::finite::opt::reports::subs_add_entry(long nev, double max_olap, double min_olap, double eps, double eig_time, double ham_time, double lu_time,
                                                 long iter, long mv, long pc) {
    if(tools::log->level() > spdlog::level::debug) return;
    subs_log.push_back({nev, max_olap, min_olap, eps, eig_time, ham_time, lu_time, iter, mv, pc});
}

void tools::finite::opt::reports::eigs_add_entry(const opt_mps &mps, spdlog::level::level_enum level) {
    if(level < tools::log->level()) return;
    std::string description = fmt::format("{:<24}", mps.get_name());
    eigs_log.push_back(eigs_entry{description, std::string(mps.get_eigs_ritz()), mps.get_tensor().size(), mps.get_eigs_idx(), mps.get_eigs_nev(),
                                  mps.get_eigs_ncv(), mps.get_energy_per_site(), mps.get_eigs_eigval(), mps.get_variance(), mps.get_overlap(), mps.get_norm(),
                                  mps.get_eigs_rnorm(), mps.get_eigs_tol(), mps.get_grad_max(), mps.get_iter(), mps.get_mv(), mps.get_pc(), mps.get_time(),
                                  level});
}
