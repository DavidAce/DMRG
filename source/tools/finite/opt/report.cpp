
#include "report.h"
#include "../opt_mps.h"
#include <general/iter.h>
#include <tid/tid.h>
#include <tools/common/log.h>

/* clang-format off */
void tools::finite::opt::reports::print_bfgs_report(){
    if (tools::log->level() > spdlog::level::debug) return;
    if (bfgs_log.empty()) return;
    tools::log->debug(FMT_STRING("{:<52} {:<7} {:<7} {:<22} {:<8} {:<18} {:<18} {:<5} {:<7} {:<8} {:<8} {:<10} {:<10}"),
                      "Optimization report",
                      "size",
                      "space",
                      "E/L",
                      "σ²H", // Special characters are counted properly in fmt 1.7.0
                      "overlap",
                      "norm",
                      "iter",
                      "ops",
                      "|Δf|",
                      "∇fᵐᵃˣ",
                      "time [s]",
                      "avg [s/op]");
    for(auto &entry : bfgs_log){
        tools::log->debug(FMT_STRING("- {:<50} {:<7} {:<7} {:<22.16f} {:<8.2e} {:<18.15f} {:<18.15f} {:<5} {:<7} {:<8.2e} {:<8.2e} {:<10.2e} {:<10.2e}"),
        entry.description, entry.size,
        entry.space, entry.energy,
        entry.variance,
        entry.overlap,entry.norm,
        entry.iter, entry.counter,
        entry.delta_f, entry.max_grad_norm,
        entry.time,
        entry.time/std::max(1.0,static_cast<double>(entry.counter)));
    }
    bfgs_log.clear();
}

void tools::finite::opt::reports::print_eigs_report(){
    if (tools::log->level() > spdlog::level::debug) return;
    if (eigs_log.empty()) return;
    tools::log->debug(FMT_STRING("- {:<5} {:<18} {:<18} {:<18} {:<11} {:<11} {:<11} {:<6}"),
                       "nev",
                       "max <φ_i|ψ>",
                       "min <φ_i|ψ>",
                       "ε:(1-Σ|<φ_i|ψ>|²)",  // Special characters are counted properly in spdlog 1.7.0
                       "eig time[s]",
                       "ham time[s]",
                       "lu Time[s]",
                       "steps" );

    for(auto &entry : eigs_log){
        tools::log->debug(FMT_STRING("- {:<5} {:<18.16f} {:<18.16f} {:<18.2e} {:<11.2e} {:<11.2e} {:<11.2e} {:<6}"),
                          entry.nev,
                          entry.max_olap,
                          entry.min_olap,
                          entry.eps ,
                          entry.eig_time,
                          entry.ham_time,
                          entry.lu_time ,
                          entry.steps
                            );
    }
    eigs_log.clear();
}




void tools::finite::opt::reports::print_time_report(){
    if (tools::log->level() > spdlog::level::trace) return;
    if(time_log.empty()) return;
    tools::log->trace(FMT_STRING("LBFGS Time report [s] {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}"),
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

void tools::finite::opt::reports::print_krylov_report(std::optional<size_t> max_entries){
    if (tools::log->level() > spdlog::level::debug) return;
    if (krylov_log.empty()) return;
    tools::log->debug(FMT_STRING("{:<52} {:<7} {:<4} {:<4} {:<4} {:<4} {:<8} {:<8} {:<8} {:<22} {:<22} {:<8} {:<18} {:<18} {:<5} {:<7} {:<10} {:<10}"),
                      "Optimization report",
                      "size",
                      "ritz",
                      "idx",
                      "nev",
                      "ncv",
                      "tol",
                      "res",
                      "∇fᵐᵃˣ",
                      "E/L",
                      "λ",
                      "σ²H", // Special characters are counted properly in fmt 1.7.0
                      "overlap",
                      "norm",
                      "iter",
                      "mv",
                      "time [s]",
                      "avg [s/op]");

    for(const auto &[idx,entry] : iter::enumerate(krylov_log)){
        if(max_entries and max_entries.value() <= idx) break;
        tools::log->debug(FMT_STRING("- {:<50} {:<7} {:<4} {:<4} {:<4} {:<4} {:<8.2e} {:<8.2e} {:<8.2e} {:<22.15f} {:<22.15f} {:<8.2e} {:<18.15f} {:<18.15f} {:<5} {:<7} {:<10.2e} {:<10.2e}"),
                          entry.description,
                          entry.size, entry.ritz,entry.idx, entry.nev, entry.ncv, entry.tol, entry.resid, entry.grad,
                          entry.energy,entry.eigval,
                          entry.variance,
                          entry.overlap,entry.norm,
                          entry.iter, entry.counter,
                          entry.time,
                          entry.time/std::max(1.0,static_cast<double>(entry.counter)));
    }
    krylov_log.clear();
}

/* clang-format on */
void tools::finite::opt::reports::bfgs_add_entry(const std::string &description, long size, long space, double energy, double variance, double overlap,
                                                 double norm, double delta_f, double grad_norm, size_t iter, size_t counter, double time) {
    if(tools::log->level() > spdlog::level::debug) return;
    bfgs_log.push_back({description, size, space, energy, variance, overlap, norm, delta_f, grad_norm, iter, counter, time});
}
void tools::finite::opt::reports::bfgs_add_entry(std::string_view mode, std::string_view tag, const opt_mps &mps, std::optional<long> space) {
    if(tools::log->level() > spdlog::level::debug) return;
    if(not space) space = mps.get_tensor().size();
    std::string description = fmt::format("{:<8} {:<16} {}", mode, mps.get_name(), tag);
    bfgs_log.push_back(bfgs_entry{description, mps.get_tensor().size(), space.value(), mps.get_energy_per_site(), mps.get_variance(), mps.get_overlap(),
                                  mps.get_norm(), mps.get_delta_f(), mps.get_max_grad(), mps.get_iter(), mps.get_mv(), mps.get_time()});
}

void tools::finite::opt::reports::time_add_opt_entry() {
    if(tools::log->level() > spdlog::level::trace) return;
    time_log.push_back({tid::get("vH2v").get_last_interval(), tid::get("vHv").get_last_interval(), tid::get("vH2").get_last_interval(),
                        tid::get("vH").get_last_interval(), tid::get("lbfgs").get_last_interval(), tid::get("step").get_last_interval()});
}

void tools::finite::opt::reports::eigs_add_entry(long nev, double max_olap, double min_olap, double eps, double eig_time, double ham_time, double lu_time,
                                                 size_t steps) {
    if(tools::log->level() > spdlog::level::debug) return;
    eigs_log.push_back({nev, max_olap, min_olap, eps, eig_time, ham_time, lu_time, steps});
}

void tools::finite::opt::reports::krylov_add_entry(const opt_mps &mps) {
    if(tools::log->level() > spdlog::level::debug) return;
    std::string description = fmt::format("krylov {:<24}", mps.get_name());
    krylov_log.push_back(krylov_entry{description, std::string(mps.get_krylov_ritz()), mps.get_tensor().size(), mps.get_krylov_idx(), mps.get_krylov_nev(),
                                      mps.get_krylov_ncv(), mps.get_energy_per_site(), mps.get_krylov_eigval(), mps.get_variance(), mps.get_overlap(),
                                      mps.get_norm(), mps.get_krylov_tol(), mps.get_krylov_resid(), mps.get_max_grad(), mps.get_iter(), mps.get_mv(),
                                      mps.get_time()});
}
