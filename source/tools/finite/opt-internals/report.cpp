
#include <tools/finite/opt.h>
#include <tools/common/log.h>
void tools::finite::opt::internal::reports::print_bfgs_report(){
    if (tools::log->level() > spdlog::level::debug) return;
    if (bfgs_log.empty()) return;
    std::string format_hdr = "{:<24} {:<7} {:<7} {:<20} {:<12} {:<18} {:<18} {:<5} {:<7} {:<18} {:<18}";
    std::string format_num = "- {:<22} {:<7} {:<7} {:<20.16f} {:<12.8f} {:<18.16f} {:<18.16f} {:<5} {:<7} {:<18.4f} {:<18.4f}";
    tools::log->debug(format_hdr.c_str(),
                      "Algorithm",
                      "size",
                      "rank",
                      "energy",
                      "variance",
                      "overlap",
                      "norm",
                      "iter",
                      "counter",
                      "Elapsed time [ms]",
                      "Time per count [ms]");

    for(auto &item : bfgs_log){
        tools::log->debug(format_num.c_str(),
        std::get<0>(item),
        std::get<1>(item),
        std::get<2>(item),
        std::get<3>(item),
        std::real(std::get<4>(item)),
        std::get<5>(item),
        std::get<6>(item),
        std::get<7>(item),
        std::get<8>(item),
        std::get<9>(item)*1000,
        std::get<9>(item)*1000 / (double)std::get<8>(item));
    }
    bfgs_log.clear();
}

void tools::finite::opt::internal::reports::print_eigs_report(){
    if (tools::log->level() > spdlog::level::debug) return;
    if (eigs_log.empty()) return;
    std::string format_hdr = "- {:<5} {:<22} {:<22} {:<23} {:<12} {:<12} {:<12}"; //Thetas are not counted
    std::string format_num = "- {:<5} {:<20.16f} {:<20.16f} {:<21.8f} {:<12.3f} {:<12.3f} {:<12.3f}";
    tools::log->debug(format_hdr.c_str(),
                       "nev",
                       "max <θ_i|θ>",
                       "min <θ_i|θ>",
                       "log₁₀(1-Σ|<θ_i|θ>|²)",
                       "Eig Time[ms]",
                       "Ham Time[ms]",
                       "LU Time[ms]");

    for(auto &item : eigs_log){
        tools::log->debug(format_num.c_str(),
                          std::get<0>(item),
                          std::get<1>(item),
                          std::get<2>(item),
                          std::get<3>(item) ,
                          std::get<4>(item) * 1000,
                          std::get<5>(item) * 1000,
                          std::get<6>(item) * 1000);
    }
    eigs_log.clear();
}




void tools::finite::opt::internal::reports::print_time_report(){
    if (tools::log->level() > spdlog::level::debug) return;
    if(time_log.empty()) return;
    std::string format_hdr = "LBFGS Time report [ms] {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}";
    std::string format_num = "                       {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f}";
    tools::log->debug(format_hdr.c_str(),
                      "vH2v",
                      "vHv",
                      "vH2",
                      "vH",
                      "tot",
                      "op");
    for(auto &item : time_log){
    tools::log->debug(format_num.c_str(),
                     1000 * std::get<0>(item),
                     1000 * std::get<1>(item),
                     1000 * std::get<2>(item),
                     1000 * std::get<3>(item),
                     1000 *(std::get<0>(item)
                          + std::get<1>(item)
                          + std::get<2>(item)
                          + std::get<3>(item)),
                     1000 * std::get<4>(item));
    }
    time_log.clear();
}