
#include <tools/finite/opt.h>
#include <tools/common/log.h>
void tools::finite::opt::internal::reports::print_report(const std::vector<direct_opt_tuple> &opt_log){
    if (tools::log->level() > spdlog::level::debug) return;
    std::string format_hdr = "{:<18} {:<7} {:<20} {:<12} {:<18} {:<18} {:<5} {:<7} {:<18} {:<18}";
    std::string format_num = "- {:<16} {:<7} {:<20.16f} {:<12.8f} {:<18.16f} {:<18.16f} {:<5} {:<7} {:<18.4f} {:<18.4f}";
    tools::log->debug(format_hdr,
                      "Algorithm",
                      "size",
                      "energy",
                      "variance",
                      "overlap",
                      "norm",
                      "iter",
                      "counter",
                      "Elapsed time [ms]",
                      "Time per count [ms]");

    for(auto &item : opt_log){
        tools::log->debug(format_num,
        std::get<0>(item),
        std::get<1>(item),
        std::get<2>(item),
        std::real(std::get<3>(item)),
        std::get<4>(item),
        std::get<5>(item),
        std::get<6>(item),
        std::get<7>(item),
        std::get<8>(item)*1000,
        std::get<8>(item)*1000 / (double)std::get<7>(item));
    }
}


void tools::finite::opt::internal::reports::print_report(const std::vector<subspc_opt_tuple> &opt_log){
    if (tools::log->level() > spdlog::level::debug) return;
    std::string format_hdr = "{:<18} {:<7} {:<20} {:<12} {:<18} {:<18} {:<5} {:<7} {:<18} {:<18}";
    std::string format_num = "- {:<16} {:<7} {:<20.16f} {:<12.8f} {:<18.16f} {:<18.16f} {:<5} {:<7} {:<18.4f} {:<18.4f}";
    tools::log->debug(format_hdr,
                     "Algorithm",
                     "size",
                     "energy",
                     "variance",
                     "overlap",
                     "norm",
                     "iter",
                     "counter",
                     "Elapsed time [ms]",
                     "Time per count [ms]");
    for(auto &item : opt_log){
    tools::log->debug(format_num,
                      std::get<0>(item),
                      std::get<1>(item),
                      std::get<2>(item),
                      std::real(std::get<3>(item)),
                      std::get<4>(item),
                      std::get<5>(item),
                      std::get<6>(item),
                      std::get<7>(item),
                      std::get<8>(item) * 1000,
                      std::get<8>(item) * 1000 / (double)std::get<7>(item));
    }
}



void tools::finite::opt::internal::reports::print_report(const std::vector<eig_tuple> &eig_log){
    if (tools::log->level() > spdlog::level::debug) return;
    std::string format_hdr = "- {:<5} {:<22} {:<22} {:<23} {:<12} {:<12}"; //Thetas are not counted
    std::string format_num = "- {:<5} {:<20.16f} {:<20.16f} {:<21.8f} {:<12.3f} {:<12.3f}";
    tools::log->debug(format_hdr,
                       "nev",
                       "max <θ_i|θ>",
                       "min <θ_i|θ>",
                       "log10(1-Σ|<θ_i|θ>|^2)",
                       "Eig Time[ms]",
                       "LU Time[ms]");

    for(auto &item : eig_log){
        tools::log->debug(format_num,
                          std::get<0>(item),
                          std::get<1>(item),
                          std::get<2>(item),
                          std::get<3>(item) ,
                          std::get<4>(item) * 1000,
                          std::get<5>(item) * 1000);
    }
}




void tools::finite::opt::internal::reports::print_report(const lbfgs_tuple lbfgs_log){
    if (tools::log->level() > spdlog::level::debug) return;
    std::string format_hdr = "LBFGS Time report [ms] {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}";
    std::string format_num = "                       {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f} {:<10.3f}";
    tools::log->debug(format_hdr,
                      "vH2v",
                      "vHv",
                      "vH2",
                      "vH",
                      "tot",
                      "op");

    tools::log->debug(format_num,
                     1000 * std::get<0>(lbfgs_log),
                     1000 * std::get<1>(lbfgs_log),
                     1000 * std::get<2>(lbfgs_log),
                     1000 * std::get<3>(lbfgs_log),
                     1000 *(std::get<0>(lbfgs_log)
                          + std::get<1>(lbfgs_log)
                          + std::get<2>(lbfgs_log)
                          + std::get<3>(lbfgs_log)),
                     1000 * std::get<4>(lbfgs_log));

}