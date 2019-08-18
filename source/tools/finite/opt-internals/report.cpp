//
// Created by david on 2019-05-31.
//

#include <tools/finite/opt.h>

void tools::finite::opt::internals::reports::print_report(const std::vector<direct_opt_tuple> &opt_log){
    std::stringstream report;
    report    << std::setprecision(16) << '\n'
              <<"    "<< std::setw(24) << std::left << "Algorithm"
              <<"    "<< std::setw(8)  << std::left << "size"
              <<"    "<< std::setw(24) << std::left << "energy"
              <<"    "<< std::setw(44) << std::left << "variance"
              <<"    "<< std::setw(24) << std::left << "overlap"
              <<"    "<< std::setw(24) << std::left << "norm"
              <<"    "<< std::setw(8)  << std::left << "iter"
              <<"    "<< std::setw(8)  << std::left << "counter"
              <<"    "<< std::setw(20) << std::left << "Elapsed time [ms]"
              <<"    "<< std::setw(20) << std::left << "Time per count [ms]"
              << '\n';
    for(auto &item : opt_log){
        report << std::setprecision(16) << std::fixed
               << "    " << std::setw(24) << std::left << std::fixed << std::get<0>(item)
               << "    " << std::setw(8)  << std::left << std::fixed << std::get<1>(item)
               << "    " << std::setw(24) << std::left << std::fixed << std::get<2>(item)
               << "    " << std::setw(44) << std::left << std::fixed << std::get<3>(item)
               << "    " << std::setw(24) << std::left << std::fixed << std::get<4>(item)
               << "    " << std::setw(24) << std::left << std::fixed << std::get<5>(item)
               << "    " << std::setw(8)  << std::left << std::fixed << std::get<6>(item) << std::setprecision(3)
               << "    " << std::setw(8)  << std::left << std::fixed << std::get<7>(item) << std::setprecision(3)
               << "    " << std::setw(20) << std::left << std::fixed << std::get<8>(item) * 1000
               << "    " << std::setw(20) << std::left << std::fixed << std::get<8>(item) * 1000 / (double)std::get<7>(item)
               << '\n';
    }
    report << '\n';
    tools::log->debug(report.str());
}


void tools::finite::opt::internals::reports::print_report(const std::vector<subspc_opt_tuple> &opt_log){
    std::stringstream report;
    report    << std::setprecision(16) << '\n'
              <<"    "<< std::setw(24) << std::left << "Algorithm"
              <<"    "<< std::setw(8)  << std::left << "size"
              <<"    "<< std::setw(24) << std::left << "energy"
              <<"    "<< std::setw(44) << std::left << "variance"
              <<"    "<< std::setw(24) << std::left << "overlap"
              <<"    "<< std::setw(24) << std::left << "norm"
              <<"    "<< std::setw(8)  << std::left << "iter"
              <<"    "<< std::setw(8)  << std::left << "counter"
              <<"    "<< std::setw(20) << std::left << "Elapsed time [ms]"
              <<"    "<< std::setw(20) << std::left << "Time per count [ms]"
              << '\n';
    for(auto &item : opt_log){
        report << std::setprecision(16) << std::fixed
               << "    " << std::setw(24) << std::left << std::fixed << std::get<0>(item)
               << "    " << std::setw(8)  << std::left << std::fixed << std::get<1>(item)
               << "    " << std::setw(24) << std::left << std::fixed << std::get<2>(item)
               << "    " << std::setw(44) << std::left << std::fixed << std::get<3>(item)
               << "    " << std::setw(24) << std::left << std::fixed << std::get<4>(item)
               << "    " << std::setw(24) << std::left << std::fixed << std::get<5>(item)
               << "    " << std::setw(8)  << std::left << std::fixed << std::get<6>(item) << std::setprecision(3)
               << "    " << std::setw(8)  << std::left << std::fixed << std::get<7>(item) << std::setprecision(3)
               << "    " << std::setw(20) << std::left << std::fixed << std::get<8>(item) * 1000
               << "    " << std::setw(20) << std::left << std::fixed << std::get<8>(item) * 1000 / (double)std::get<7>(item)
               << '\n';
    }
    tools::log->debug(report.str());

}



void tools::finite::opt::internals::reports::print_report(const std::vector<eig_tuple> &eig_log){
    std::stringstream solver_report;
    solver_report << '\n'
                  << std::setw(12) << std::right << "n eigvecs"
                  << std::setw(24) << std::right << "max overlap"
                  << std::setw(24) << std::right << "min overlap"
                  << std::setw(34) << std::right << "eps = Σ_i |<state_i|old>|^2"
                  << std::setw(32) << std::right << "quality = log10(1 - eps)"
                  << std::setw(18) << std::right << "Eig Time[ms]"
                  << std::setw(18) << std::right << "LU  Time[ms]"
                  << '\n';
    for(auto &item : eig_log){
        solver_report
                << std::setprecision(16) << std::fixed
                << std::setw(12) << std::right << std::get<0>(item)
                << std::setw(24) << std::right << std::get<1>(item)
                << std::setw(24) << std::right << std::get<2>(item)
                << std::setw(34) << std::right << std::get<3>(item)
                << std::setw(32) << std::right << std::get<4>(item) << std::setprecision(3)
                << std::setw(18) << std::right << std::get<5>(item) * 1000
                << std::setw(18) << std::right << std::get<6>(item) * 1000
                << '\n';
    }
    solver_report << '\n' << std::flush;
    tools::log->debug(solver_report.str());

}




void tools::finite::opt::internals::reports::print_report(const lbfgs_tuple lbfgs_log){
    std::stringstream report;
    report
            << std::setprecision(3) << '\n'
            << "    " << std::setw(24) << std::left << "LBFGS Time report"
            << "    " << std::setw(12) << std::left << "vH2v  [ms]"
            << "    " << std::setw(12) << std::left << "vHv  [ms]"
            << "    " << std::setw(12) << std::left << "vH2  [ms]"
            << "    " << std::setw(12) << std::left << "vH  [ms]"
            << "    " << std::setw(12) << std::left << "tot  [ms]"
            << "    " << std::setw(12) << std::left << "op  [ms]"
            << '\n';
    report
            << std::setprecision(3)
            << "    " << std::setw(24) << std::left << " "
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * std::get<0>(lbfgs_log)
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * std::get<1>(lbfgs_log)
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * std::get<2>(lbfgs_log)
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * std::get<3>(lbfgs_log)
            << "    " << std::setw(12) << std::left << std::fixed << 1000 *
                                                                     (    std::get<0>(lbfgs_log)
                                                                        + std::get<1>(lbfgs_log)
                                                                        + std::get<2>(lbfgs_log)
                                                                        + std::get<3>(lbfgs_log)
                                                                     )
            << "    " << std::setw(12) << std::left << std::fixed << 1000 * std::get<4>(lbfgs_log)
            << '\n';

    report << '\n';
    tools::log->debug(report.str());

}