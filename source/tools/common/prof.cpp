#include "prof.h"
#include "log.h"
#include <algorithm>
#include <algorithms/AlgorithmStatus.h>
#include <config/enums.h>
#include <config/settings.h>
#include <fstream>
#include <sstream>
#include <tid/tid.h>

void tools::common::profile::print_profiling(const AlgorithmStatus &status) {
    if(settings::profiling::on) {
//        static double last_print_time = 0;
        auto          t_tot           = tid::get_unscoped("t_tot");
//        if(std::abs(t_tot.get_time() - last_print_time) < 5.0) return; // Do not print if there's already been a print within the last 5 seconds

        for(const auto &t : tid::get_tree(status.algo_type_sv())) {
            if(t->get_level() <= tid::level::pedant) tools::log->info("{}", t.str());
        }
//        last_print_time = t_tot.get_time();
    }
}

void tools::common::profile::print_profiling() {
    if(settings::profiling::on) {
//        static double last_print_time = 0;
        auto          t_tot           = tid::get_unscoped("t_tot");
//        if(std::abs(t_tot.get_time() - last_print_time) < 5.0) return; // Do not print if there's already been a print within the last 5 seconds

        for(const auto &t : tid::get_tree())
            if(t->get_level() <= tid::level::pedant) tools::log->info("{}", t.str());

//        last_print_time = t_tot.get_time();
    }
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
