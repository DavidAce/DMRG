#pragma once
#include <string>
#include <string_view>

namespace debug {
    extern std::string hostname();
    extern std::string cpu_info(std::string_view info = "model name");
    extern void        print_mem_usage();
    extern double      mem_usage_in_mb(std::string_view name);
    extern double      mem_rss_in_mb();
    extern double      mem_hwm_in_mb();
    extern double      mem_vm_in_mb();
}