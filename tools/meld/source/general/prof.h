#pragma once
#include "meld-io/id.h"
#include <h5pp/details/h5ppHid.h>
#include <memory>
#include <vector>
namespace tools::prof {
    void print_profiling();

    extern std::string get_mem_usage();
    extern void        print_mem_usage();
    extern void        print_mem_usage_oneliner();
    extern double      mem_usage_in_mb(std::string_view name);
    extern double      mem_rss_in_mb();
    extern double      mem_hwm_in_mb();
    extern double      mem_vm_in_mb();

}
