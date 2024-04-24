
#pragma once

#include "logger.h"
#include <h5pp/h5pp.h>

namespace tools::h5dbg {
    std::string get_hid_string_details(const hid_t &id);

    void print_dangling_ids(const h5pp::File &file, std::string_view func, int line, unsigned obj = H5F_OBJ_FILE);
    void assert_no_dangling_ids(const h5pp::File &file, std::string_view func, int line, unsigned obj = H5F_OBJ_FILE);
}