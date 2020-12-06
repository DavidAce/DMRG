#pragma once

#include <stdexcept>

namespace ex{
    class state_error: public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    class file_error: public std::runtime_error {
        using std::runtime_error::runtime_error;
    };
}

