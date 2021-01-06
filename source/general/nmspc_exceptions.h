#pragma once

#include <stdexcept>

namespace except {
    class state_error: public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    class file_error: public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    class load_error: public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    class resume_error: public std::runtime_error {
        using std::runtime_error::runtime_error;
    };
}

