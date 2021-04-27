#pragma once

#include <stdexcept>

namespace except {
    class state_error : public std::runtime_error {
        // Used for signaling that no resumable state was found
        using std::runtime_error::runtime_error;
    };

    class file_error : public std::runtime_error {
        // Used for signaling that the existing file is corrupted
        using std::runtime_error::runtime_error;
    };

    class load_error : public std::runtime_error {
        // Used for signaling an error when loading an existing file
        using std::runtime_error::runtime_error;
    };

    class resume_error : public std::runtime_error {
        // Used to signal that an error ocurred when trying to resume a simulation
        using std::runtime_error::runtime_error;
    };
}
