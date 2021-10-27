#pragma once

#include <io/fmt.h>
#include <stdexcept>

namespace except {

    template<typename... Args>
    [[nodiscard]] std::runtime_error runtime_error(std::string_view s, Args... args) {
        if(s.rfind("dmrg++: ", 0) == 0)
            return std::runtime_error(fmt::format(s, std::forward<Args>(args)...));
        else
            return std::runtime_error("dmrg++: " + fmt::format(s, std::forward<Args>(args)...));
    }
    template<typename... Args>
    [[nodiscard]] std::logic_error logic_error(std::string_view s, Args... args) {
        if(s.rfind("dmrg++: ", 0) == 0)
            return std::logic_error(fmt::format(s, std::forward<Args>(args)...));
        else
            return std::logic_error("dmrg++: " + fmt::format(s, std::forward<Args>(args)...));
    }
    //    template<typename... Args>
    //    std::runtime_error runtime_error(Args... args) {
    //        return std::runtime_error("dmrg++: " + fmt::format(std::forward<Args>(args)...));
    //    }
    //    template<typename... Args>
    //    std::logic_error logic_error(Args... args) {
    //        return std::logic_error("dmrg++: " + h5pp::format(std::forward<Args>(args)...));
    //    }

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
