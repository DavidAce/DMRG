#include "stacktrace.h"
#include <csignal>
#include <cxxabi.h>
#include <stdexcept>

void debug::signal_callback_handler(int signum) {
    switch(signum) {
        case SIGTERM: {
            std::fprintf(stderr, "Caught SIGTERM\n");
            break;
        }
        case SIGKILL: {
            std::fprintf(stderr, "Caught SIGKILL\n");
            break;
        }
        case SIGINT: {
            std::fprintf(stderr, "Caught SIGINT\n");
            break;
        }
        case SIGHUP: {
            std::fprintf(stderr, "Caught SIGHUP\n");
            break;
        }
        case SIGQUIT: {
            std::fprintf(stderr, "Caught SIGQUIT\n");
            break;
        }
        case SIGABRT: {
            std::fprintf(stderr, "Caught SIGABRT\n");
            break;
        }
        case SIGSEGV: {
            std::fprintf(stderr, "Caught SIGSEGV\n");
            break;
        }
        default: break;
    }
    debug::print_stack_trace();
    std::fprintf(stderr, "Exiting\n");
    std::quick_exit(signum);
}

void debug::register_callbacks() {
    signal(SIGTERM, signal_callback_handler);
    signal(SIGKILL, signal_callback_handler);
    signal(SIGINT, signal_callback_handler);
    signal(SIGHUP, signal_callback_handler);
    signal(SIGQUIT, signal_callback_handler);
    signal(SIGABRT, signal_callback_handler);
    signal(SIGSEGV, signal_callback_handler);
}

#if __has_include(<libunwind.h>) && defined(DMRG_HAS_UNWIND)

    #define UNW_LOCAL_ONLY
    #include <array>
    #include <climits>
    #include <cstdio>
    #include <cstdlib>
    #include <execinfo.h>
    #include <libunwind.h>
    #include <memory>
    #include <unistd.h>
    #include <vector>
std::string getexepath() {
    static char result[PATH_MAX];
    ssize_t     count = readlink("/proc/self/exe", result, PATH_MAX);
    return std::string(result, (count > 0) ? static_cast<unsigned long>(count) : 0);
}

std::string sh(std::string_view cmd) {
    std::array<char, 256> buffer{};
    std::string           result;
    std::shared_ptr<FILE> pipe(popen(cmd.data(), "r"), pclose);
    if(!pipe) return std::string();
    while(!feof(pipe.get())) {
        if(fgets(buffer.data(), 256, pipe.get()) != nullptr) { result += buffer.data(); }
    }
    return result;
}

std::string getErrorLocation(unw_word_t addr) {
    // our program need to be passed after the -e parameter
    auto exec_path = getexepath();
    // Generate a hex address for addr2line
    char hex_addr[16];
    snprintf(hex_addr, 16, "0x%lx", addr);

    // Generate the command and call it
    std::string cmd       = "addr2line --exe " + exec_path + " --functions --demangle " + hex_addr;
    auto        locString = sh(cmd);
    if(locString.empty()) return locString;

    // Format the string
    size_t start_pos = 0;
    while((start_pos = locString.find('\n', start_pos)) != std::string::npos) {
        if(start_pos >= locString.size() - 1)
            locString.replace(start_pos, 1, "");
        else
            locString.replace(start_pos, 1, " in file ");
        start_pos += 1;
    }
    if(locString.find("??") != std::string::npos)
        return std::string();
    else
        return locString;
}

void debug::print_stack_trace() {
    unw_cursor_t  cursor;
    unw_context_t context;
    // Initialize cursor to current frame for local unwinding.
    unw_getcontext(&context);
    unw_init_local(&cursor, &context);
    int count = 0;
    // Unwind frames one by one, going up the frame stack.
    while(unw_step(&cursor) > 0) {
        unw_word_t offset, ip;
        unw_get_reg(&cursor, UNW_REG_IP, &ip);
        if(ip == 0) { break; }
        char sym[1024];
        if(unw_get_proc_name(&cursor, sym, sizeof(sym), &offset) == 0) {
            std::fprintf(stderr, "#%-3d 0x%-10lx", count++, ip);
            auto estr = getErrorLocation(ip);
            if(not estr.empty())
                std::fprintf(stderr, " %s\n", estr.c_str());
            else {
                int   status;
                char *nameptr   = sym;
                char *demangled = abi::__cxa_demangle(sym, nullptr, nullptr, &status);
                if(status == 0) { nameptr = demangled; }
                std::fprintf(stderr, " %s\n", nameptr);
                std::free(demangled);
            }

        } else {
            std::fprintf(stderr, " -- error: unable to obtain symbol name for this frame\n");
        }
    }
}

#else
void debug::print_stack_trace() {}
#endif

void debug::throw_stack_trace(const std::string & msg) {
    print_stack_trace();
    throw std::runtime_error(msg);
}