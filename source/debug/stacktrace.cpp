#include "stacktrace.h"
#include <csignal>
#include <cxxabi.h>
#include <stdexcept>

void debug::signal_callback_handler(int status) {
    debug::print_stack_trace();
    switch(status) {
        case SIGTERM: {
            std::fprintf(stderr, "Exit SIGTERM: %d\n", status);
            break;
        }
        case SIGKILL: {
            std::fprintf(stderr, "Exit SIGKILL: %d\n", status);
            std::quick_exit(status); // Does not call destructors
            break;
        }
        case SIGINT: {
            std::fprintf(stderr, "Exit SIGINT: %d\n", status);
            break;
        }
        case SIGHUP: {
            std::fprintf(stderr, "Exit SIGHUP: %d\n", status);
            break;
        }
        case SIGQUIT: {
            std::fprintf(stderr, "Exit SIGQUIT: %d\n", status);
            break;
        }
        case SIGABRT: {
            std::fprintf(stderr, "Exit SIGABRT: %d\n", status);
            break;
        }
        case SIGSEGV: {
            std::fprintf(stderr, "Exit SIGSEGV: %d\n", status);
            std::quick_exit(status); // Does not call destructors
            break;
        }
        default: {
            std::fprintf(stderr, "Exit %d\n", status);
            break;
        }
    }

    std::exit(status);
}


void debug::register_callbacks() {
    signal(SIGTERM, signal_callback_handler);
    signal(SIGKILL, signal_callback_handler);
    signal(SIGINT, signal_callback_handler);
    signal(SIGHUP, signal_callback_handler);
    signal(SIGQUIT, signal_callback_handler);
    signal(SIGABRT, signal_callback_handler);
    signal(SIGSEGV, signal_callback_handler);

    /* ISO C99 signals.  */
    signal(SIGINT, signal_callback_handler);  //	 2	/* Interactive attention signal.  */
    signal(SIGILL, signal_callback_handler);  //	 4	/* Illegal instruction.  */
    signal(SIGABRT, signal_callback_handler); //	 6	/* Abnormal termination.  */
    signal(SIGFPE, signal_callback_handler);  //	 8	/* Erroneous arithmetic operation.  */
    signal(SIGSEGV, signal_callback_handler); //	 11	/* Invalid access to storage.  */
    signal(SIGTERM, signal_callback_handler); //	 15	/* Termination request.  */

    /* Historical signals specified by POSIX. */
    signal(SIGHUP, signal_callback_handler);  //	 1	/* Hangup.  */
    signal(SIGQUIT, signal_callback_handler); //	 3	/* Quit.  */
    signal(SIGTRAP, signal_callback_handler); //	 5	/* Trace/breakpoint trap.  */
    signal(SIGKILL, signal_callback_handler); //	 9	/* Killed.  */
    signal(SIGBUS, signal_callback_handler);  //	 10	/* Bus error.  */
    signal(SIGSYS, signal_callback_handler);  //	 12	/* Bad system call.  */
    signal(SIGPIPE, signal_callback_handler); //	 13	/* Broken pipe.  */
    signal(SIGALRM, signal_callback_handler); //	 14	/* Alarm clock.  */

    /* New(er) POSIX signals (1003.1-2008, 1003.1-2013).  */
    signal(SIGURG, signal_callback_handler);    //   16 /* Urgent data is available at a socket.  */
    signal(SIGSTOP, signal_callback_handler);   //   17 /* Stop, unblockable.  */
    signal(SIGTSTP, signal_callback_handler);   //   18 /* Keyboard stop.  */
    signal(SIGCONT, signal_callback_handler);   //   19 /* Continue.  */
    signal(SIGCHLD, signal_callback_handler);   //   20 /* Child terminated or stopped.  */
    signal(SIGTTIN, signal_callback_handler);   //   21 /* Background read from control terminal.  */
    signal(SIGTTOU, signal_callback_handler);   //   22 /* Background write to control terminal.  */
    signal(SIGPOLL, signal_callback_handler);   //   23 /* Pollable event occurred (System V).  */
    signal(SIGXCPU, signal_callback_handler);   //   24 /* CPU time limit exceeded.  */
    signal(SIGXFSZ, signal_callback_handler);   //   25 /* File size limit exceeded.  */
    signal(SIGVTALRM, signal_callback_handler); //   26 /* Virtual timer expired.  */
    signal(SIGPROF, signal_callback_handler);   //   27 /* Profiling timer expired.  */
    signal(SIGUSR1, signal_callback_handler);   //   30 /* User-defined signal 1.  */
    signal(SIGUSR2, signal_callback_handler);   //   31 /* User-defined signal 2.  */
    signal(137, signal_callback_handler);   //   137 /* LLDB?  */
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