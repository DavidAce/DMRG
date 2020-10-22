
#include "stack_trace.h"
#include <stdexcept>

#if __has_include(<backward.hpp2>)
    #include "backward.hpp"
namespace backward {
    backward::SignalHandling sh;
}

void debug::print_stack_trace() {
    using namespace backward;
    StackTrace st;
    st.load_here(32);
    Printer p;
    p.object     = true;
    p.color_mode = ColorMode::always;
    p.address    = true;
    p.print(st, stderr);
}

#elif __has_include(<libunwind.h>)
    #define UNW_LOCAL_ONLY

    #include <cstdio>
    #include <cxxabi.h>
    #include <execinfo.h>
    #include <libunwind.h>
    #include <unistd.h>
void debug::print_stack_trace() {
    unw_cursor_t  cursor;
    unw_context_t context;

    // Initialize cursor to current frame for local unwinding.
    unw_getcontext(&context);
    unw_init_local(&cursor, &context);

    // Unwind frames one by one, going up the frame stack.
    while(unw_step(&cursor) > 0) {
        unw_word_t offset, pc;
        unw_get_reg(&cursor, UNW_REG_IP, &pc);
        if(pc == 0) { break; }
        std::printf("0x%-10lx", pc);

        char sym[256];
        if(unw_get_proc_name(&cursor, sym, sizeof(sym), &offset) == 0) {
            char *nameptr = sym;
            int   status;
            char *demangled = abi::__cxa_demangle(sym, nullptr, nullptr, &status);
            if(status == 0) { nameptr = demangled; }
            std::printf(" +0x%-10lx: %s\n", offset,nameptr);
            std::free(demangled);
        } else {
            std::printf(" -- error: unable to obtain symbol name for this frame\n");
        }
    }

    // get void*'s for all entries on the stack
    void * array[10];
    size_t size = backtrace(array, 10);

    // print out all the frames to stderr
    backtrace_symbols_fd(array, size, STDERR_FILENO);
}

#else
void debug::print_stack_trace() {}
#endif

void debug::throw_stack_trace(const std::string &msg) {
    print_stack_trace();
    throw std::runtime_error(msg);
}