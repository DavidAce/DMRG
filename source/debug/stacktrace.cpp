#include "stacktrace.h"
#include <csignal>
#include <cstdio>
#include <cstdlib>

void debug::signal_callback_handler(int status) {
    switch(status) {
        case SIGTERM: {
            //            debug::print_stack_trace();
            std::fprintf(stderr, "Exit SIGTERM: %d\n", status);
            break;
        }
        case SIGSTOP: {
            std::fprintf(stderr, "Exit SIGSTOP: %d\n", status);
            break;
        }
        case SIGCHLD: {
            std::fprintf(stderr, "Exit SIGCHLD: %d\n", status);
            break;
        }
        case SIGKILL: {
            std::fprintf(stderr, "Exit SIGKILL: %d\n", status);
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
            //            debug::print_stack_trace();
            std::fprintf(stderr, "Exit SIGABRT: %d\n", status);
            break;
        }
        case SIGILL: {
            //            debug::print_stack_trace();
            std::fprintf(stderr, "Exit SIGILL: %d\n", status);
            break;
        }
        case SIGFPE: {
            //            debug::print_stack_trace();
            std::fprintf(stderr, "Exit SIGFPE: %d\n", status);
            break;
        }
        case SIGSEGV: {
            //            debug::print_stack_trace();
            std::fprintf(stderr, "Exit SIGSEGV: %d\n", status);
            break;
        }
        default: {
            //            debug::print_stack_trace();
            std::fprintf(stderr, "Exit %d\n", status);
            break;
        }
    }
    if(status != 0) debug::print_stack_trace();
    std::exit(status);
}

void debug::register_callbacks() {
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
    signal(SIGKILL, signal_callback_handler); //	 9	/* Killed. This one is special and can't be caught  */
    signal(SIGBUS, signal_callback_handler);  //	 10	/* Bus error.  */
    signal(SIGSYS, signal_callback_handler);  //	 12	/* Bad system call.  */
    signal(SIGPIPE, signal_callback_handler); //	 13	/* Broken pipe.  */
    signal(SIGALRM, signal_callback_handler); //	 14	/* Alarm clock.  */

    /* New(er) POSIX signals (1003.1-2008, 1003.1-2013).  */
    signal(SIGURG, signal_callback_handler);    //   16 /* Urgent data is available at a socket.  */
    signal(SIGSTOP, signal_callback_handler);   //   17 /* Stop, unblockable.  This one is special and can't be caught  */
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
    signal(137, signal_callback_handler);       //   137 /* LLDB?  */
}

#if __has_include(<backward.hpp>)
    #include <backward.hpp>
void debug::print_stack_trace() {
    backward::StackTrace st;
    st.load_here(128);
    // Skip this scope (1) , as well as the signal_callback_handler scope (2)
    st.skip_n_firsts(2);
    backward::Printer p;
    p.print(st);
}
#else
void debug::print_stack_trace() {}
#endif
