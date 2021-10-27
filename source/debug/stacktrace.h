#pragma once

namespace debug {
    void register_callbacks();
    void signal_callback_handler(int signum);
    void print_stack_trace();
}