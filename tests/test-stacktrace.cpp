
#include <general/stack_trace.h>
#include <stdexcept>


int main() {
    // Register termination codes and what to do in those cases
    debug::register_callbacks();
    throw std::runtime_error("Testing");
}