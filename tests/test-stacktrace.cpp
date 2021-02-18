
#include <general/stack_trace.h>
#include <stdexcept>
#include <iostream>

int main() {
    // Register termination codes and what to do in those cases
    debug::register_callbacks();
    try{
        throw std::runtime_error("Testing");
    }catch (const std::exception & ex){
        std::cout << "Caught error: " << ex.what() << std::endl;
    }
    return 0;
}