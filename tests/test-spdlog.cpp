#include <h5pp/details/h5ppFormat.h>

#include <h5pp/details/h5ppLogger.h>

int main (){
    fmt::print("Hello world");
    std::runtime_error test ("");
    std::runtime_error::operator=(test);
}