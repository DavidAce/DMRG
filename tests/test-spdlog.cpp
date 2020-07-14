#include <tools/common/log.h>
int main (){
    tools::log = tools::Logger::setLogger("test-log",0,true);
    tools::log->info("This is a test log");
    fmt::print("Hello world");
}