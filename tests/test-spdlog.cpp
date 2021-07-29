#include <io/fmt.h>
#include <tools/common/log.h>
int main (){
    tools::log = tools::Logger::setLogger("test-log",0,true);
    tools::log->info("This is a test log");
    fmt::print("Hello world\n");

    std::vector<bool> test1 = {true,false,true};
    std::vector<double> test2 = {1,2,3};
//    tools::log->info("This tests logging a vector<bool>   {}", test1); // Broken
    tools::log->info("This tests logging a vector<bool>   {}", fmt::join(test1,","));
    tools::log->info("This tests logging a vector<double> {}", test2);
    tools::log->info("This tests logging a vector<double> {:+.3f}", fmt::join(test2,","));
}