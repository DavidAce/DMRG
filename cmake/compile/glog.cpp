#include <glog/logging.h>
int main(int argc, char* argv[]) {
    // Initialize Google's logging library.
    google::InitGoogleLogging(argv[0]);
    LOG(INFO) << "Found " << 0 << " cookies";
    return 0;
}