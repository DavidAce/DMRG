//
// Created by david on 2018-01-16.
//

#include "class_custom_cout.h"
class_custom_cout&  class_custom_cout::operator<<( std::ostream&(*f)(std::ostream&) ) {
    if(settings::console::timestamp){
        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        timestamp = std::ctime(&end_time);
        timestamp = timestamp.substr(0, timestamp.size()-1) + ": ";

    }else{
        timestamp = "";
    }
    if (settings::console::verbosity >= verbosity){
        std::cout << timestamp << oss.str() << f;
    }
    oss.str("");

    verbosity = settings::console::verbosity;
    return *this;
}


class_custom_cout& class_custom_cout::operator()(int level) {
    verbosity = level;
    return *this;
}