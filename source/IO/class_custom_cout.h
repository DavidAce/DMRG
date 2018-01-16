//
// Created by david on 2018-01-16.
//

#ifndef DMRG_CLASS_CONSOLE_PRINTER_H
#define DMRG_CLASS_CONSOLE_PRINTER_H

#include <chrono>
#include <sim_parameters/n_sim_settings.h>
#include <sstream>
#include <iostream>

class class_custom_cout {
private:
    std::ostringstream oss;
    std::string timestamp = "";
    int verbosity = settings::console::verbosity;
public:
    class_custom_cout& operator<<( std::ostream&(*f)(std::ostream&) );

    class_custom_cout& operator()(int level);
    // this overload receive the single values to append via <<
    template<typename T>
    class_custom_cout& operator<<(T&& value) {
        oss << value;
        return *this;
    }
};



#endif //DMRG_CLASS_CONSOLE_PRINTER_H
