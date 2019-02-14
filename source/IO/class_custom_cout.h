//
// Created by david on 2018-01-16.
//

#ifndef DMRG_CLASS_CONSOLE_PRINTER_H
#define DMRG_CLASS_CONSOLE_PRINTER_H

#include <chrono>
#include <iostream>

class class_custom_cout
{
public:
    int verbosity = 0;
    int temp_verbosity = 0;
    std::ostream &os;

    explicit class_custom_cout(std::ostream &o = std::cout):os(o){};
    explicit class_custom_cout(int verbosity_, std::ostream &o = std::cout):verbosity(verbosity_),temp_verbosity(verbosity_),os(o){};

    void set_verbosity(int new_verbosity){verbosity = new_verbosity; temp_verbosity = new_verbosity;}


    template <typename T>
    class_custom_cout &operator<<(const T &a) {
        if (temp_verbosity >= verbosity) {
            os << a;
        }else{
            os.clear();
        }
        return *this;
    }

    class_custom_cout &operator<<(std::ostream& (*pf) (std::ostream&)){
        if (temp_verbosity <= verbosity) {
            os << pf;
        }else{
            os.clear();
        }
        temp_verbosity = verbosity;
        return *this;
    }
    class_custom_cout& operator()(int level) {
        temp_verbosity = level;
        return *this;
    }
};


#endif //DMRG_CLASS_CONSOLE_PRINTER_H
