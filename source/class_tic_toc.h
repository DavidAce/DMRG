//
// Created by david on 2016-08-14.
//

#ifndef CLASS_TIC_TOC_H
#define CLASS_TIC_TOC_H


#include <chrono>
#include <iostream>

using namespace std;
using namespace chrono;
class class_profiling {
private:
    high_resolution_clock::time_point delta_tic;
    high_resolution_clock::time_point delta_toc;
    int profiling;           //Whether we are profiling or not.
    int print_precision;
    int padding = 5;
    string name;
public:
    class_profiling(int on_off, int prec, string output_text);                 //Constructor
    high_resolution_clock::duration delta_time;
    high_resolution_clock::duration total_time;

    void tic();
    void toc();
    void print_total();
    void print_total(high_resolution_clock::duration total_runtime);
    void print_delta();
    void print_total_reset();
    void reset();
    friend std::ostream &operator<<(std::ostream &, const class_profiling &);
};


#endif //CLASS_TIC_TOC_H