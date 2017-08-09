//
// Created by david on 2016-08-14.
//

#ifndef CLASS_TIC_TOC_H
#define CLASS_TIC_TOC_H


#include <chrono>
#include <iostream>

using namespace std;
using namespace chrono;



class class_tic_toc {
private:
    high_resolution_clock::time_point tic_timepoint;
    high_resolution_clock::time_point start_timepoint;
    bool profiling           = false;           //Whether we are profiling_on or not.
    int print_precision      = 5 ;
    string name              = "";
    int padding              = 5;

public:
    class_tic_toc(bool on_off, int prec, string output_text);                 //Constructor
    class_tic_toc(){};
    high_resolution_clock::duration delta_time;
    high_resolution_clock::duration measured_time;
    void set_properties(bool on_off, int prec, string output_text);
    void tic();
    void toc();
    void print_time();
//    void print_time(high_resolution_clock::duration total_runtime);
    void print_time_w_percent();
    void print_delta();
    void print_total_reset();
    void reset();
    friend std::ostream &operator<<(std::ostream &, const class_tic_toc &);
};


#endif //CLASS_TIC_TOC_H