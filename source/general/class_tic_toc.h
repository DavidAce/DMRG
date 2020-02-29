//
// Created by david on 2016-08-14.
//

#pragma once

#include <chrono>
#include <iostream>
#include <string>


class class_tic_toc {
private:
    std::chrono::high_resolution_clock::time_point tic_timepoint;
    std::chrono::high_resolution_clock::time_point start_timepoint;
    bool profiling           = false;           //Whether we are profiling or not.
    int print_precision      = 5 ;
    std::string name         = "";
    int padding              = 5;

public:
    class_tic_toc(bool on_off,int prec, std::string output_text);                 //Constructor
    class_tic_toc()= default;
    std::chrono::high_resolution_clock::duration delta_time;
    std::chrono::high_resolution_clock::duration measured_time;
    void set_properties(bool on_off,int prec, std::string output_text);
    void set_label(std::string output_text);
    std::string  get_name();
    double get_measured_time();


    void tic();
    void toc();
    void print_age();
    void print_time();
    double get_age();
    double get_last_time_interval();
//    void print_time(high_resolution_clock::duration total_runtime);
    void print_time_w_percent();
    void print_time_w_percent(class_tic_toc &parent);
    void print_time_w_percent_if_nonzero(class_tic_toc &parent);
    void print_delta();
    void print_total_reset();
    void reset();
    friend std::ostream &operator<<(std::ostream &, const class_tic_toc &);
};

