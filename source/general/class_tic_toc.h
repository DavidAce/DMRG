//
// Created by david on 2016-08-14.
//

#pragma once

#include <chrono>
#include <iostream>
#include <string>

class class_tic_toc {
    private:
    using hresclock = std::chrono::high_resolution_clock;
    hresclock::time_point tic_timepoint;
    hresclock::time_point start_timepoint;

    std::string name            = "";
    bool        enable          = false; // Whether we are profiling or not.
    int         print_precision = 5;
    int         padding         = 5;

    public:
    class_tic_toc(bool on_off, int prec, std::string output_text); // Constructor
    class_tic_toc() = default;
    void tic();
    void toc();

    hresclock::duration delta_time;
    hresclock::duration measured_time;

    void        set_properties(bool on_off, int prec, std::string output_text);
    void        set_label(std::string output_text);
    void        set_time(double measured_time);
    std::string get_name() const;
    double      get_age() const;
    double      get_last_time_interval() const ;
    double      get_measured_time() const;
    double      get_measured_time_and_reset();

    void print_age() const;
    void print_measured_time() const;
    void print_measured_time_and_reset();
    void print_last_time_interval() const;
    void print_measured_time_w_percent() const;
    void print_measured_time_w_percent(class_tic_toc &parent) const;
    void print_measured_time_w_percent_if_nonzero(class_tic_toc &parent) const;
    void reset();

    class_tic_toc & operator= (double measured_time);
    friend std::ostream &operator<<(std::ostream &, const class_tic_toc &);
};
