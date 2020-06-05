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
    hresclock::time_point toc_timepoint;
    hresclock::time_point start_timepoint;
    hresclock::duration   delta_time    = std::chrono::high_resolution_clock::duration::zero();
    hresclock::duration   measured_time = std::chrono::high_resolution_clock::duration::zero();

    std::string name            = "";
    bool        enable          = true; // Whether we are profiling or not.
    bool        is_measuring    = false;
    int         print_precision = 5;
    int         padding         = 5;

    public:
    class_tic_toc(bool on_off, int prec, std::string output_text); // Constructor
    class_tic_toc();
    void tic();
    void toc();

    void                      set_properties(bool on_off, int prec, std::string output_text);
    void                      set_label(std::string output_text);
    void                      set_time(double measured_time);
    void                      print_age() const;
    void                      print_measured_time() const;
    void                      print_last_time_interval() const;
    void                      print_measured_time_w_percent(double cmp = std::numeric_limits<double>::quiet_NaN()) const;
    [[nodiscard]] std::string get_name() const;
    [[nodiscard]] double      get_age() const;
    [[nodiscard]] double      get_measured_time() const;
    [[nodiscard]] double      get_last_time_interval() const;
    [[nodiscard]] std::string string(double tgt = std::numeric_limits<double>::quiet_NaN(), double cmp = std::numeric_limits<double>::quiet_NaN()) const;
    [[nodiscard]] std::string string_age() const;
    [[nodiscard]] std::string string_measured_time() const;
    [[nodiscard]] std::string string_last_time_interval() const;
    [[nodiscard]] std::string string_measured_time_w_percent(double cmp = std::numeric_limits<double>::quiet_NaN()) const;

    void reset();

    class_tic_toc &      operator=(double other_time_in_seconds);
    class_tic_toc &      operator+=(double other_time_in_seconds);
    class_tic_toc &      operator-=(double other_time_in_seconds);
    class_tic_toc &      operator+=(const class_tic_toc & rhs);
    class_tic_toc &      operator-=(const class_tic_toc & rhs);
    friend std::ostream &operator<<(std::ostream &, const class_tic_toc &);
};
