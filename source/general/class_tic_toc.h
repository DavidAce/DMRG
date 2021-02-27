//
// Created by david on 2016-08-14.
//

#pragma once

#include <chrono>
#include <string>

class class_tic_toc {
    private:
    using hresclock = std::chrono::high_resolution_clock;
    hresclock::time_point tic_timepoint;
    hresclock::time_point toc_timepoint;
    hresclock::time_point lap_timepoint;
    hresclock::time_point start_timepoint;
    hresclock::duration   delta_time    = std::chrono::high_resolution_clock::duration::zero();
    hresclock::duration   measured_time = std::chrono::high_resolution_clock::duration::zero();
    hresclock::duration   lap_time;

    std::string name;
    bool        enable          = true; // Whether we are profiling or not.
    int         print_precision = 5;

    public:
    bool is_measuring = false;
    class_tic_toc(bool on_off, int prec, std::string output_text); // Constructor
    class_tic_toc();
    void                      tic();
    void                      toc();

    void                      reset();
    void                      set_properties(bool on_off, int prec, std::string output_text);
    void                      set_label(std::string output_text);
    void                      set_measured_time(double measured_time);
    void                      start_lap();
    [[nodiscard]] std::string get_name() const;
    [[nodiscard]] double      get_age() const;
    [[nodiscard]] double      get_measured_time() const;
    [[nodiscard]] double      get_last_interval() const;
    [[nodiscard]] double      restart_lap();
    [[nodiscard]] double      get_lap() const;

    class_tic_toc &      operator=(double other_time_in_seconds);
    class_tic_toc &      operator+=(double other_time_in_seconds);
    class_tic_toc &      operator-=(double other_time_in_seconds);
    class_tic_toc &      operator+=(const class_tic_toc &rhs);
    class_tic_toc &      operator-=(const class_tic_toc &rhs);

    class token {
        class_tic_toc & t;
        public:
        explicit token(class_tic_toc & t_);
        ~token() noexcept;
        void tic();
        void toc();
    };
    [[nodiscard]] token tic_token(); // Gives a token RAII-style tic-toc: toc() is called by destructor!

};
