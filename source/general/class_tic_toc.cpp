//
// Created by david on 2016-08-14.
//

#include "class_tic_toc.h"
#include <iomanip>

class_tic_toc::class_tic_toc(bool on_off, int prec, string output_text)
        : profiling(on_off),
          print_precision(prec),
          name(output_text)
{
    if (profiling) {
        if (!name.empty()){
            name = " " + name + ": ";
        }
        measured_time   = measured_time.zero();
        delta_time      = delta_time.zero();
        start_timepoint = high_resolution_clock::now();
    }
}

void class_tic_toc::set_properties(bool on_off, int prec, string output_text){
    *this = class_tic_toc(on_off, prec, output_text);
}

void class_tic_toc::tic(){
    if (profiling) {
        tic_timepoint = high_resolution_clock::now();
    }
}

void class_tic_toc::toc(){
    if (profiling) {
        delta_time       = high_resolution_clock::now() - tic_timepoint;
        measured_time   += delta_time;
    }
}

void class_tic_toc::print_delta(){
    if (profiling) {
        cout << setprecision(print_precision) << fixed << setw(print_precision + padding)
             << name
             << duration_cast<duration<double>>(delta_time).count();
    }
}

void class_tic_toc::print_time(){
    if (profiling) {
        cout << name
             << fixed << setprecision(print_precision) << setw(print_precision + padding) << left
             << duration_cast<duration<double>>(measured_time).count()<< " s\n";
    }
}


void class_tic_toc::print_time_w_percent(){
    if (profiling) {
        cout << name
             << fixed << setprecision(print_precision) << setw(print_precision + padding) << left
             << duration_cast<duration<double>>(measured_time).count() << " s |"
             << fixed << setprecision(print_precision) << setw(print_precision + padding) << right
             << 100.0*measured_time.count() / (high_resolution_clock::now() - start_timepoint).count() << " %\n";
    }
}

void class_tic_toc::print_total_reset(){
    if (profiling) {
        cout << setprecision(print_precision)  << fixed << setw(print_precision + padding)
             << name
             << duration_cast<duration<double>>(measured_time).count() << " s\n";
        reset();
    }
}

double class_tic_toc::get_age() {
    return duration_cast<duration<double>> (high_resolution_clock::now() - start_timepoint).count();
}

void class_tic_toc::reset() {
    if (profiling) {
        measured_time = measured_time.zero();
    }
}


std::ostream &operator<<(std::ostream &os, const class_tic_toc &t) {
    if (t.profiling) {
        os  << setprecision(t.print_precision)  << fixed << setw(t.print_precision + t.padding)
            << t.name
            << t.measured_time.count();
    }
    return os;
}
