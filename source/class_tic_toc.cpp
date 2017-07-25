//
// Created by david on 2016-08-14.
//

#include <class_tic_toc.h>
#include <iomanip>

class_profiling::class_profiling(int on_off, int prec, string output_text)
        : profiling(on_off),
          print_precision(prec),
          name(output_text)
{
    if (profiling) {
        if (!name.empty()){
            name = " " + name + ": ";
        }
        total_time = total_time.zero();
        delta_time = delta_time.zero();
    }
}

void class_profiling::tic(){
    if (profiling) {
        delta_tic = high_resolution_clock::now();
    }
}

void class_profiling::toc(){
    if (profiling) {
        delta_toc   = high_resolution_clock::now();
        delta_time  = delta_toc - delta_tic;
        total_time += delta_time;
    }
}

void class_profiling::print_delta(){
    if (profiling) {
        cout << setprecision(print_precision) << fixed << setw(print_precision + padding)
             << name
             << duration_cast<duration<double>>(delta_time).count();
    }
}

void class_profiling::print_total(){
    if (profiling) {
        cout << setprecision(print_precision)  << fixed << setw(print_precision + padding)
             << name
             << duration_cast<duration<double>>(total_time).count();
    }
}

void class_profiling::print_total(high_resolution_clock::duration total_runtime){
    if (profiling) {
        cout << setprecision(print_precision)  << fixed << setw(print_precision + padding)
             << name
             << duration_cast<duration<double>>(total_time).count()
             <<"       | " << 100.0*total_time.count() / total_runtime.count() << " %";
    }
}

void class_profiling::print_total_reset(){
    if (profiling) {
        cout << setprecision(print_precision)  << fixed << setw(print_precision + padding)
             << name
             << duration_cast<duration<double>>(total_time).count();
        reset();
    }
}

void class_profiling::reset() {
    if (profiling) {
        total_time = total_time.zero();
    }
}


std::ostream &operator<<(std::ostream &os, const class_profiling &t) {
    if (t.profiling) {
        os  << setprecision(t.print_precision)  << fixed << setw(t.print_precision + t.padding)
            << t.name
            << t.total_time.count();
    }
    return os;
}
