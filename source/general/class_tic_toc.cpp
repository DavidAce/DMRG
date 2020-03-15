//
// Created by david on 2016-08-14.
//

#include "class_tic_toc.h"
#include <iomanip>


class_tic_toc::class_tic_toc(bool on_off, int prec, std::string output_text)
        : name(output_text),
          enable(on_off),
          print_precision(prec)
{
    if (enable) {
        if (!name.empty()){
            name = name + ": ";
        }
        reset();
    }
}

void class_tic_toc::set_properties(bool on_off,int prec, std::string output_text){
    *this = class_tic_toc(on_off, prec, output_text);
}

void class_tic_toc::set_label(std::string output_text){
    *this = class_tic_toc(enable, print_precision, output_text);
}

void class_tic_toc::set_time(double other_time_in_seconds) {
    measured_time = std::chrono::duration_cast<hresclock::duration>(std::chrono::duration<double>(other_time_in_seconds));
}

std::string class_tic_toc::get_name() const {
    return name;
}

double class_tic_toc::get_measured_time() const {
    return std::chrono::duration_cast<std::chrono::duration<double>>(measured_time).count();
}

double class_tic_toc::get_measured_time_and_reset(){
    auto result = std::chrono::duration_cast<std::chrono::duration<double>>(measured_time).count();
    reset();
    return result;
}

void class_tic_toc::tic(){
    if (enable) tic_timepoint = std::chrono::high_resolution_clock::now();
}

void class_tic_toc::toc(){
    if (enable) {
        delta_time     = std::chrono::high_resolution_clock::now() - tic_timepoint;
        measured_time += delta_time;
    }
}

void class_tic_toc::print_last_time_interval() const {
    if (enable) {
        std::cout << std::setprecision(print_precision) << std::fixed << std::setw(print_precision + padding)
             << name
             << std::chrono::duration_cast<std::chrono::duration<double>>(delta_time).count();
    }
}

void class_tic_toc::print_age() const {
    if (enable) {
        std::cout
            << name
            << std::fixed << std::setprecision(print_precision) << std::setw(print_precision + padding) << std::left
            << get_age() << " s\n";
    }
}


void class_tic_toc::print_measured_time() const {
    if (enable) {
        std::cout
            << name
            << std::fixed << std::setprecision(print_precision)   << std::setw(print_precision + padding) << std::left
            << std::chrono::duration_cast<std::chrono::duration<double>>(measured_time).count() << " s\n";
    }
}

void class_tic_toc::print_measured_time_and_reset(){
    if (enable) {
        std::cout
            << name
            << std::fixed << std::setprecision(print_precision)   << std::setw(print_precision + padding) << std::left
            << std::chrono::duration_cast<std::chrono::duration<double>>(measured_time).count() << " s\n";
        reset();
    }
}

void class_tic_toc::print_measured_time_w_percent() const{
    if (enable) {
        std::cout << name
             << std::fixed << std::setprecision(print_precision) << std::setw(print_precision + padding) << std::left
             << std::chrono::duration_cast<std::chrono::duration<double>>(measured_time).count() << " s |"
             << std::fixed << std::setprecision(print_precision) << std::setw(print_precision + padding) << std::right
             << 100.0 << " %\n";
    }
}

void class_tic_toc::print_measured_time_w_percent(class_tic_toc &parent) const {
    if (enable) {
        std::cout << name
             << std::fixed << std::setprecision(print_precision) << std::setw(print_precision + padding) << std::left
             << std::chrono::duration_cast<std::chrono::duration<double>>(measured_time).count() << " s |"
             << std::fixed << std::setprecision(print_precision) << std::setw(print_precision + padding) << std::right
             << 100.0*measured_time.count() / parent.measured_time.count() << " % \n";
    }
}

void class_tic_toc::print_measured_time_w_percent_if_nonzero(class_tic_toc &parent) const {
    if (measured_time.count() > 0.0){
        print_measured_time_w_percent(parent);
    }
}



double class_tic_toc::get_age() const {
    return std::chrono::duration_cast<std::chrono::duration<double>> (std::chrono::high_resolution_clock::now() - start_timepoint).count();
}

double class_tic_toc::get_last_time_interval() const {
    return std::chrono::duration_cast<std::chrono::duration<double>> (delta_time).count();
}

void class_tic_toc::reset() {
    if (enable) {
        measured_time   = std::chrono::high_resolution_clock::duration::zero();
        delta_time      = std::chrono::high_resolution_clock::duration::zero();
        start_timepoint = std::chrono::high_resolution_clock::now();
    }
}


class_tic_toc & class_tic_toc::operator=(double other_time_in_seconds) {
    this->measured_time = std::chrono::duration_cast<hresclock::duration>(std::chrono::duration<double>(other_time_in_seconds));
    return *this;
}

std::ostream &operator<<(std::ostream &os, const class_tic_toc &t) {
    if (t.enable) {
        os  << std::setprecision(t.print_precision)  << std::fixed << std::setw(t.print_precision + t.padding)
            << t.name
            << t.measured_time.count();
    }
    return os;
}
