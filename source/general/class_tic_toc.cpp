//
// Created by david on 2016-08-14.
//

#include "class_tic_toc.h"
#include <stdexcept>

class_tic_toc::class_tic_toc() : class_tic_toc(true, 5, "") {}

class_tic_toc::class_tic_toc(bool on_off, int prec, std::string output_text) : name(std::move(output_text)), enable(on_off), print_precision(prec) {
    if(enable) {
        if(!name.empty()) { name += ": "; }
        reset();
    }
}

void class_tic_toc::tic() {
    if(enable) {
        if(is_measuring) throw std::runtime_error("Called tic() twice: this timer is already measuring: " + name);
        tic_timepoint = std::chrono::high_resolution_clock::now();
        is_measuring  = true;
    }
}

void class_tic_toc::toc() {
    if(enable) {
        if(not is_measuring) throw std::runtime_error("Called toc() twice or without prior tic()");
        toc_timepoint = std::chrono::high_resolution_clock::now();
        delta_time    = toc_timepoint - tic_timepoint;
        measured_time += delta_time;
        lap_time += delta_time;
        is_measuring = false;
    }
}

class_tic_toc::token class_tic_toc::tic_token() {return class_tic_toc::token(*this);}


void class_tic_toc::set_properties(bool on_off, int prec, std::string output_text) { *this = class_tic_toc(on_off, prec, std::move(output_text)); }

void class_tic_toc::set_label(std::string output_text) { *this = class_tic_toc(enable, print_precision, std::move(output_text)); }

void class_tic_toc::set_measured_time(double other_time_in_seconds) {
    measured_time = std::chrono::duration_cast<hresclock::duration>(std::chrono::duration<double>(other_time_in_seconds));
}

void class_tic_toc::start_lap() {
    lap_time      = std::chrono::high_resolution_clock::duration::zero();
    lap_timepoint = std::chrono::high_resolution_clock::now();
}

std::string class_tic_toc::get_name() const { return name; }

double class_tic_toc::get_age() const {
    return std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - start_timepoint).count();
}

double class_tic_toc::get_measured_time() const {
    if(is_measuring) {
        return std::chrono::duration_cast<std::chrono::duration<double>>(measured_time + std::chrono::high_resolution_clock::now() - tic_timepoint).count();
    } else
        return std::chrono::duration_cast<std::chrono::duration<double>>(measured_time).count();
}

double class_tic_toc::get_last_interval() const {

    if(is_measuring){
        auto delta_temp  = std::chrono::high_resolution_clock::now() - tic_timepoint;
        return std::chrono::duration_cast<std::chrono::duration<double>>(delta_temp).count();
    }else
        return std::chrono::duration_cast<std::chrono::duration<double>>(delta_time).count();
}

double class_tic_toc::get_lap() const {
    if(is_measuring) {
        // From the measured time, subtract the time since the lap time point
        auto measured_time_since_tic = std::chrono::high_resolution_clock::now() - tic_timepoint;
        auto measured_time_until_lap = lap_timepoint >= tic_timepoint ? lap_timepoint - tic_timepoint : std::chrono::high_resolution_clock::duration::zero();
        auto measured_time_since_lap =
            std::chrono::duration_cast<std::chrono::duration<double>>(lap_time + measured_time_since_tic - measured_time_until_lap).count();
        return measured_time_since_lap;
    } else
        return std::chrono::duration_cast<std::chrono::duration<double>>(lap_time).count();
}

double class_tic_toc::restart_lap() {
    double lap = get_lap();
    start_lap();
    return lap;
}


void class_tic_toc::reset() {
    if(enable) {
        measured_time   = std::chrono::high_resolution_clock::duration::zero();
        delta_time      = std::chrono::high_resolution_clock::duration::zero();
        lap_time        = std::chrono::high_resolution_clock::duration::zero();
        start_timepoint = std::chrono::high_resolution_clock::now();
        lap_timepoint   = std::chrono::high_resolution_clock::now();
        is_measuring    = false;
    }
}

class_tic_toc &class_tic_toc::operator=(double other_time_in_seconds) {
    if(enable){
        this->measured_time = std::chrono::duration_cast<hresclock::duration>(std::chrono::duration<double>(other_time_in_seconds));
        this->delta_time    = this->measured_time;
    }
    return *this;
}

class_tic_toc &class_tic_toc::operator+=(double other_time_in_seconds) {
    if(enable){
        this->measured_time += std::chrono::duration_cast<hresclock::duration>(std::chrono::duration<double>(other_time_in_seconds));
        this->delta_time = std::chrono::duration_cast<hresclock::duration>(std::chrono::duration<double>(other_time_in_seconds));
    }
    return *this;
}

class_tic_toc &class_tic_toc::operator-=(double other_time_in_seconds) {
    if(enable){
        this->measured_time -= std::chrono::duration_cast<hresclock::duration>(std::chrono::duration<double>(other_time_in_seconds));
        this->delta_time = -std::chrono::duration_cast<hresclock::duration>(std::chrono::duration<double>(other_time_in_seconds));
    }
    return *this;
}

class_tic_toc &class_tic_toc::operator+=(const class_tic_toc &rhs) {
    if(enable) {
        this->measured_time += rhs.measured_time;
        this->delta_time = rhs.measured_time;
    }
    return *this;
}

class_tic_toc &class_tic_toc::operator-=(const class_tic_toc &rhs) {
    if(enable){
        this->measured_time -= rhs.measured_time;
        this->delta_time = -rhs.measured_time;
    }
    return *this;
}



class_tic_toc::token::token(class_tic_toc &t_) :t(t_){
    t.tic();
}

class_tic_toc::token::~token() noexcept {
    try{
        if(t.is_measuring) t.toc();
    }catch(const std::exception & ex){
        fprintf(stderr,"Exception in token destructor: %s", ex.what());
    }
}

void class_tic_toc::token::tic() {t.tic();}
void class_tic_toc::token::toc() {t.toc();}