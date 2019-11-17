//
// Created by david on 2018-01-18.
//

#include <fstream>
#include <complex>
#include "class_algorithm_base.h"
#include <io/class_h5table_buffer.h>
#include <io/nmspc_logger.h>
#include <tools/nmspc_tools.h>
#include <h5pp/h5pp.h>
#include <simulation/nmspc_settings.h>

using Scalar = class_algorithm_base::Scalar;




class_algorithm_base::class_algorithm_base(std::shared_ptr<h5pp::File> h5ppFile_,
                                           std::string sim_name_,
                                           SimulationType sim_type_)
        : h5pp_file       (std::move(h5ppFile_)),
          sim_name       (std::move(sim_name_)),
          sim_type       (sim_type_) {

    log        = Logger::setLogger(sim_name,settings::console::verbosity,settings::console::timestamp);
    tools::log = Logger::setLogger(sim_name,settings::console::verbosity,settings::console::timestamp);
    log->trace("Constructing class_algorithm_base");
    set_profiling_labels();
    tools::common::profile::init_profiling();
    if (settings::output::storage_level >= StorageLevel::NORMAL){
        log->trace("Constructing table buffers in base");
        h5tbuf_profiling  = std::make_unique<class_h5table_buffer<class_h5table_profiling>>        (h5pp_file, sim_name + "/journal/profiling");
        h5tbuf_sim_status = std::make_unique<class_h5table_buffer<class_h5table_simulation_status>>(h5pp_file, sim_name + "/journal/sim_status");
    }


    if(h5pp_file) log->trace("Writing input file");
    if(h5pp_file) h5pp_file->writeDataset(settings::input::input_filename  , "common/input_filename");
    if(h5pp_file) h5pp_file->writeDataset(settings::input::input_file_raw  , "common/input_file");
}

std::list<double> get_cumsum(std::list<int> &X_vec,std::list<double> &Y_vec){
    if (X_vec.size() != Y_vec.size())throw std::logic_error("Lists must be equal in size!");
    auto Y_abs = Y_vec;
    std::for_each(Y_abs.begin(),Y_abs.end(), std::abs<double>);
    double minval = *std::min_element(Y_abs.begin(),Y_abs.end());
    std::for_each(Y_abs.begin(),Y_abs.end(),[minval](auto &val){val-=minval;});
    double maxval = *std::max_element(Y_abs.begin(),Y_abs.end());
    std::for_each(Y_abs.begin(),Y_abs.end(),[maxval](auto &val){val/=maxval;});
    //Now data should be normalized between 0 and 1.

    double sum = 0.0;
    std::list<double> cumsum;
    std::list<double> deltaX(X_vec.size());
    std::adjacent_difference(X_vec.begin(),X_vec.end(),deltaX.begin());
    deltaX.front() = 0;
    auto d_it = deltaX.begin();
    auto y_it = Y_abs.begin();
    while(d_it != deltaX.end()){
        sum += *d_it * *y_it;
        cumsum.emplace_back(sum);
        d_it++;
        y_it++;
    }
    return cumsum;
}

std::list<std::list<double>> get_slopes(std::list<int> &X_vec,std::list<double> &Y_vec){
    if (X_vec.size() != Y_vec.size())throw std::logic_error("Lists must be equal in size!");
    size_t size = X_vec.size();

    std::list<std::list<double>> slope_lists;
    for(size_t win_size = size; win_size > 2; win_size --){
        auto x_it = X_vec.begin();
        auto y_it = Y_vec.begin();
        std::list<double> slopes;
        while(std::distance(x_it,X_vec.end()) >= (int) win_size){
            std::vector<int>    x_win(x_it, std::next(x_it,win_size) );
            std::vector<double> y_win(y_it, std::next(y_it,win_size) );
            double avgX = accumulate(x_win.begin(), x_win.end(), 0.0) / x_win.size();
            double avgY = accumulate(y_win.begin(), y_win.end(), 0.0) / y_win.size();
            double numerator   = 0.0;
            double denominator = 0.0;
            for(size_t i = 0; i < win_size; i++){
                numerator   += (x_win[i] - avgX) * (y_win[i] - avgY);
                denominator += (x_win[i] - avgX) * (x_win[i] - avgX);
            }
            double slope = std::abs(numerator / denominator) / avgY * 100;
            slope        = std::isnan(slope) ? 0.0 : slope;
            slopes.push_back(slope);
            std::advance(x_it, 1);
            std::advance(y_it, 1);
        }
        slope_lists.push_back(slopes);
    }
    return slope_lists;

}

class_algorithm_base::SaturationReport
class_algorithm_base::check_saturation_using_slope(
//        std::list<bool>  & B_vec,
        std::list<double> &Y_vec,
        std::list<int> &X_vec,
        double new_data,
        int iter,
        int rate)
/*! \brief Checks convergence based on slope.
 * We want to check once every "rate" steps. First, check the sim_state.iteration number when you last measured.
 * If the measurement happened less than rate iterations ago, return.
 * Otherwise, compute the slope of the last 25% of the measurements that have been made.
 * The slope here is defined as the relative slope, i.e. \f$ \frac{1}{ \langle y\rangle} * \frac{dy}{dx} \f$.
 */

{
    SaturationReport report;
    int last_measurement = X_vec.empty() ? 0 : X_vec.back();
    if (iter - last_measurement < rate){return report;}

    // It's time to check. Insert current numbers
//    B_vec.push_back(false);
    Y_vec.push_back(new_data);
    X_vec.push_back(iter);
    size_t min_data_points = 2;
    if (Y_vec.size() < min_data_points){return report;}
    size_t max_data_points = std::max(min_data_points,size_t(0.5*Y_vec.size()) ) ;

    size_t check_from = X_vec.size() - max_data_points;
    double n = X_vec.size() - check_from;
    double numerator   = 0.0;
    double denominator = 0.0;

    auto x_it = X_vec.begin();
    auto y_it = Y_vec.begin();
    std::advance(x_it, check_from);
    std::advance(y_it, check_from);
    double avgX = accumulate(x_it, X_vec.end(), 0.0) / n;
    double avgY = accumulate(y_it, Y_vec.end(), 0.0) / n;
    auto v_end = Y_vec.end();
    while(y_it != v_end){
        numerator   += (*x_it - avgX) * (*y_it - avgY);
        denominator += (*x_it - avgX) * (*x_it - avgX);
        y_it++;
        x_it++;

    }
    //Scale the slope so that it can be interpreted as change in percent, just as the tolerance.
    double slope = std::abs(numerator / denominator) / avgY * 100;
    slope       = std::isnan(slope) ? 0.0 : slope;
    report.slope = slope;
    report.avgY  = avgY;
    report.has_computed  = true;
    return report;
}


class_algorithm_base::SaturationReport2
class_algorithm_base::check_saturation_using_slope2(
        std::list<double> &Y_vec,
        std::list<int>    &X_vec,
        double new_data,
        int iter,
        int rate,
        double tolerance)
/*! \brief Checks convergence based on slope.
 * NOTE! THIS FUNCTION REQUIRES MONOTONICALLY DECREASING Y-elements
 * We want to check once every "rate" steps. First, check the sim_state.iteration number when you last measured.
 * If the last measurement happened less than rate iterations ago, return.
 * Starting from the last measurement, and including at last 2 data points, check how far back you can go before
 * the slope to grows larger than the threshold.
 * The slope here is defined as the relative slope, i.e. \f$ \frac{1}{ \langle y\rangle} * \frac{dy}{dx} \f$.
 */

{
    SaturationReport2 report;
    int last_measurement = X_vec.empty() ? 0 : X_vec.back();
    if (iter - last_measurement < rate){return report;}

    // It's time to check. Insert current numbers
    Y_vec.push_back(new_data);
    X_vec.push_back(iter);
    unsigned long data_points = 0;
    while(data_points <= Y_vec.size()){
        auto x_it = X_vec.end();
        auto y_it = Y_vec.end();
        std::advance(x_it, -data_points);
        std::advance(y_it, -data_points);
        if (data_points >= 2){
            double numerator   = 0.0;
            double denominator = 0.0;
            auto v_end = Y_vec.end();
            double avgX = accumulate(x_it, X_vec.end(), 0.0) / (double)data_points;
            double avgY = accumulate(y_it, Y_vec.end(), 0.0) / (double)data_points;
            while(y_it != v_end){
                numerator   += (*x_it - avgX) * (*y_it - avgY);
                denominator += (*x_it - avgX) * (*x_it - avgX);
                y_it++;
                x_it++;
            }

            double slope = std::abs(numerator / denominator) / avgY * 100;
            slope        = std::isnan(slope) ? 0.0 : slope;

            report.has_computed  = true;
            report.slopes.push_back(slope);
            report.avgY.push_back(avgY);
        }
        if(x_it == X_vec.begin()) break;
        if(y_it == Y_vec.begin()) break;
        data_points++;
    }

    if(report.has_computed){
//        auto first_greater_than_tolerance = std::distance(report.slopes.begin(), std::upper_bound(report.slopes.begin(),report.slopes.end(),tolerance));
        auto first_greater_than_tolerance = std::distance(report.slopes.begin(),
                std:: find_if(report.slopes.begin(), report.slopes.end(),[tolerance](const double & x) { return x > tolerance; }));
        report.saturated_for = first_greater_than_tolerance;
        report.has_saturated = report.saturated_for > 0;
        std::reverse(report.slopes.begin(),report.slopes.end()); //Reverse looks better on print
    }
    return report;
}


void class_algorithm_base::print_profiling(){
    if (settings::profiling::on) {
        log->trace("Printing profiling information (tot)");
        t_tot.print_time_w_percent();
        t_run.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_con.print_time_w_percent(t_tot);
        tools::common::profile::print_profiling(t_tot);
    }
}



double class_algorithm_base::process_memory_in_mb(std::string name){
    std::ifstream filestream("/proc/self/status");
    std::string line;
    while (std::getline(filestream, line)){
        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, ':')){
            if (key == name){
                std::string value_str;
                if (std::getline(is_line, value_str)) {
                    // Extract the number
                    std::string::size_type sz;   // alias of size_t
                    int value = std::stoi (value_str,&sz);
                    // Now we have the value in kb
                    return value/1024.0;
//                    auto pos = value.find_first_not_of(" \t");
//                    auto trimmed_value = value.substr(pos != std::string::npos ? pos : 0);
//                    return trimmed_value;
                }
            }
        }
    }

    return -1.0;
}

void class_algorithm_base::set_profiling_labels() {
    using namespace settings::profiling;
    t_tot.set_properties(true, precision,"+Total Time              ");
    t_prt.set_properties(on,   precision,"↳ Printing to console    ");
    t_con.set_properties(on,   precision,"↳ Convergence checks     ");
    t_run.set_properties(on,   precision, "↳+Simulation             ");
//    t_obs.set_properties(on,   precision,"↳ Computing observables  ");

//    t_sto.set_properties(on,   precision,"↳ Store to file          ");
//    t_ste.set_properties(on,   precision,"↳ finite state storage   ");

//    t_evo.set_properties(on,   precision,"↳ Time Evolution         ");
//    t_opt.set_properties(on,   precision,"↳+Optimize MPS           ");
//    t_eig.set_properties(on,   precision," ↳ Eigenvalue solver     ");
//    t_ham.set_properties(on,   precision," ↳ Build Hamiltonian     ");
//    t_svd.set_properties(on,   precision,"↳ SVD Truncation         ");
//    t_udt.set_properties(on,   precision,"↳ Update Timestep        ");
//    t_env.set_properties(on,   precision,"↳ Update Environments    ");
}
