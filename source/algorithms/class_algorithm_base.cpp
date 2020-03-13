//
// Created by david on 2018-01-18.
//

#include <fstream>
#include <sstream>
#include <complex>
#include "class_algorithm_base.h"
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <h5pp/h5pp.h>
#include <simulation/nmspc_settings.h>
#include <math/nmspc_math.h>

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
    tools::common::profile::init_profiling();
//    if (settings::output::storage_level >= StorageLevel::NORMAL){
//        log->trace("Constructing table buffers in base");
//        h5tbuf_profiling  = std::make_unique<class_h5table_buffer<class_h5table_profiling>>        (h5pp_file, sim_name + "/journal/profiling");
//        h5tbuf_sim_status = std::make_unique<class_h5table_buffer<class_h5table_simulation_status>>(h5pp_file, sim_name + "/journal/sim_status");
//    }


    if(h5pp_file) log->trace("Writing input file");
    if(h5pp_file) h5pp_file->writeDataset(settings::input::config_filename, "common/input_filename");
    if(h5pp_file) h5pp_file->writeDataset(settings::input::config_file_contents, "common/input_file");
}




class_algorithm_base::SaturationReport
class_algorithm_base::check_saturation_using_slope(
//        std::list<bool>  & B_vec,
        std::list<double> &Y_vec,
        std::list<int> &X_vec,
        double new_data,
        int iter,
        int rate,
        double tolerance)
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
    size_t start_point = 0;
    double band_size   = 2.0 + 2.0*tolerance;  // Between 2 and  4 standard deviations away

    size_t recent_point   = std::floor(0.75*Y_vec.size());
    recent_point = std::min(Y_vec.size()-min_data_points , recent_point);
    double recent_point_std = math::stdev(Y_vec, recent_point); //Computes the standard dev of Y_vec from recent_point to end
    for(size_t some_point = 0; some_point < Y_vec.size(); some_point++){
        double some_point_std = math::stdev(Y_vec, some_point); //Computes the standard dev of Y_vec from some_point to end
        std::string arrow = "";
        if(some_point_std < band_size * recent_point_std and start_point == 0){
            start_point = some_point;
            break;
        }
    }
    //Scale the slope so that it can be interpreted as change in percent, just as the tolerance.
    double avgY          = math::mean(Y_vec,start_point);
    double slope         = math::slope(X_vec,Y_vec,start_point)/avgY * 100 / std::sqrt(Y_vec.size()-start_point); //TODO: Is dividing by sqrt(elems) reasonable?
    slope                = std::isnan(slope) ? 0.0 : slope;
    report.slope         = slope;
    report.check_from    = start_point;
    report.avgY          = avgY;
    report.has_computed  = true;
    return report;
}

//
//class_algorithm_base::SaturationReport2
//class_algorithm_base::check_saturation_using_slope2(
//        std::list<double> &Y_vec,
//        std::list<int>    &X_vec,
//        double new_data,
//        int iter,
//        int rate,
//        double tolerance)
///*! \brief Checks convergence based on slope.
// * NOTE! THIS FUNCTION REQUIRES MONOTONICALLY DECREASING Y-elements
// * We want to check once every "rate" steps. First, check the sim_state.iteration number when you last measured.
// * If the last measurement happened less than rate iterations ago, return.
// * Starting from the last measurement, and including at last 2 data points, check how far back you can go before
// * the slope to grows larger than the threshold.
// * The slope here is defined as the relative slope, i.e. \f$ \frac{1}{ \langle y\rangle} * \frac{dy}{dx} \f$.
// */
//
//{
//    SaturationReport2 report;
//    int last_measurement = X_vec.empty() ? 0 : X_vec.back();
//    if (iter - last_measurement < rate){return report;}
//
//    // It's time to check. Insert current numbers
//    Y_vec.push_back(new_data);
//    X_vec.push_back(iter);
//    unsigned long data_points = 0;
//    while(data_points <= Y_vec.size()){
//        auto x_it = X_vec.end();
//        auto y_it = Y_vec.end();
//        std::advance(x_it, -data_points);
//        std::advance(y_it, -data_points);
//        if (data_points >= 2){
//            double numerator   = 0.0;
//            double denominator = 0.0;
//            auto v_end = Y_vec.end();
//            double avgX = accumulate(x_it, X_vec.end(), 0.0) / (double)data_points;
//            double avgY = accumulate(y_it, Y_vec.end(), 0.0) / (double)data_points;
//            while(y_it != v_end){
//                numerator   += (*x_it - avgX) * (*y_it - avgY);
//                denominator += (*x_it - avgX) * (*x_it - avgX);
//                y_it++;
//                x_it++;
//            }
//
//            double slope = std::abs(numerator / denominator) / avgY * 100;
//            slope        = std::isnan(slope) ? 0.0 : slope;
//
//            report.has_computed  = true;
//            report.slopes.push_back(slope);
//            report.avgY.push_back(avgY);
//        }
//        if(x_it == X_vec.begin()) break;
//        if(y_it == Y_vec.begin()) break;
//        data_points++;
//    }
//
//    if(report.has_computed){
////        auto first_greater_than_tolerance = std::distance(report.slopes.begin(), std::upper_bound(report.slopes.begin(),report.slopes.end(),tolerance));
//        auto first_greater_than_tolerance = std::distance(report.slopes.begin(),
//                std:: find_if(report.slopes.begin(), report.slopes.end(),[tolerance](const double & x) { return x > tolerance; }));
//        report.saturated_for = first_greater_than_tolerance;
//        report.has_saturated = report.saturated_for > 0;
//        std::reverse(report.slopes.begin(),report.slopes.end()); //Reverse looks better on print
//    }
//    return report;
//}



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

