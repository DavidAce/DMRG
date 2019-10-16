//
// Created by david on 2018-01-18.
//

#include <fstream>
#include <complex>
#include "class_algorithm_base.h"
#include <io/class_hdf5_log_buffer.h>
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
    log->trace("Constructing log buffers in base");
    log_profiling  = std::make_unique<class_hdf5_log<class_log_profiling>>        (h5pp_file, sim_name + "/logs", "profiling", sim_name);
    log_sim_status = std::make_unique<class_hdf5_log<class_log_simulation_status>>(h5pp_file, sim_name + "/logs", "status"   , sim_name);



    log->trace("Writing input file");
    h5pp_file->writeDataset(settings::input::input_file, "common/input_file");
    h5pp_file->writeDataset(settings::input::input_filename, "common/input_filename");
}





class_algorithm_base::SaturationReport
class_algorithm_base::check_saturation_using_slope(
        std::list<bool>  & B_vec,
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
    B_vec.push_back(false);
    Y_vec.push_back(new_data);
    X_vec.push_back(iter);
    unsigned long min_data_points = 2;
    if (Y_vec.size() < min_data_points){return report;}
//    [2019-09-04 09:41:16][xDMRG][  info  ] Entanglement Entropies  = {-0, 0.107435, 0.0755767, 0.689875, 0.692682, 0.709075, 0.936858, 0.771467, 0.609487, 0.637618, 0.708048, 0.703904, 0.716228, 0.131454 , 0.0982134, 0.165784, -0}
//    [2019-09-04 14:09:38][xDMRG][  info  ] Entanglement Entropies  = {-0, 0.107185, 0.0750735, 0.689648, 0.69319 , 0.705737, 0.764265, 0.71776 , 0.640145, 0.658814, 0.706533, 0.687854, 0.659679, 0.0969757, 0.0707726, 0.14884 , -0}
//    [2019-09-04 15:50:25][xDMRG][  info  ] Entanglement Entropies  = {-0, 0.109655, 0.0794533, 0.690183, 0.633393, 0.638498, 0.852765, 0.680258, 0.604823, 0.636992, 0.707849, 0.693525, 0.664779, 0.0980222, 0.191848 , 0.249013, -0} xDMRG Iter: 7     E: 1.1380981061517343    ε: 0.5688  log₁₀ σ²(E): -13.0853826559
//    [2019-09-04 16:12:27][xDMRG][  info  ] Entanglement Entropies  = {-0, 0.260162, 0.0767326, 0.691241, 0.693091, 0.704597, 0.744191, 0.683444, 0.603736, 0.63557 , 0.705953, 0.693286, 0.664763, 0.0973455, 0.0711479, 0.149059, -0} xDMRG Iter: 20    E: 0.3989677607237975    ε: 0.5241  log₁₀ σ²(E): -13.0536331913 (svd 1e-10)


    auto check_from =  (unsigned long)(X_vec.size()*0.75); //Check from last part of the measurements in Y_vec.
    while (X_vec.size() - check_from < min_data_points and check_from > 0){
        check_from -=1; //Decrease check from if out of bounds.
    }


    double n = X_vec.size() - check_from;
    double numerator = 0.0;
    double denominator = 0.0;


    auto x_it = X_vec.begin();
    auto y_it = Y_vec.begin();
    std::advance(x_it, check_from);
    std::advance(y_it, check_from);

    auto v_end = Y_vec.end();
    double avgX = accumulate(x_it, X_vec.end(), 0.0) / n;
    double avgY = accumulate(y_it, Y_vec.end(), 0.0) / n;

    while(y_it != v_end){
        numerator   += (*x_it - avgX) * (*y_it - avgY);
        denominator += (*x_it - avgX) * (*x_it - avgX);
        y_it++;
        x_it++;

    }

    double slope = std::abs(numerator / denominator) / avgY * 100;
    slope       = std::isnan(slope) ? 0.0 : slope;
    //Scale the slope so that it can be interpreted as change in percent, just as the tolerance.
    bool has_saturated;
    if (slope < tolerance){
        B_vec.back()  = true;
        has_saturated = true;
    }else{
        B_vec.clear();
        has_saturated = false;
    }
    report.has_computed  = true;
    report.has_saturated = has_saturated;
    report.slope         = slope;
    report.avgY          = avgY;
    report.check_from    = check_from;
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
 * We want to check once every "rate" steps. First, check the sim_state.iteration number when you last measured.
 * If the last measurement happened less than rate iterations ago, return.
 * Starting from the last measurement, and including at least 2 data points, check how far back you can go before
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
    t_run.set_properties(on, precision, "↳+Simulation             ");
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
