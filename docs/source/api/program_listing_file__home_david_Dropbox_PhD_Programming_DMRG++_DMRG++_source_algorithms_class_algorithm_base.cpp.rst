
.. _program_listing_file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_algorithm_base.cpp:

Program Listing for File class_algorithm_base.cpp
=================================================

|exhale_lsh| :ref:`Return to documentation for file <file__home_david_Dropbox_PhD_Programming_DMRG++_DMRG++_source_algorithms_class_algorithm_base.cpp>` (``/home/david/Dropbox/PhD/Programming/DMRG++/DMRG++/source/algorithms/class_algorithm_base.cpp``)

.. |exhale_lsh| unicode:: U+021B0 .. UPWARDS ARROW WITH TIP LEFTWARDS

.. code-block:: cpp

   //
   // Created by david on 2018-01-18.
   //
   
   #include <fstream>
   #include <complex>
   #include "class_algorithm_base.h"
   #include <io/class_hdf5_log_buffer.h>
   #include <io/nmspc_logger.h>
   #include <state/class_infinite_state.h>
   #include <state/class_environment.h>
   #include <state/class_finite_state.h>
   #include <state/tools/nmspc_tools.h>
   #include <state/class_mps_2site.h>
   #include <math/nmspc_math.h>
   #include <general/nmspc_random_numbers.h>
   #include <general/nmspc_quantum_mechanics.h>
   #include <io/log_types.h>
   #include <h5pp/h5pp.h>
   
   namespace s = settings;
   using namespace std;
   using namespace Textra;
   //using namespace std::complex_literals;
   //using namespace eigsolver_properties;
   using Scalar = class_algorithm_base::Scalar;
   
   
   
   
   class_algorithm_base::class_algorithm_base(std::shared_ptr<h5pp::File> h5ppFile_,
                                              std::string sim_name_,
                                              SimulationType sim_type_)
           : h5pp_file       (std::move(h5ppFile_)),
             sim_name       (std::move(sim_name_)),
             sim_type       (sim_type_) {
   
       log = Logger::setLogger(sim_name,settings::console::verbosity,settings::console::timestamp);
       tools::log = Logger::setLogger(sim_name,settings::console::verbosity,settings::console::timestamp);
       log->trace("Constructing class_algorithm_base");
       set_profiling_labels();
       tools::finite::profile::init_profiling();
       log_profiling  = std::make_unique<class_hdf5_log<class_log_profiling>>        (h5pp_file, sim_name + "/logs", "profiling", sim_name);
       log_sim_status = std::make_unique<class_hdf5_log<class_log_simulation_status>>(h5pp_file, sim_name + "/logs", "status"   , sim_name);
   
   
   
       log->trace("Writing input file");
       h5pp_file->writeDataset(settings::input::input_file, "common/input_file");
       h5pp_file->writeDataset(settings::input::input_filename, "common/input_filename");
   }
   
   
   
   
   
   
   bool class_algorithm_base::check_saturation_using_slope(
           std::list<bool>  & B_vec,
           std::list<double> &Y_vec,
           std::list<int> &X_vec,
           double new_data,
           int iter,
           int rate,
           double tolerance,
           double &slope)
   {
   
       int last_measurement = X_vec.empty() ? 0 : X_vec.back();
       if (iter - last_measurement < rate){return false;}
   
       // It's time to check. Insert current numbers
       B_vec.push_back(false);
       Y_vec.push_back(new_data);
       X_vec.push_back(iter);
       unsigned long min_data_points = 2;
       if (Y_vec.size() < min_data_points){return false;}
       auto check_from =  (unsigned long)(X_vec.size()*0.75); //Check the last quarter of the measurements in Y_vec.
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
   
       slope = std::abs(numerator / denominator) / avgY;
       slope = std::isnan(slope) ? 0.0 : slope;
       //Scale the slope so that it can be interpreted as change in percent, just as the tolerance.
       bool has_saturated;
       if (slope < tolerance){
           B_vec.back() = true;
           has_saturated = true;
       }else{
           B_vec.clear();
           has_saturated = false;
       }
       log->debug("Slope details:");
       log->debug(" -- relative slope  = {} %", slope);
       log->debug(" -- tolerance       = {} ", tolerance);
       log->debug(" -- avgY            = {} ", avgY);
       log->debug(" -- has saturated   = {} ", has_saturated);
       log->debug(" -- check from      = {} of {}", check_from, X_vec.size());
       return has_saturated;
   }
   
   void class_algorithm_base::write_status(bool force){
       if (not force){
           if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
           if (write_freq() == 0){return;}
           if (settings::hdf5::storage_level <= StorageLevel::NONE){return;}
       }
       log->trace("Writing simulation status to file");
       h5pp_file->writeDataset(false, sim_name + "/simOK");
       tools::common::io::write_simulation_status(sim_status, *h5pp_file, sim_name);
       h5pp_file->writeDataset(true, sim_name + "/simOK");
   }
   
   
   void class_algorithm_base::update_bond_dimension(){
       sim_status.chi_max = chi_max();
       if(not chi_grow() or sim_status.bond_dimension_has_reached_max or sim_status.chi_temp == chi_max() ){
           sim_status.chi_temp = chi_max();
           sim_status.bond_dimension_has_reached_max = true;
       }
       if(not sim_status.simulation_has_converged
          and sim_status.simulation_has_saturated
          and sim_status.chi_temp < chi_max()){
           log->trace("Updating bond dimension");
           sim_status.chi_temp = std::min(chi_max(), sim_status.chi_temp * 2);
           log->info("New chi = {}", sim_status.chi_temp);
           clear_saturation_status();
       }
       if(sim_status.chi_temp == chi_max()){
           sim_status.bond_dimension_has_reached_max = true;
       }
   }
   
   
   
   
   
   
   //void class_algorithm_base::store_profiling_deltas(bool force) {
   //    if(not force){
   //        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
   //        if (not settings::profiling::on or not settings::hdf5::save_profiling){return;}
   //        if (settings::hdf5::storage_level < StorageLevel::NORMAL){return;}
   //    }
   //
   //    log->trace("Storing profiling deltas");
   //}
   //
   //void class_algorithm_base::store_profiling_totals(bool force) {
   //    if(not force){
   //        if (math::mod(sim_status.iteration, write_freq()) != 0) {return;}
   //        if (not settings::profiling::on or not settings::hdf5::save_profiling){return;}
   //        if (settings::hdf5::storage_level < StorageLevel::NORMAL){return;}
   //    }
   //
   //    log->trace("Storing profiling totals");
   //
   //}
   
   
   
   double class_algorithm_base::process_memory_in_mb(std::string name){
       ifstream filestream("/proc/self/status");
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
       t_sim.set_properties(on,   precision,"↳+Simulation             ");
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
