//
// Created by david on 2018-01-18.
//

#include <fstream>
#include <complex>
#include "class_algorithm_base.h"
#include <io/class_hdf5_table_buffer2.h>
#include <io/nmspc_logger.h>
#include <state/class_infinite_state.h>
#include <state/class_environment.h>
#include <state/class_finite_state.h>
#include <state/tools/nmspc_tools.h>
#include <state/class_mps_2site.h>
#include <math/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/nmspc_quantum_mechanics.h>
#include <algorithms/table_types.h>

#include <h5pp/h5pp.h>
/*! \brief Prints the content of a vector nicely */
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::list<T> &v) {
    if (!v.empty()) {
        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}


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

    log->trace("Writing input file");
    h5pp_file->writeDataset(settings::input::input_file, "common/input_file");
    h5pp_file->writeDataset(settings::input::input_filename, "common/input_filename");
}






void class_algorithm_base::check_saturation_using_slope(
        std::list<bool>  & B_vec,
        std::list<double> &Y_vec,
        std::list<int> &X_vec,
        double new_data,
        int iter,
        int rate,
        double tolerance,
        double &slope,
        bool   &has_saturated){
    //Check convergence based on slope.
    log->trace("Checking saturation using slope");


    // We want to check once every "rate" steps
    // Get the sim_status.iteration number when you last measured.
    // If the measurement happened less than rate iterations ago, return.
    int last_measurement = X_vec.empty() ? 0 : X_vec.back();
    if (iter - last_measurement < rate){return;}

    // It's time to check. Insert current numbers
    B_vec.push_back(false);
    Y_vec.push_back(new_data);
    X_vec.push_back(iter);
    unsigned long min_data_points = 2;
    if (Y_vec.size() < min_data_points){return;}
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

    slope = std::abs(numerator / denominator);
    //Scale the slope so that it can be interpreted as change in percent, just as the tolerance.
    double relative_slope     = slope / avgY;

    if (relative_slope < tolerance){
        B_vec.back() = true;
        has_saturated = true;
    }else{
        B_vec.clear();
        has_saturated = false;
    }
    log->debug("Slope details:");
    log->debug(" -- change per step = {}              | log₁₀ = {}", slope, std::log10(slope));
    log->debug(" -- relative_slope  = {} ", relative_slope);
    log->debug(" -- tolerance       = {} ", tolerance);
    log->debug(" -- avgY            = {} ", avgY);
    log->debug(" -- has saturated   = {} ", has_saturated);
    log->debug(" -- check from      = {}              | {} ", check_from, X_vec.size());
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




void class_algorithm_base::store_algorithm_state_to_file(){
    if (settings::hdf5::storage_level < StorageLevel::LIGHT){return;}
    log->trace("Storing simulation state to file");
    t_sto.tic();
    tools::common::io::write_algorithm_state(sim_status, *h5pp_file, sim_name);
    t_sto.toc();
}


void class_algorithm_base::store_profiling_deltas(bool force) {
    if(not force){
        if (math::mod(sim_status.iteration, store_freq()) != 0) {return;}
        if (not settings::profiling::on or not settings::hdf5::store_profiling){return;}
        if (settings::hdf5::storage_level < StorageLevel::NORMAL){return;}
    }

    log->trace("Storing profiling deltas");
    t_sto.tic();
    table_profiling->append_record(
            sim_status.iteration,
            t_tot.get_last_time_interval(),
            t_opt.get_last_time_interval(),
            t_sim.get_last_time_interval(),
            t_svd.get_last_time_interval(),
            t_env.get_last_time_interval(),
            t_evo.get_last_time_interval(),
            t_udt.get_last_time_interval(),
            t_sto.get_last_time_interval(),
            t_ste.get_last_time_interval(),
            t_prt.get_last_time_interval(),
            t_obs.get_last_time_interval(),
            t_mps.get_last_time_interval(),
            t_con.get_last_time_interval()
    );
    t_sto.toc();
}

void class_algorithm_base::store_profiling_totals(bool force) {
    if(not force){
        if (math::mod(sim_status.iteration, store_freq()) != 0) {return;}
        if (not settings::profiling::on or not settings::hdf5::store_profiling){return;}
        if (settings::hdf5::storage_level < StorageLevel::NORMAL){return;}
    }

    log->trace("Storing profiling totals");
    table_profiling->append_record(
            sim_status.iteration,
            t_tot.get_measured_time(),
            t_opt.get_measured_time(),
            t_sim.get_measured_time(),
            t_svd.get_measured_time(),
            t_env.get_measured_time(),
            t_evo.get_measured_time(),
            t_udt.get_measured_time(),
            t_sto.get_measured_time(),
            t_ste.get_measured_time(),
            t_prt.get_measured_time(),
            t_obs.get_measured_time(),
            t_mps.get_measured_time(),
            t_con.get_measured_time()
    );

}

//void class_algorithm_base::initialize_superblock(std::string initial_state) {
//    log->trace("Initializing state: {}", initial_state);
//    //Set the size and initial values for the MPS and environments
//    //Choose between GHZ, W, Random, Product state (up, down, etc), None, etc...
//    long d    = state->HA->get_spin_dimension();
//    long chiA = state->MPS->chiA();
//    long chiB = state->MPS->chiB();
//    Eigen::Tensor<Scalar,1> LA;
//    Eigen::Tensor<Scalar,3> GA;
//    Eigen::Tensor<Scalar,1> LC;
//    Eigen::Tensor<Scalar,3> GB;
//    Eigen::Tensor<Scalar,1> LB;
//    Eigen::Tensor<Scalar,4> theta;
//
//    GA.setZero();
//    GB.setZero();
//
//    if(initial_state == "upup"){
//        log->info("Initializing Up-Up-state  |up,up>");
//        GA.resize(array3{d,1,1});
//        GB.resize(array3{d,1,1});
//        LA.resize(array1{1});
//        LB.resize(array1{1});
//        LC.resize(array1{1});
//        LA.setConstant(1.0);
//        LB.setConstant(1.0);
//        LC.setConstant(1.0);
//        GA(0, 0, 0) = 1;
//        GB(0, 0, 0) = 1;
//        state->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = state->get_theta();
//        state->truncate_MPS(theta, 1, settings::precision::SVDThreshold);
//    }else if(initial_state == "updown"){
//        log->info("Initializing Up down -state  |up,down>");
//        GA.resize(array3{d,1,1});
//        GB.resize(array3{d,1,1});
//        LA.resize(array1{1});
//        LB.resize(array1{1});
//        LC.resize(array1{1});
//        LA.setConstant(1.0);
//        LB.setConstant(1.0);
//        LC.setConstant(1.0);
//        GA(0  , 0, 0) = 1;
//        GB(d-1, 0, 0) = 1;
//        state->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = state->get_theta();
//        state->truncate_MPS(theta, 1, settings::precision::SVDThreshold);
//    }else if(initial_state == "ghz"){
//        log->info("Initializing GHZ-statee");
//        // GHZ state (|up,up> + |down, down > ) /sqrt(2)
//        GA.resize(array3{d,1,2});
//        GB.resize(array3{d,2,1});
//        LA.resize(array1{1});
//        LB.resize(array1{1});
//        LC.resize(array1{2});
//        LA.setConstant(1.0);
//        LB.setConstant(1.0);
//        LC.setConstant(1.0/std::sqrt(2));
//
//        // GA^0 = (1,0)
//        // GA^1 = (0,1)
//        // GB^0 = (1,0)^Scalar_
//        // GB^1 = (0,1)^Scalar_
//        GA(0, 0, 0) = 1;
//        GA(0, 0, 1) = 0;
//        GA(1, 0, 0) = 0;
//        GA(1, 0, 1) = 1;
//        GB(0, 0, 0) = 1;
//        GB(0, 1, 0) = 0;
//        GB(1, 0, 0) = 0;
//        GB(1, 1, 0) = 1;
//        state->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = state->get_theta();
//        state->truncate_MPS(theta, 2, settings::precision::SVDThreshold);
//    }else if(initial_state == "lambda"){
//        log->info("Initializing W-state");
//        // W state (|up,down> + |down, up > ) /sqrt(2)
//        GA.resize(array3{d,1,2});
//        GB.resize(array3{d,2,1});
//        LA.resize(array1{1});
//        LB.resize(array1{1});
//        LC.resize(array1{2});
//        LA.setConstant(1.0);
//        LB.setConstant(1.0);
//        LC.setConstant(1.0/std::sqrt(2));
//        // GA^0 = (1,0)
//        // GA^1 = (0,1)
//        // GB^0 = (0,1)^Scalar_
//        // GB^1 = (1,0)^Scalar_
//        GA(0, 0, 0) = 1;
//        GA(0, 0, 1) = 0;
//        GA(1, 0, 0) = 0;
//        GA(1, 0, 1) = 1;
//        GB(0, 0, 0) = 0;
//        GB(0, 1, 0) = 1;
//        GB(1, 0, 0) = 1;
//        GB(1, 1, 0) = 0;
//        state->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = state->get_theta();
//        state->truncate_MPS(theta, 2, settings::precision::SVDThreshold);
//    }
//
//    else if (initial_state == "rps"){
//        // Random product state
//        log->info("Initializing random product state");
//
//        //Initialize as spinors
//        theta = Textra::Matrix_to_Tensor(Eigen::MatrixXcd::Random(d*chiA,d*chiB),d,chiA,d,chiB);
//        state->truncate_MPS(theta, 1, settings::precision::SVDThreshold);
//
//    }else if (initial_state == "random_chi" ){
//        // Random state
//        log->info("Initializing random state with bond dimension chi = {}", chi_max());
//        sim_status.chi_temp = chi_max();
//        GA.resize(array3{d,chi_max(),chi_max()});
//        GB.resize(array3{d,chi_max(),chi_max()});
//        LA.resize(array1{chi_max()});
//        LB.resize(array1{chi_max()});
//        LC.resize(array1{chi_max()});
//        LA.setConstant(1.0/sqrt(chi_max()));
//        LB.setConstant(1.0/sqrt(chi_max()));
//        LC.setConstant(1.0/sqrt(chi_max()));
//        state->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = Textra::Matrix_to_Tensor(Eigen::MatrixXcd::Random(d*chi_max(),d*chi_max()),d,chi_max(),d,chi_max());
//        state->truncate_MPS(theta, chi_max(), settings::precision::SVDThreshold);
//
//    }else{
//        std::cerr << "Invalid state given for initialization. Check 'model::initial_state' your input file. Please choose one of: " << std::endl;
//        std::cerr << "  upup" << std::endl;
//        std::cerr << "  updown" << std::endl;
//        std::cerr << "  GHZ" << std::endl;
//        std::cerr << "  W" << std::endl;
//        std::cerr << "  rps" << std::endl;
//        std::cerr << "  random_chi (only for iDMRG!)" << std::endl;
//        exit(1);
//    }
//
//
//
//
//    //Reset the environment blocks to the correct dimensions
//    state->Lblock->set_edge_dims(*state->MPS, state->HA->MPO());
//    state->Rblock->set_edge_dims(*state->MPS, state->HB->MPO());
//    state->Lblock2->set_edge_dims(*state->MPS, state->HA->MPO());
//    state->Rblock2->set_edge_dims(*state->MPS, state->HB->MPO());
//
//    state->environment_size = state->Lblock->size + state->Rblock->size;
//
//    assert(state->Lblock->block.dimension(0) == state->MPS->chiA());
//    assert(state->Rblock->block.dimension(0) == state->MPS->chiB());
//
//
//
//    if(sim_type == SimulationType::fDMRG or sim_type == SimulationType::xDMRG ){
//        tools::finite::state::insert_superblock_to_state(*state, *state);
//    }else{
//    }
//
//    enlarge_environment();
//
//    if (sim_type == SimulationType::iDMRG){
//        sim_status.iteration = (int)state->Lblock->size;
//    }
//    swap();
//}



//
//
//
//void class_algorithm_base::insert_superblock_to_chain() {
//    log->trace("Insert state into state");
//    t_sim.tic();
//    t_ste.tic();
//    auto new_position = tools::finite::state::insert_superblock_to_state(*state, *state);
//    state->set_positions(new_position);
//    t_ste.toc();
//    t_sim.toc();
//}
//
//void class_algorithm_base::copy_superblock_mps_to_chain(){
//    log->trace("Copy state mps to state");
//    t_sim.tic();
//    t_ste.tic();
//    tools::finite::state::copy_superblock_mps_to_state(*state, *state);
//    t_ste.toc();
//    t_sim.toc();
//}
//
//void class_algorithm_base::copy_superblock_mpo_to_chain(){
//    log->trace("Copy state mpo to state");
//    t_sim.tic();
//    t_ste.tic();
//    tools::finite::state::copy_superblock_mpo_to_state(*state, *state);
//    t_ste.toc();
//    t_sim.toc();
//}
//
//void class_algorithm_base::copy_superblock_env_to_chain(){
//    log->trace("Copy state env to state");
//    t_sim.tic();
//    t_ste.tic();
//    tools::finite::state::copy_superblock_env_to_state(*state, *state);
//    t_ste.toc();
//    t_sim.toc();
//}
//
//void class_algorithm_base::copy_superblock_to_chain(){
//    log->trace("Copy state to state");
//    t_sim.tic();
//    t_ste.tic();
//    tools::finite::state::copy_superblock_to_state(*state, *state);
//    t_ste.toc();
//    t_sim.toc();
//}
//



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
    t_sto.set_properties(on,   precision,"↳ Store to file          ");
    t_ste.set_properties(on,   precision,"↳ finite state storage   ");
    t_prt.set_properties(on,   precision,"↳ Printing to console    ");
    t_obs.set_properties(on,   precision,"↳ Computing observables  ");
    t_sim.set_properties(on,   precision,"↳+Simulation             ");
    t_evo.set_properties(on,   precision,"↳ Time Evolution         ");
    t_opt.set_properties(on,   precision,"↳+Optimize MPS           ");
    t_eig.set_properties(on,   precision," ↳ Eigenvalue solver     ");
    t_ham.set_properties(on,   precision," ↳ Build Hamiltonian     ");
    t_svd.set_properties(on,   precision,"↳ SVD Truncation         ");
    t_udt.set_properties(on,   precision,"↳ Update Timestep        ");
    t_env.set_properties(on,   precision,"↳ Update Environments    ");
    t_con.set_properties(on,   precision,"↳ Check Convergence      ");
}
