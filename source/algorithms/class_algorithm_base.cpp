//
// Created by david on 2018-01-18.
//

#include <fstream>
#include <complex>
#include "class_algorithm_base.h"
#include <IO/class_hdf5_file.h>
#include <IO/class_hdf5_table_buffer2.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_finite_chain_sweeper.h>
#include <mps_routines/class_mps_2site.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/class_svd_wrapper.h>
#include <algorithms/table_types.h>


namespace s = settings;
using namespace std;
using namespace Textra;
using namespace std::complex_literals;
using namespace eigsolver_properties;
using Scalar = class_algorithm_base::Scalar;

class_algorithm_base::class_algorithm_base(std::shared_ptr<class_hdf5_file> hdf5_,
                                           std::string sim_name_,
                                           SimulationType sim_type_)
        :hdf5           (std::move(hdf5_)),
         sim_name       (std::move(sim_name_)),
         sim_type       (sim_type_) {
    set_profiling_labels();
    table_profiling = std::make_unique<class_hdf5_table<class_table_profiling>>(hdf5, sim_name, "profiling");
    superblock   = std::make_shared<class_superblock>();

    //Default constructed objects
    env_storage  = std::make_shared<class_finite_chain_sweeper>();
}


void class_algorithm_base::single_DMRG_step(long chi_max, Ritz ritz){
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();
    t_opt.tic();
    Eigen::Tensor<Scalar,4> theta = superblock->MPS->get_theta();
//    superblock->MPS->theta = superblock->MPS->get_theta();
    theta = superblock->optimize_MPS(theta, ritz);
    t_opt.toc();
    t_svd.tic();
    superblock->truncate_MPS(theta, chi_max, s::precision::SVDThreshold);
    t_svd.toc();
    //Reduce the hamiltonians if you are doing infinite systems:
    if(sim_type == SimulationType::iDMRG){
        superblock->E_optimal /= 2.0;
        superblock->HA->set_reduced_energy(superblock->E_optimal);
        superblock->HB->set_reduced_energy(superblock->E_optimal);
    }
    measurement->set_not_measured();
    t_sim.toc();
}




void class_algorithm_base::check_convergence_all(){
    t_con.tic();
    check_convergence_entanglement();
    check_convergence_variance_mpo();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    update_bond_dimension();
    if(entanglement_has_converged and
       variance_mpo_has_converged and
       variance_ham_has_converged and
       variance_mom_has_converged and
       bond_dimension_has_reached_max)
    {
        simulation_has_converged = true;
    }
    t_con.toc();
}

void class_algorithm_base::check_saturation_using_slope(
        std::list<double> &Y_vec,
        std::list<int> &X_vec,
        double new_data,
        int rate,
        double tolerance,
        double &slope,
        bool &has_saturated){
    //Check convergence based on slope.


    // We want to check once every "rate" steps
    // Get the iteration number when you last measured.
    // If the measurement happened less than rate iterations ago, return.
    int last_measurement = X_vec.empty() ? 0 : X_vec.back();
    if (iteration - last_measurement < rate){return;}

    // It's time to check. Insert current numbers and remove old ones
    Y_vec.push_back(new_data);
    X_vec.push_back(iteration);
    unsigned long min_data_points = 2;
    if (Y_vec.size() < min_data_points){return;}
    auto check_from =  (unsigned long)(X_vec.size()*0.25);
    while (X_vec.size() - check_from < min_data_points and check_from > 0){
        check_from -=1;
    }


    double n = X_vec.size() - check_from;
    double numerator = 0.0;
    double denominator = 0.0;

    double diffsq = 0.0;
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
        diffsq      += (*y_it - avgY) * (*y_it - avgY);
        y_it++;
        x_it++;

    }
    double var = diffsq / n;
    double standard_deviation = std::sqrt(var);
    slope = std::abs(numerator / denominator);
    bool reason1 = false;
    bool reason2 = false;
//    bool reason3 = false;
    //Scale the tolerance so that it is relative to the size of the values in question.
    double relative_tolerance = tolerance * avgY;

    if (slope < relative_tolerance ) {
        //If the average slope is very close to zero
        has_saturated = true;
        reason1 = true;
    }

    if(slope  <  10*relative_tolerance and standard_deviation > 0.25 * avgY and check_from > 0){
        // If the slope is not too large but still larger than the tolerance.
        // AND the "noise" is large.
        // This happens when chi is too small. Then we say that the quantity has saturated
        // because "it can't get better with given chi".
        has_saturated = true;
        reason2 = true;

    }
//    if(slope < 0.01*avgY   and standard_deviation > 0.25 * avgY ){
//        // If the average change per iteration is about 1%,
//        // AND the "noise" is large.
//        // This happens when chi is too small. Then we say that the quantity has saturated
//        // because "it can't get better with given chi".
//        has_saturated = true;
//        reason3 = true;
//    }
    if (settings::console::verbosity >= 2 and has_saturated) {
        std::cout << setprecision(16);
        std::cout << "change per step      : " << slope                 << std::endl
                  << "relative_tolerance   : " << relative_tolerance    << std::endl
                  << "standard_deviation   : " << standard_deviation    << std::endl
                  << "avgY                 : " << avgY                  << std::endl
                  << "check from           : "<< check_from             << std::endl
                  << "because"                                          << std::endl
                  << "has_saturated        : " << has_saturated         << std::endl
                  << "reason 1             : " << std::boolalpha << reason1 << std::endl
                  << "reason 2             : " << std::boolalpha << reason2 << std::endl
//                  << " \nreason 3             :" << std::boolalpha << reason3
                  << std::endl;
    }
}

void class_algorithm_base::check_convergence_variance_mpo(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(V_mpo_vec,
                                 X_mpo_vec,
                                 measurement->get_variance_mpo(),
                                 1,
                                 slope_threshold,
                                 V_mpo_slope,
                                 variance_mpo_has_saturated);
    variance_mpo_has_converged = measurement->get_variance_mpo() < threshold;
}

void class_algorithm_base::check_convergence_variance_ham(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.

    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(V_ham_vec,
                                 X_ham_vec,
                                 measurement->get_variance_ham(),
                                 1,
                                 slope_threshold,
                                 V_ham_slope,
                                 variance_ham_has_saturated);
    variance_ham_has_converged = measurement->get_variance_ham() < threshold;
}

void class_algorithm_base::check_convergence_variance_mom(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.

    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(V_mom_vec,
                                 X_mom_vec,
                                 measurement->get_variance_mom(),
                                 1,
                                 slope_threshold,
                                 V_mom_slope,
                                 variance_mom_has_saturated);
    variance_mom_has_converged = measurement->get_variance_mom() < threshold;
}

void class_algorithm_base::check_convergence_entanglement(double slope_threshold) {
    //Based on the the slope of entanglement entanglement_entropy
    // This one is cheap to compute.

    slope_threshold = std::isnan(slope_threshold) ? settings::precision::EntEntrSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(S_vec,
                                 XS_vec,
                                 measurement->get_entanglement_entropy(),
                                 1,
                                 slope_threshold,
                                 S_slope,
                                 entanglement_has_saturated);
    entanglement_has_converged = entanglement_has_saturated;
}

void class_algorithm_base::update_bond_dimension(){
    if(not chi_grow or bond_dimension_has_reached_max or chi_temp == chi_max ){
        chi_temp = chi_max;
        bond_dimension_has_reached_max = true;
    }
    if(not variance_mpo_has_converged
       and variance_mpo_has_saturated
       and chi_temp < chi_max
       and measurement->get_chi() == chi_temp){
        chi_temp = std::min(chi_max, (long)(chi_temp * 2));
        std::cout << "New chi = " << chi_temp << std::endl;
        clear_saturation_status();
    }
    if(chi_temp == chi_max){
        bond_dimension_has_reached_max = true;
    }
}

void class_algorithm_base::clear_saturation_status(){
    S_vec.clear();
    XS_vec.clear();

    V_mpo_vec.clear();
    X_mpo_vec.clear();
    V_ham_vec.clear();
    X_ham_vec.clear();
    V_mom_vec.clear();
    X_mom_vec.clear();

    entanglement_has_saturated      = false;
    variance_mpo_has_saturated      = false;
    variance_ham_has_saturated      = false;
    variance_mom_has_saturated      = false;


}

void class_algorithm_base::store_profiling_to_file_delta(bool force) {
//    if (Math::mod(iteration, store_freq) != 0) {return;}
//    t_sto.tic();

    if (force or (settings::profiling::on and settings::hdf5::store_profiling))
    table_profiling->append_record(
            iteration,
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
}

void class_algorithm_base::store_profiling_to_file_total(bool force) {
//    if (Math::mod(iteration, store_freq) != 0) {return;}
//    t_sto.tic();

    if (force or (settings::profiling::on and settings::hdf5::store_profiling))
        table_profiling->append_record(
                iteration,
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



void class_algorithm_base::initialize_state(std::string initial_state ) {
    //Set the size and initial values for the MPS and environments
    //Choose between GHZ, W, Random, Product state (up, down, etc), None, etc...
    long d    = superblock->HA->get_spin_dimension();
    long chiA = superblock->MPS->chiA();
    long chiB = superblock->MPS->chiB();
    std::srand((unsigned int) 1);
    Eigen::Tensor<Scalar,1> LA;
    Eigen::Tensor<Scalar,3> GA;
    Eigen::Tensor<Scalar,1> LC;
    Eigen::Tensor<Scalar,3> GB;
    Eigen::Tensor<Scalar,1> LB;
    Eigen::Tensor<Scalar,4> theta;

    GA.setZero();
    GB.setZero();

    if(initial_state == "upup"){
        std::cout << "Initializing Up-Up-state  |up,up>" << std::endl;
        GA.resize(array3{d,1,1});
        GB.resize(array3{d,1,1});
        LA.resize(array1{1});
        LB.resize(array1{1});
        LC.resize(array1{1});
        LA.setConstant(1.0);
        LB.setConstant(1.0);
        LC.setConstant(1.0);
        GA(0, 0, 0) = 1;
        GB(0, 0, 0) = 1;
        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
        theta = superblock->MPS->get_theta();
        superblock->truncate_MPS(theta, 1, settings::precision::SVDThreshold);
    }else if(initial_state == "updown"){
        std::cout << "Initializing Up down -state  |up,down>" << std::endl;
        GA.resize(array3{d,1,1});
        GB.resize(array3{d,1,1});
        LA.resize(array1{1});
        LB.resize(array1{1});
        LC.resize(array1{1});
        LA.setConstant(1.0);
        LB.setConstant(1.0);
        LC.setConstant(1.0);
        GA(0  , 0, 0) = 1;
        GB(d-1, 0, 0) = 1;
        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
        theta = superblock->MPS->get_theta();
        superblock->truncate_MPS(theta, 1, settings::precision::SVDThreshold);
    }else if(initial_state == "ghz"){
        std::cout << "Initializing GHZ-state" << std::endl;
        // GHZ state (|up,up> + |down, down > ) /sqrt(2)
        GA.resize(array3{d,1,2});
        GB.resize(array3{d,2,1});
        LA.resize(array1{1});
        LB.resize(array1{1});
        LC.resize(array1{2});
        LA.setConstant(1.0);
        LB.setConstant(1.0);
        LC.setConstant(1.0/std::sqrt(2));

        // GA^0 = (1,0)
        // GA^1 = (0,1)
        // GB^0 = (1,0)^T
        // GB^1 = (0,1)^T
        GA(0, 0, 0) = 1;
        GA(0, 0, 1) = 0;
        GA(1, 0, 0) = 0;
        GA(1, 0, 1) = 1;
        GB(0, 0, 0) = 1;
        GB(0, 1, 0) = 0;
        GB(1, 0, 0) = 0;
        GB(1, 1, 0) = 1;
        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
        theta = superblock->MPS->get_theta();
        superblock->truncate_MPS(theta, 2, settings::precision::SVDThreshold);
    }else if(initial_state == "lambda"){
        std::cout << "Initializing W-state" << std::endl;
        // W state (|up,down> + |down, up > ) /sqrt(2)
        GA.resize(array3{d,1,2});
        GB.resize(array3{d,2,1});
        LA.resize(array1{1});
        LB.resize(array1{1});
        LC.resize(array1{2});
        LA.setConstant(1.0);
        LB.setConstant(1.0);
        LC.setConstant(1.0/std::sqrt(2));
        // GA^0 = (1,0)
        // GA^1 = (0,1)
        // GB^0 = (0,1)^T
        // GB^1 = (1,0)^T
        GA(0, 0, 0) = 1;
        GA(0, 0, 1) = 0;
        GA(1, 0, 0) = 0;
        GA(1, 0, 1) = 1;
        GB(0, 0, 0) = 0;
        GB(0, 1, 0) = 1;
        GB(1, 0, 0) = 1;
        GB(1, 1, 0) = 0;
        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
        theta = superblock->MPS->get_theta();
        superblock->truncate_MPS(theta, 2, settings::precision::SVDThreshold);
    }

    else if (initial_state == "rps"){
        // Random product state
        std::cout << "Initializing random product state" << std::endl;
        //Initialize as spinors
        theta = Textra::Matrix_to_Tensor(Eigen::MatrixXcd::Random(d*chiA,d*chiB),d,chiA,d,chiB);
        superblock->truncate_MPS(theta, 1, settings::precision::SVDThreshold);

    }else if (initial_state == "random_chi" ){
        // Random state
        std::cout << "Initializing random state with bond dimension chi = " << chi_max << std::endl;
        chi_temp = chi_max;
        GA.resize(array3{d,chi_max,chi_max});
        GB.resize(array3{d,chi_max,chi_max});
        LA.resize(array1{chi_max});
        LB.resize(array1{chi_max});
        LC.resize(array1{chi_max});
        LA.setConstant(1.0/sqrt(chi_max));
        LB.setConstant(1.0/sqrt(chi_max));
        LC.setConstant(1.0/sqrt(chi_max));
        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
        theta = Textra::Matrix_to_Tensor(Eigen::MatrixXcd::Random(d*chi_max,d*chi_max),d,chi_max,d,chi_max);
        superblock->truncate_MPS(theta, chi_max, settings::precision::SVDThreshold);

    }else{
        std::cerr << "Invalid state given for initialization. Check 'model::initial_state' your input file. Please choose one of: " << std::endl;
        std::cerr << "  upup" << std::endl;
        std::cerr << "  updown" << std::endl;
        std::cerr << "  GHZ" << std::endl;
        std::cerr << "  W" << std::endl;
        std::cerr << "  rps" << std::endl;
        std::cerr << "  random_chi (only for iDMRG!)" << std::endl;
        exit(1);
    }




    //Reset the environment blocks to the correct dimensions
    superblock->Lblock->set_edge_dims(superblock->MPS, superblock->HA->MPO);
    superblock->Rblock->set_edge_dims(superblock->MPS, superblock->HB->MPO);
    superblock->Lblock2->set_edge_dims(superblock->MPS, superblock->HA->MPO);
    superblock->Rblock2->set_edge_dims(superblock->MPS, superblock->HB->MPO);

    superblock->environment_size = superblock->Lblock->size + superblock->Rblock->size;

    assert(superblock->Lblock->block.dimension(0) == superblock->MPS->chiA());
    assert(superblock->Rblock->block.dimension(0) == superblock->MPS->chiB());



    if(sim_type == SimulationType::fDMRG or sim_type == SimulationType::xDMRG ){
        env_storage_insert();
    }else{
    }

    enlarge_environment();

    if (sim_type == SimulationType::iDMRG){
        iteration = (int)superblock->Lblock->size;
    }
    swap();
}


void class_algorithm_base::compute_observables(){
    t_obs.tic();
    measurement->compute_all_observables_from_superblock();
    t_obs.toc();
}


void class_algorithm_base::enlarge_environment(){
    t_env.tic();
    superblock->enlarge_environment();
    t_env.toc();
}

void class_algorithm_base::enlarge_environment(int direction){
    t_env.tic();
    superblock->enlarge_environment(direction);
    t_env.toc();
}

void class_algorithm_base::swap(){
    superblock->swap_AB();
}

void class_algorithm_base::env_storage_insert() {
    t_ste.tic();
    env_storage->insert();
    t_ste.toc();
}

void class_algorithm_base::env_storage_overwrite_local_MPS(){
    t_ste.tic();
    env_storage->overwrite_local_MPS();
    t_ste.toc();
}

void class_algorithm_base::env_storage_overwrite_local_MPO(){
    t_ste.tic();
    env_storage->overwrite_local_MPO();
    t_ste.toc();
}

void class_algorithm_base::env_storage_overwrite_local_ENV(){
    t_ste.tic();
    env_storage->overwrite_local_ENV();
    t_ste.toc();
}

void class_algorithm_base::env_storage_overwrite_local_ALL(){
    t_ste.tic();
    env_storage->overwrite_local_MPS();
    env_storage->overwrite_local_MPO();
    env_storage->overwrite_local_ENV();
    t_ste.toc();
}


void class_algorithm_base::env_storage_move(){
    t_ste.tic();
    env_storage->move();
    t_ste.toc();
}


double process_memory_in_mb(std::string name){
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
void class_algorithm_base::print_status_update() {
    if (Math::mod(iteration, print_freq) != 0) {return;}
    if (not env_storage->position_is_the_middle()) {return;}
    if (print_freq == 0) {return;}
    compute_observables();
    t_prt.tic();
    std::cout << setprecision(16) << fixed << left;
    ccout(1) << left  << sim_name << " ";
    ccout(1) << left  << "Step: "                       << setw(10) << iteration;
    ccout(1) << left  << "E: ";
    switch(sim_type) {
        case SimulationType::iDMRG:
            ccout(1) << setw(21) << setprecision(16)    << fixed   << measurement->get_energy_mpo();
            ccout(1) << setw(21) << setprecision(16)    << fixed   << measurement->get_energy_ham();
            ccout(1) << setw(21) << setprecision(16)    << fixed   << measurement->get_energy_mom();
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(1) << setw(21) << setprecision(16)    << fixed   << measurement->get_energy_mpo();
            break;
        case SimulationType::iTEBD:
            ccout(1) << setw(21) << setprecision(16)    << fixed   << measurement->get_energy_ham();
            ccout(1) << setw(21) << setprecision(16)    << fixed   << measurement->get_energy_mom();
            break;
    }

    ccout(1) << left  << "log₁₀ σ²(E): ";
    switch(sim_type) {
        case SimulationType::iDMRG:
            ccout(1) << setw(12) << setprecision(4)    << fixed   << std::log10(measurement->get_variance_mpo());
            ccout(1) << setw(12) << setprecision(4)    << fixed   << std::log10(measurement->get_variance_ham());
            ccout(1) << setw(12) << setprecision(4)    << fixed   << std::log10(measurement->get_variance_mom());
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(1) << setw(12) << setprecision(4)    << fixed   << std::log10(measurement->get_variance_mpo());
            break;
        case SimulationType::iTEBD:
            ccout(1) << setw(12) << setprecision(4)    << fixed   << std::log10(measurement->get_variance_ham());
            ccout(1) << setw(12) << setprecision(4)    << fixed   << std::log10(measurement->get_variance_mom());
            break;
    }


    ccout(1) << left  << "S: "                          << setw(21) << setprecision(16)    << fixed   << measurement->get_entanglement_entropy();
    ccout(1) << left  << "χmax: "                       << setw(4)  << setprecision(3)     << fixed   << chi_max;
    ccout(1) << left  << "χ: "                          << setw(4)  << setprecision(3)     << fixed   << measurement->get_chi();
    ccout(1) << left  << "log₁₀ truncation: "           << setw(10) << setprecision(4)     << fixed   << log10(measurement->get_truncation_error());
    ccout(1) << left  << "Chain length: "               << setw(6)  << setprecision(1)     << fixed   << measurement->get_chain_length();
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(1) << left  << "@ site: "                    << setw(5)  << env_storage->get_position();
            //ccout(1) << left  << "Dir: "                    << setw(3)  << env_storage->get_direction();
            //ccout(1) << left  << "Sweep: "                  << setw(4)  << env_storage->get_sweeps();
            break;
        case SimulationType::iTEBD:
//            ccout(1) << left  << "δt: "               << setw(13) << setprecision(12)    << fixed   << delta_t;
            break;
        default:
            break;
    }
    ccout(1) << left  << " Convergence [";
    ccout(1) << left  << " S-"  << std::boolalpha << setw(6) << entanglement_has_converged;
    switch(sim_type){
        case SimulationType::iDMRG:
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(1) << left  << " σ²-"  << std::boolalpha << setw(6) << variance_mpo_has_converged;
            break;
        case SimulationType::iTEBD:
            break;
    }
    ccout(1) << left  << "]";
    ccout(1) << left  << " Saturation [";
    ccout(1) << left  << " S-"  << std::boolalpha << setw(6) << entanglement_has_saturated;
    switch(sim_type){
        case SimulationType::iDMRG:
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(1) << left  << " σ²-"  << std::boolalpha << setw(6) << variance_mpo_has_saturated;
            break;
        case SimulationType::iTEBD:
            break;
    }
    ccout(1) << left  << "]";

    ccout(1) << left  << " Time: "                          << setw(10) << setprecision(2)    << fixed   << t_tot.get_age() ;

    ccout(1) << left << " Memory [";
    ccout(1) << left << "Now: "   << process_memory_in_mb("VmSize")<< " MB ";
    ccout(1) << left << "Peak: "  << process_memory_in_mb("VmPeak")<< " MB";
    ccout(1) << left << "]";


    ccout(1) << std::endl;
    t_prt.toc();
}

void class_algorithm_base::print_status_full(){
    compute_observables();
    t_prt.tic();
    std::cout << std::endl;
    std::cout << " -- Final results -- " << sim_name << std::endl;
    ccout(0)  << setw(20) << "Iterations               = " << setprecision(16) << fixed      << iteration     << std::endl;
    switch(sim_type){
        case SimulationType::iDMRG:
            ccout(0)  << setw(20) << "Energy MPO           = " << setprecision(16) << fixed      << measurement->get_energy_mpo()     << std::endl;
            ccout(0)  << setw(20) << "Energy HAM           = " << setprecision(16) << fixed      << measurement->get_energy_ham()     << std::endl;
            ccout(0)  << setw(20) << "Energy MOM           = " << setprecision(16) << fixed      << measurement->get_energy_mom()     << std::endl;
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(0)  << setw(20) << "Energy MPO           = " << setprecision(16) << fixed      << measurement->get_energy_mpo()     << std::endl;
            break;
        case SimulationType::iTEBD:
            ccout(0)  << setw(20) << "Energy HAM           = " << setprecision(16) << fixed      << measurement->get_energy_ham()     << std::endl;
            ccout(0)  << setw(20) << "Energy MOM           = " << setprecision(16) << fixed      << measurement->get_energy_mom()     << std::endl;
            break;
    }
    switch(sim_type){
        case SimulationType::iDMRG:
            ccout(0)  << setw(20) << "log₁₀ σ²(E) MPO:     = " << setprecision(6) << fixed      << log10(measurement->get_variance_mpo())  << std::endl;
            ccout(0)  << setw(20) << "log₁₀ σ²(E) HAM:     = " << setprecision(6) << fixed      << log10(measurement->get_variance_ham())  << std::endl;
            ccout(0)  << setw(20) << "log₁₀ σ²(E) MOM:     = " << setprecision(6) << fixed      << log10(measurement->get_variance_mom())  << std::endl;
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(0)  << setw(20) << "log₁₀ σ²(E) MPO:     = " << setprecision(6) << fixed      << log10(measurement->get_variance_mpo())  << std::endl;
            break;
        case SimulationType::iTEBD:
            ccout(0)  << setw(20) << "log₁₀ σ²(E) HAM:     = " << setprecision(6) << fixed      << log10(measurement->get_variance_ham())  << std::endl;
            ccout(0)  << setw(20) << "log₁₀ σ²(E) MOM:     = " << setprecision(6) << fixed      << log10(measurement->get_variance_mom())  << std::endl;
            break;
    }
    ccout(0)  << setw(20) << "Entanglement Entropy = " << setprecision(16) << fixed      << measurement->get_entanglement_entropy()   << std::endl;
    ccout(0)  << setw(20) << "χmax                 = " << setprecision(4)  << fixed      << chi_max                                   << std::endl;
    ccout(0)  << setw(20) << "χ                    = " << setprecision(4)  << fixed      << measurement->get_chi()                    << std::endl;
    ccout(0)  << setw(20) << "log₁₀ truncation:    = " << setprecision(4)  << fixed      << log10(measurement->get_truncation_error())<< std::endl;

    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(0)  << setw(20) << "Chain length         = " << setprecision(1)  << fixed      << measurement->get_chain_length()           << std::endl;
            ccout(0)  << setw(20) << "Sweep                = " << setprecision(1)  << fixed      << env_storage->get_sweeps() << std::endl;
            break;
        case SimulationType::iTEBD:
//    ccout(0)  << setw(20) << "δt:                  = " << setprecision(16) << fixed      << superblock->H->step_size << std::endl;
            break;
        default:
            break;
    }

    ccout(0)  << setw(20) << "Simulation converged : " << std::boolalpha << simulation_has_converged << std::endl;
    ccout(0)  << setw(20) << "S slope              = " << setprecision(16) << fixed      << S_slope
                          << "Converged: " << std::boolalpha << entanglement_has_converged
                          << "Saturated: " << std::boolalpha << entanglement_has_saturated
                          << std::endl;
    switch(sim_type){
        case SimulationType::iDMRG:
            ccout(0)  << setw(20) << "σ² MPO slope         = " << setprecision(16) << fixed      << V_mpo_slope  << " " << std::boolalpha << variance_mpo_has_converged << std::endl;
            ccout(0)  << setw(20) << "σ² HAM slope         = " << setprecision(16) << fixed      << V_ham_slope  << " " << std::boolalpha << variance_ham_has_converged << std::endl;
            ccout(0)  << setw(20) << "σ² MOM slope         = " << setprecision(16) << fixed      << V_mom_slope  << " " << std::boolalpha << variance_mom_has_converged << std::endl;
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
//            ccout(0)  << setw(20) << "σ² MPO slope         = " << setprecision(16) << fixed      << V_mpo_slope  << " " << std::boolalpha << variance_mpo_has_converged << std::endl;
            ccout(0)  << setw(20) << "σ² MPO slope         = " << setprecision(16) << fixed      << V_mpo_slope
                      << "Converged: " << std::boolalpha << variance_mpo_has_converged
                      << "Saturated: " << std::boolalpha << variance_mpo_has_saturated
                      << std::endl;
            break;
        case SimulationType::iTEBD:
            ccout(0)  << setw(20) << "σ² HAM slope         = " << setprecision(16) << fixed      << V_ham_slope  << " " << std::boolalpha << variance_ham_has_converged << std::endl;
            ccout(0)  << setw(20) << "σ² MOM slope         = " << setprecision(16) << fixed      << V_mom_slope  << " " << std::boolalpha << variance_mom_has_converged << std::endl;
            break;
    }

    ccout(0) << setw(20) << "Time                 = " << setw(10) << setprecision(2)    << fixed   << t_tot.get_age()  << std::endl;
    ccout(0) << setw(20) << "Peak memory          = " << process_memory_in_mb("VmPeak") << " MB" << std::endl;

    std::cout << std::endl;
    t_prt.toc();
}


void class_algorithm_base::set_profiling_labels() {
    using namespace settings::profiling;
    t_tot.set_properties(on, precision,"+Total Time              ");
    t_sto.set_properties(on, precision,"↳ Store to file          ");
    t_ste.set_properties(on, precision,"↳ Finite chain storage   ");
    t_prt.set_properties(on, precision,"↳ Printing to console    ");
    t_obs.set_properties(on, precision,"↳ Computing observables  ");
    t_sim.set_properties(on, precision,"↳+Simulation             ");
    t_evo.set_properties(on, precision,"↳ Time Evolution         ");
    t_opt.set_properties(on, precision,"↳+Optimize MPS           ");
    t_eig.set_properties(on, precision," ↳ Eigenvalue solver     ");
    t_ham.set_properties(on, precision," ↳ Build Hamiltonian     ");
    t_svd.set_properties(on, precision,"↳ SVD Truncation         ");
    t_udt.set_properties(on, precision,"↳ Update Timestep        ");
    t_env.set_properties(on, precision,"↳ Update Environments    ");
    t_con.set_properties(on, precision,"↳ Check Convergence      ");
}
