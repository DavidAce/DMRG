//
// Created by david on 2018-01-18.
//

#include <fstream>
#include <complex>
#include "class_algorithm_base.h"
#include <io/class_hdf5_table_buffer2.h>
#include <io/nmspc_logger.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_environment.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_state/nmspc_mps_tools.h>
#include <mps_state/class_mps_2site.h>
#include <general/nmspc_math.h>
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
    MPS_Tools::log = Logger::setLogger(sim_name,settings::console::verbosity,settings::console::timestamp);
    log->trace("Constructing class_algorithm_base");
//    ccout.set_verbosity(settings::console::verbosity);
    set_profiling_labels();
    MPS_Tools::Common::Prof::init_profiling(settings::profiling::on,settings::profiling::precision);
    log->trace("Constructing default state");
    state             = std::make_shared<class_finite_chain_state>();
    log->trace("Constructing table_profiling");
    table_profiling = std::make_unique<class_hdf5_table<class_table_profiling>>(h5pp_file, sim_name + "/measurements", "profiling", sim_name);
    log->trace("Constructing superblock");
    superblock      = std::make_shared<class_superblock>(sim_type,sim_name);
    log->trace("Writing input_file");
    h5pp_file->writeDataset(settings::input::input_file, "common/input_file");
    h5pp_file->writeDataset(settings::input::input_filename, "common/input_filename");
}



void class_algorithm_base::single_DMRG_step(eigutils::eigSetting::Ritz ritz){
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    log->trace("Starting single DMRG step");
    t_sim.tic();
    t_opt.tic();
    Eigen::Tensor<Scalar,4> theta = superblock->get_theta();
//    superblock->MPS->theta = superblock->get_theta();
    theta = superblock->optimize_MPS(theta, ritz);
    t_opt.toc();
    t_svd.tic();
    superblock->truncate_MPS(theta, sim_state.chi_temp, s::precision::SVDThreshold);
    t_svd.toc();
    //Reduce the hamiltonians if you are doing infinite systems:
//    if(sim_type == SimulationType::iDMRG){
//        superblock->E_optimal /= 2.0;
//        superblock->HA->set_reduced_energy(superblock->E_optimal);
//        superblock->HB->set_reduced_energy(superblock->E_optimal);
//    }
//    measurement->unset_measurements();
    superblock->unset_measurements();
    t_sim.toc();
    sim_state.wall_time = t_tot.get_age();
    sim_state.simu_time = t_sim.get_age();
}


void class_algorithm_base::check_convergence(){
    log->trace("Checking convergence");
    t_con.tic();
    check_convergence_entg_entropy();
    check_convergence_variance_mpo();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    update_bond_dimension();
    if(sim_state.entanglement_has_converged and
       sim_state.variance_mpo_has_converged and
       sim_state.variance_ham_has_converged and
       sim_state.variance_mom_has_converged and
       sim_state.bond_dimension_has_reached_max)
    {
        sim_state.simulation_has_converged = true;
    }
    t_con.toc();
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
        bool &has_saturated){
    //Check convergence based on slope.
    log->trace("Checking saturation using slope");


    // We want to check once every "rate" steps
    // Get the sim_state.iteration number when you last measured.
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

void class_algorithm_base::check_convergence_variance_mpo(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    log->debug("Checking convergence of variance mpo");
    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
//    compute_observables();
    check_saturation_using_slope(B_mpo_vec,
                                 V_mpo_vec,
                                 X_mpo_vec,
                                 MPS_Tools::Common::Measure::energy_variance_per_site_mpo(*superblock),
                                 sim_state.step,
                                 1,
                                 slope_threshold,
                                 V_mpo_slope,
                                 sim_state.variance_mpo_has_saturated);
//    int rewind = std::min((int)B_mpo_vec.size(), (int)measurement->get_chain_length()/4 );
//    auto b_it = B_mpo_vec.begin();
//    std::advance(b_it, -rewind);
    sim_state.variance_mpo_saturated_for = (int) count(B_mpo_vec.begin(), B_mpo_vec.end(), true);
    sim_state.variance_mpo_has_converged =  MPS_Tools::Common::Measure::energy_variance_per_site_mpo(*superblock) < threshold;

}

void class_algorithm_base::check_convergence_variance_ham(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    log->trace("Checking convergence of variance ham");

    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(B_ham_vec,
                                 V_ham_vec,
                                 X_ham_vec,
                                 MPS_Tools::Common::Measure::energy_variance_per_site_ham(*superblock),
                                 sim_state.step,
                                 1,
                                 slope_threshold,
                                 V_ham_slope,
                                 sim_state.variance_ham_has_saturated);
    sim_state.variance_ham_has_converged = MPS_Tools::Common::Measure::energy_variance_per_site_ham(*superblock) < threshold;
}

void class_algorithm_base::check_convergence_variance_mom(double threshold,double slope_threshold){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    log->trace("Checking convergence of variance mom");

    threshold       = std::isnan(threshold)       ? settings::precision::VarConvergenceThreshold : threshold;
    slope_threshold = std::isnan(slope_threshold) ? settings::precision::VarSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(B_mom_vec,
                                 V_mom_vec,
                                 X_mom_vec,
                                 MPS_Tools::Common::Measure::energy_variance_per_site_mom(*superblock),
                                 sim_state.step,
                                 1,
                                 slope_threshold,
                                 V_mom_slope,
                                 sim_state.variance_mom_has_saturated);
    sim_state.variance_mom_has_converged = MPS_Tools::Common::Measure::energy_variance_per_site_mom(*superblock) < threshold;
}

void class_algorithm_base::check_convergence_entg_entropy(double slope_threshold) {
    //Based on the the slope of entanglement middle_entanglement_entropy
    // This one is cheap to compute.
    log->debug("Checking convergence of entanglement");

    slope_threshold = std::isnan(slope_threshold) ? settings::precision::EntEntrSaturationThreshold  : slope_threshold;
    check_saturation_using_slope(BS_vec,
                                 S_vec,
                                 XS_vec,
                                 MPS_Tools::Common::Measure::current_entanglement_entropy(*superblock),
                                 sim_state.iteration,
                                 1,
                                 slope_threshold,
                                 S_slope,
                                 sim_state.entanglement_has_saturated);
    sim_state.entanglement_has_converged = sim_state.entanglement_has_saturated;
}

void class_algorithm_base::update_bond_dimension(size_t min_saturation_length){
    sim_state.chi_max = chi_max();
    if(not chi_grow() or sim_state.bond_dimension_has_reached_max or sim_state.chi_temp == chi_max() ){
        sim_state.chi_temp = chi_max();
        sim_state.bond_dimension_has_reached_max = true;
    }
    if(not sim_state.variance_mpo_has_converged
       and sim_state.variance_mpo_has_saturated
       and sim_state.variance_mpo_saturated_for >= min_saturation_length
       and sim_state.chi_temp < chi_max()){
        log->trace("Updating bond dimension");
        sim_state.chi_temp = std::min(chi_max(), sim_state.chi_temp * 2);
        log->info("New chi = {}", sim_state.chi_temp);
        clear_saturation_status();
    }
    if(sim_state.chi_temp == chi_max()){
        sim_state.bond_dimension_has_reached_max = true;
    }
}

void class_algorithm_base::clear_saturation_status(){
    log->trace("Clearing saturation status");

    BS_vec.clear();
    S_vec.clear();
    XS_vec.clear();

    B_mpo_vec.clear();
    V_mpo_vec.clear();
    X_mpo_vec.clear();
    B_ham_vec.clear();
    V_ham_vec.clear();
    X_ham_vec.clear();
    B_mom_vec.clear();
    V_mom_vec.clear();
    X_mom_vec.clear();

    sim_state.entanglement_has_saturated      = false;
    sim_state.variance_mpo_has_saturated      = false;
    sim_state.variance_ham_has_saturated      = false;
    sim_state.variance_mom_has_saturated      = false;

    sim_state.variance_mpo_saturated_for = 0;
    sim_state.variance_ham_saturated_for = 0;
    sim_state.variance_mom_saturated_for = 0;



    sim_state.entanglement_has_converged = false;
    sim_state.variance_mpo_has_converged = false;
    sim_state.variance_ham_has_converged = false;
    sim_state.variance_mom_has_converged = false;

    sim_state.bond_dimension_has_reached_max = false;
    sim_state.simulation_has_to_stop         = false;
}


void class_algorithm_base::store_algorithm_state_to_file(){
    if (settings::hdf5::storage_level < StorageLevel::LIGHT){return;}
    log->trace("Storing simulation state to file");
    t_sto.tic();
    MPS_Tools::Common::H5pp::write_algorithm_state(sim_state, *h5pp_file, sim_name);
    t_sto.toc();
}


void class_algorithm_base::store_profiling_deltas(bool force) {
    if(not force){
        if (Math::mod(sim_state.iteration, store_freq()) != 0) {return;}
        if (not state->position_is_any_edge()) {return;}
        if (not settings::profiling::on or not settings::hdf5::store_profiling){return;}
        if (settings::hdf5::storage_level < StorageLevel::NORMAL){return;}
    }

    log->trace("Storing profiling deltas");
    t_sto.tic();
    table_profiling->append_record(
            sim_state.iteration,
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
        if (Math::mod(sim_state.iteration, store_freq()) != 0) {return;}
        if (not state->position_is_any_edge()) {return;}
        if (not settings::profiling::on or not settings::hdf5::store_profiling){return;}
        if (settings::hdf5::storage_level < StorageLevel::NORMAL){return;}
    }



    log->trace("Storing profiling totals");
    table_profiling->append_record(
            sim_state.iteration,
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
//    long d    = superblock->HA->get_spin_dimension();
//    long chiA = superblock->MPS->chiA();
//    long chiB = superblock->MPS->chiB();
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
//        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = superblock->get_theta();
//        superblock->truncate_MPS(theta, 1, settings::precision::SVDThreshold);
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
//        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = superblock->get_theta();
//        superblock->truncate_MPS(theta, 1, settings::precision::SVDThreshold);
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
//        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = superblock->get_theta();
//        superblock->truncate_MPS(theta, 2, settings::precision::SVDThreshold);
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
//        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = superblock->get_theta();
//        superblock->truncate_MPS(theta, 2, settings::precision::SVDThreshold);
//    }
//
//    else if (initial_state == "rps"){
//        // Random product state
//        log->info("Initializing random product state");
//
//        //Initialize as spinors
//        theta = Textra::Matrix_to_Tensor(Eigen::MatrixXcd::Random(d*chiA,d*chiB),d,chiA,d,chiB);
//        superblock->truncate_MPS(theta, 1, settings::precision::SVDThreshold);
//
//    }else if (initial_state == "random_chi" ){
//        // Random state
//        log->info("Initializing random state with bond dimension chi = {}", chi_max());
//        sim_state.chi_temp = chi_max();
//        GA.resize(array3{d,chi_max(),chi_max()});
//        GB.resize(array3{d,chi_max(),chi_max()});
//        LA.resize(array1{chi_max()});
//        LB.resize(array1{chi_max()});
//        LC.resize(array1{chi_max()});
//        LA.setConstant(1.0/sqrt(chi_max()));
//        LB.setConstant(1.0/sqrt(chi_max()));
//        LC.setConstant(1.0/sqrt(chi_max()));
//        superblock->MPS->set_mps(LA,GA,LC,GB,LB);
//        theta = Textra::Matrix_to_Tensor(Eigen::MatrixXcd::Random(d*chi_max(),d*chi_max()),d,chi_max(),d,chi_max());
//        superblock->truncate_MPS(theta, chi_max(), settings::precision::SVDThreshold);
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
//    superblock->Lblock->set_edge_dims(*superblock->MPS, superblock->HA->MPO());
//    superblock->Rblock->set_edge_dims(*superblock->MPS, superblock->HB->MPO());
//    superblock->Lblock2->set_edge_dims(*superblock->MPS, superblock->HA->MPO());
//    superblock->Rblock2->set_edge_dims(*superblock->MPS, superblock->HB->MPO());
//
//    superblock->environment_size = superblock->Lblock->size + superblock->Rblock->size;
//
//    assert(superblock->Lblock->block.dimension(0) == superblock->MPS->chiA());
//    assert(superblock->Rblock->block.dimension(0) == superblock->MPS->chiB());
//
//
//
//    if(sim_type == SimulationType::fDMRG or sim_type == SimulationType::xDMRG ){
//        MPS_Tools::Finite::Chain::insert_superblock_to_state(*state, *superblock);
//    }else{
//    }
//
//    enlarge_environment();
//
//    if (sim_type == SimulationType::iDMRG){
//        sim_state.iteration = (int)superblock->Lblock->size;
//    }
//    swap();
//}

void class_algorithm_base::reset_full_mps_to_random_product_state(const std::string parity) {
    log->trace("Resetting MPS to random product state");
    if (state->get_length() != (size_t)num_sites()) throw std::range_error("System size mismatch");
    sim_state.iteration = state->reset_sweeps();

    // Randomize chain
    *state = MPS_Tools::Finite::Ops::set_random_product_state(*state,parity);
    MPS_Tools::Finite::Ops::rebuild_superblock(*state,*superblock);
    clear_saturation_status();
}


void class_algorithm_base::compute_observables(const class_superblock & superblock){
    log->trace("Starting all measurements on current superblock");
    t_sim.tic();
    t_obs.tic();
    superblock.do_all_measurements();
    t_obs.toc();
    t_sim.toc();
}

void class_algorithm_base::compute_observables(const class_finite_chain_state & state){
    log->trace("Starting all measurements on current superblock");
    t_sim.tic();
    t_obs.tic();
    state.do_all_measurements();
    t_obs.toc();
    t_sim.toc();

}

void class_algorithm_base::enlarge_environment(){
    log->trace("Enlarging environments");
    t_sim.tic();
    t_env.tic();
    superblock->enlarge_environment();
    t_env.toc();
    t_sim.toc();
}

void class_algorithm_base::enlarge_environment(int direction){
    log->trace("Enlarging environment in direction: {}",direction );
    t_sim.tic();
    t_env.tic();
    superblock->enlarge_environment(direction);
    t_env.toc();
    t_sim.toc();
}

void class_algorithm_base::swap(){
    log->trace("Swap AB sites on superblock");
    t_sim.tic();
    superblock->swap_AB();
    t_sim.toc();
}

void class_algorithm_base::insert_superblock_to_chain() {
    log->trace("Insert superblock into chain");
    t_sim.tic();
    t_ste.tic();
    auto new_position = MPS_Tools::Finite::Chain::insert_superblock_to_state(*state, *superblock);
    superblock->set_positions(new_position);
    t_ste.toc();
    t_sim.toc();
}

void class_algorithm_base::copy_superblock_mps_to_chain(){
    log->trace("Copy superblock mps to chain");
    t_sim.tic();
    t_ste.tic();
    MPS_Tools::Finite::Chain::copy_superblock_mps_to_state(*state, *superblock);
    t_ste.toc();
    t_sim.toc();
}

void class_algorithm_base::copy_superblock_mpo_to_chain(){
    log->trace("Copy superblock mpo to chain");
    t_sim.tic();
    t_ste.tic();
    MPS_Tools::Finite::Chain::copy_superblock_mpo_to_state(*state, *superblock);
    t_ste.toc();
    t_sim.toc();
}

void class_algorithm_base::copy_superblock_env_to_chain(){
    log->trace("Copy superblock env to chain");
    t_sim.tic();
    t_ste.tic();
    MPS_Tools::Finite::Chain::copy_superblock_env_to_state(*state, *superblock);
    t_ste.toc();
    t_sim.toc();
}

void class_algorithm_base::copy_superblock_to_chain(){
    log->trace("Copy superblock to chain");
    t_sim.tic();
    t_ste.tic();
    MPS_Tools::Finite::Chain::copy_superblock_to_state(*state, *superblock);
    t_ste.toc();
    t_sim.toc();
}


void class_algorithm_base::move_center_point(){
    log->trace("Moving center point ");
    t_sim.tic();
    t_ste.tic();
    MPS_Tools::Finite::Chain::move_center_point(*state,*superblock);
    t_ste.toc();
    t_sim.toc();
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
    if (Math::mod(sim_state.iteration, print_freq()) != 0) {return;}
//    if (not state->position_is_the_middle()) {return;}
    if (print_freq() == 0) {return;}
    compute_observables(*superblock);
    t_prt.tic();
    std::stringstream report;
    report << setprecision(16) << fixed << left;
    report << left  << sim_name << " ";
    report << left  << "Iter: "                       << setw(6) << sim_state.iteration;
    report << left  << "E: ";

    switch(sim_type) {
        case SimulationType::iDMRG:
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_mpo.value();
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_ham.value();
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_mom.value();
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_mpo.value();
            break;
        case SimulationType::iTEBD:
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_ham.value();
            report << setw(21) << setprecision(16)    << fixed   << superblock->measurements.energy_per_site_mom.value();
            break;
    }

    if (sim_type == SimulationType::xDMRG){
        report << left  << " ε: "<< setw(8) << setprecision(4) << fixed << sim_state.energy_dens;
    }

    report << left  << "log₁₀ σ²(E): ";
    switch(sim_type) {
        case SimulationType::iDMRG:
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_mpo.value());
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_ham.value());
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_mom.value());
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << setw(14) << setprecision(6)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_mpo.value());
            break;
        case SimulationType::iTEBD:
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_ham.value());
            report << setw(12) << setprecision(4)    << fixed   << std::log10(superblock->measurements.energy_variance_per_site_mom.value());
            break;
    }


    report << left  << "S: "                          << setw(21) << setprecision(16)    << fixed   << superblock->measurements.current_entanglement_entropy.value();
    report << left  << "χmax: "                       << setw(4)  << setprecision(3)     << fixed   << chi_max();
    report << left  << "χ: "                          << setw(4)  << setprecision(3)     << fixed   << superblock->measurements.bond_dimension.value();
    report << left  << "log₁₀ trunc: "                << setw(10) << setprecision(4)     << fixed   << std::log10(superblock->measurements.truncation_error.value());
    report << left  << "Sites: "                      << setw(6)  << setprecision(1)     << fixed   << superblock->measurements.length.value();
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << left  << "@ site: "                    << setw(5)  << state->get_position();
            //ccout(1) << left  << "Dir: "                    << setw(3)  << state->get_direction();
            //ccout(1) << left  << "Sweep: "                  << setw(4)  << state->get_sweeps();
            break;
        case SimulationType::iTEBD:
//            ccout(1) << left  << "δt: "               << setw(13) << setprecision(12)    << fixed   << delta_t;
            break;
        default:
            break;
    }
    report << left  << " Convergence [";
    switch(sim_type){
        case SimulationType::iDMRG:
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_state.entanglement_has_converged;
            report << left  << " σ²-"  << std::boolalpha << setw(6) << sim_state.variance_mpo_has_converged;
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << left  << " σ²-"  << std::boolalpha << setw(6) << sim_state.variance_mpo_has_converged;
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_state.entanglement_has_converged;
            break;
        case SimulationType::iTEBD:
            report << left  << " S-"  << std::boolalpha << setw(6) << sim_state.entanglement_has_converged;
            break;
    }
    report << left  << "]";
    report << left  << " Saturation [";
    switch(sim_type){
        case SimulationType::iDMRG:
            report << left  << " σ²- " << setw(2) << sim_state.variance_mpo_saturated_for << " steps";
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_state.entanglement_has_saturated;
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            report << left  << " σ²-" << setw(2) << sim_state.variance_mpo_saturated_for << " steps";
            break;
        case SimulationType::iTEBD:
            report << left  << " S-"   << std::boolalpha << setw(6) << sim_state.entanglement_has_saturated;
            break;
    }
    report << left  << "]";
    report << left  << " Time: "                          << setw(10) << setprecision(2)    << fixed   << t_tot.get_age() ;
    report << left << " Memory [";
    report << left << "Rss: "     << process_memory_in_mb("VmRSS")<< " MB ";
    report << left << "RssPeak: "  << process_memory_in_mb("VmHWM")<< " MB ";
    report << left << "VmPeak: "  << process_memory_in_mb("VmPeak")<< " MB";
    report << left << "]";
    log->info(report.str());
    t_prt.toc();
}

void class_algorithm_base::print_status_full(){
    compute_observables(*superblock);

    using namespace MPS_Tools::Common::Measure;
    t_prt.tic();
    log->info("--- Final results  --- {} ---", sim_name);
    log->info("Iterations            = {:<16d}"    , sim_state.iteration);
    switch(sim_type){
        case SimulationType::iDMRG:
            log->info("Energy MPO            = {:<16.16f}" , superblock->measurements.energy_per_site_mpo.value());
            log->info("Energy HAM            = {:<16.16f}" , superblock->measurements.energy_per_site_ham.value());
            log->info("Energy MOM            = {:<16.16f}" , superblock->measurements.energy_per_site_mom.value());
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("Energy MPO            = {:<16.16f}" , superblock->measurements.energy_per_site_mpo.value());
            break;
        case SimulationType::iTEBD:
            log->info("Energy HAM            = {:<16.16f}" , superblock->measurements.energy_per_site_ham.value());
            log->info("Energy MOM            = {:<16.16f}" , superblock->measurements.energy_per_site_mom.value());
            break;
    }
    switch(sim_type){
        case SimulationType::iDMRG:
            log->info("log₁₀ σ²(E) MPO       = {:<16.16f}" , superblock->measurements.energy_per_site_mpo.value());
            log->info("log₁₀ σ²(E) HAM       = {:<16.16f}" , superblock->measurements.energy_per_site_ham.value());
            log->info("log₁₀ σ²(E) MOM       = {:<16.16f}" , superblock->measurements.energy_per_site_mom.value());
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("log₁₀ σ²(E) MPO       = {:<16.16f}" , log10(superblock->measurements.energy_variance_per_site_mpo.value()));
            break;
        case SimulationType::iTEBD:
            log->info("log₁₀ σ²(E) HAM       = {:<16.16f}" , superblock->measurements.energy_per_site_ham.value());
            log->info("log₁₀ σ²(E) MOM       = {:<16.16f}" , superblock->measurements.energy_per_site_mom.value());
            break;
    }

    log->info("Entanglement Entropy  = {:<16.16f}" , superblock->measurements.current_entanglement_entropy.value());
    log->info("χmax                  = {:<16d}"    , chi_max()                                            );
    log->info("χ                     = {:<16d}"    , superblock->measurements.bond_dimension.value()      );
    log->info("log₁₀ truncation:     = {:<16.16f}" , log10(superblock->measurements.truncation_error.value()));

    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("Chain length          = {:<16d}"    , superblock->measurements.length.value());
            log->info("Sweep                 = {:<16d}"    , state->get_sweeps());
            break;
        case SimulationType::iTEBD:
            log->info("δt                    = {:<16.16f}" , sim_state.delta_t);
            break;
        default:
            break;
    }

    log->info("Simulation converged  = {:<}"    , sim_state.simulation_has_converged);

    switch(sim_type){
        case SimulationType::iDMRG:
            log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_state.entanglement_has_converged, sim_state.entanglement_has_saturated);
            log->info("σ² MPO slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mpo_slope ,sim_state.variance_mpo_has_converged, sim_state.variance_mpo_has_saturated);
            log->info("σ² HAM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_ham_slope ,sim_state.variance_ham_has_converged, sim_state.variance_ham_has_saturated);
            log->info("σ² MOM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mom_slope ,sim_state.variance_mom_has_converged, sim_state.variance_mom_has_saturated);
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            log->info("σ² MPO slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mpo_slope ,sim_state.variance_mpo_has_converged, sim_state.variance_mpo_has_saturated);
            break;
        case SimulationType::iTEBD:
            log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_state.entanglement_has_converged, sim_state.entanglement_has_saturated);
            log->info("σ² HAM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_ham_slope ,sim_state.variance_ham_has_converged, sim_state.variance_ham_has_saturated);
            log->info("σ² MOM slope          = {:<16.16f} | Converged : {} \t\t Saturated: {}" , V_mom_slope ,sim_state.variance_mom_has_converged, sim_state.variance_mom_has_saturated);
            break;
    }
    log->info("S slope               = {:<16.16f} | Converged : {} \t\t Saturated: {}" , S_slope,sim_state.entanglement_has_converged, sim_state.entanglement_has_saturated);
    log->info("Time                  = {:<16.16f}" , t_tot.get_age());
    log->info("Peak memory           = {:<6.1f} MB" , process_memory_in_mb("VmPeak"));
    t_prt.toc();
}


void class_algorithm_base::set_profiling_labels() {
    using namespace settings::profiling;
    t_tot.set_properties(true, precision,"+Total Time              ");
    t_sto.set_properties(on,   precision,"↳ Store to file          ");
    t_ste.set_properties(on,   precision,"↳ Finite chain storage   ");
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
