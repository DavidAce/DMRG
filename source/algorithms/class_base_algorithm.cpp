//
// Created by david on 2018-01-18.
//

#include <complex>
#include "class_base_algorithm.h"
#include <IO/class_hdf5_file.h>
#include <IO/class_hdf5_table_buffer2.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_hamiltonian.h>
#include <mps_routines/class_finite_chain_sweeper.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_mps.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/class_svd_wrapper.h>
#include <algorithms/table_types.h>

namespace s = settings;
using namespace std;
using namespace Textra;
using namespace std::complex_literals;
using Scalar = class_base_algorithm::Scalar;

class_base_algorithm::class_base_algorithm(std::shared_ptr<class_hdf5_file> hdf5_,
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
};


void class_base_algorithm::single_DMRG_step(long chi_max){
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();

    superblock->set_current_dimensions();
    t_opt.tic();
    superblock->MPS->theta = superblock->MPS->get_theta();
    superblock->MPS->theta = superblock->optimize_MPS(superblock->MPS->theta);
    t_opt.toc();
    t_svd.tic();
    superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, chi_max, s::precision::SVDThreshold);
    t_svd.toc();
    //Reduce the hamiltonians if you are doing infinite systems:
    if(sim_type == SimulationType::iDMRG){
        superblock->HA->update_site_energy(superblock->E_one_site);
        superblock->HB->update_site_energy(superblock->E_one_site);
    }
    superblock->set_current_dimensions();
    measurement->set_not_measured();
    t_sim.toc();
}


void class_base_algorithm::single_TEBD_step(long chi_max){
/*!
 * \fn single_iTEBD_step(class_superblock &superblock)
 * \brief infinite Time evolving block decimation.
 * \param superblock A class containing MPS, environment and Hamiltonian MPO objects.
 */
    t_sim.tic();
    for (auto &U: superblock->H->U){
        t_evo.tic();
        superblock->set_current_dimensions();
        superblock->MPS->theta = superblock->evolve_MPS(superblock->MPS->get_theta() ,U);
        t_evo.toc();

        t_svd.tic();
        superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, chi_max, s::precision::SVDThreshold);
        t_svd.toc();

        if (&U != &superblock->H->U.back()) {
            superblock->swap_AB();        }
    }
    superblock->set_current_dimensions();
    measurement->set_not_measured();
    t_sim.toc();
}

void class_base_algorithm::check_convergence_overall(){
    t_con.tic();
    check_convergence_entanglement();
    check_convergence_variance_mpo();
    check_convergence_variance_ham();
    check_convergence_variance_mom();
    check_convergence_bond_dimension();
    if(entanglement_has_converged and
       variance_mpo_has_converged and
       variance_ham_has_converged and
       variance_mom_has_converged and
       bond_dimension_has_converged)
    {
        simulation_has_converged = true;
    }
    t_con.toc();
}

void class_base_algorithm::check_convergence_using_slope(
                                   std::list<double> &Y_vec,
                                   std::list<int> &X_vec,
                                   double new_data,
                                   int    rate,
                                   double tolerance,
                                   double &slope,
                                   bool &has_converged){
    //Check convergence based on slope.
    if (not measurement->has_been_measured()){return;}
    // We want to check once every "rate" steps
    unsigned long max_data_points = 10;
    unsigned long min_data_points = 4;
    // Get the iteration number when you last measured.
    int last_measurement = X_vec.empty() ? 0 : X_vec.back();
    // If the measurement happened less than rate iterations ago, return.
    if (iteration - last_measurement < rate){return;}

    // It's time to check. Insert current numbers and remove old ones
    Y_vec.push_back(new_data);
    X_vec.push_back(iteration);
    if (Y_vec.size() < min_data_points){return;}
    if (Y_vec.size() > max_data_points){
        Y_vec.pop_front();
        X_vec.pop_front();
    }

    double n = Y_vec.size();
    double avgX = accumulate(X_vec.begin(), X_vec.end(), 0.0) / n;
    double avgY = accumulate(Y_vec.begin(), Y_vec.end(), 0.0) / n;

    double numerator = 0.0;
    double denominator = 0.0;

    auto x_it = X_vec.begin();
    auto y_it = Y_vec.begin();
    auto v_end = Y_vec.end();
    while(y_it != v_end){
        numerator   += (*x_it - avgX) * (*y_it - avgY);
        denominator += (*x_it - avgX) * (*x_it - avgX);
        y_it++;
        x_it++;
    }

    slope = numerator / denominator;
    if (std::abs(slope) < tolerance) {
        has_converged = true;
    }
}

void class_base_algorithm::check_convergence_variance_mpo(){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    check_convergence_using_slope(V_mpo_vec,
                                  X_mpo_vec,
                                  measurement->get_variance_mpo(),
                                  1,
                                  settings::precision::eigThreshold,
                                  V_mpo_slope,
                                  variance_mpo_has_converged);
}

void class_base_algorithm::check_convergence_variance_ham(){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    check_convergence_using_slope(V_ham_vec,
                                  X_ham_vec,
                                  measurement->get_variance_ham(),
                                  1,
                                  settings::precision::eigThreshold,
                                  V_ham_slope,
                                  variance_ham_has_converged);

}

void class_base_algorithm::check_convergence_variance_mom(){
    //Based on the the slope of the variance
    // We want to check every time we can because the variance is expensive to compute.
    check_convergence_using_slope(V_mom_vec,
                                  X_mom_vec,
                                  measurement->get_variance_mom(),
                                  1,
                                  settings::precision::eigThreshold,
                                  V_mom_slope,
                                  variance_mom_has_converged);

}

void class_base_algorithm::check_convergence_entanglement() {
    //Based on the the slope of entanglement entanglement_entropy
    // This one is cheap to compute but doesnt need to be done every step.
    check_convergence_using_slope(S_vec,
                                  X2_vec,
                                  measurement->get_entanglement_entropy(),
                                  1,
                                  settings::precision::SVDThreshold,
                                  S_slope,
                                  entanglement_has_converged);
}

void class_base_algorithm::check_convergence_bond_dimension(){
    if(not env_storage->position_is_the_middle()){return;}
    if(not chi_grow or bond_dimension_has_converged or chi_temp == chi_max ){
        chi_temp = chi_max;
        bond_dimension_has_converged = true;
    }else{
        if(variance_mpo_has_converged or entanglement_has_converged){
            chi_temp = std::min(chi_max, chi_temp + 4);
            clear_convergence_checks();
        }
        if(chi_temp == chi_max){
            bond_dimension_has_converged = true;
        }
    }
}

void class_base_algorithm::clear_convergence_checks(){
    S_vec.clear();
    V_mpo_vec.clear();
    X_mpo_vec.clear();
    X2_vec.clear();
    simulation_has_converged        = false;
    bond_dimension_has_converged    = false;
    entanglement_has_converged      = false;
    variance_mpo_has_converged      = false;
}


void class_base_algorithm::store_profiling_to_file() {
//    if (Math::mod(iteration, store_freq) != 0) {return;}
//    t_sto.tic();
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


void class_base_algorithm::initialize_state(std::string initial_state ) {
    //Set the size and initial values for the MPS and environments
    //Choose between GHZ, W, Random, Product state (up, down, etc), None, etc...
    long d  = superblock->d;
    std::srand((unsigned int) 1);

    if(initial_state == "upup"){
        std::cout << "Initializing Up-Up-state  |up,up>" << std::endl;
        superblock->MPS->GA.resize(array3{d,1,1});
        superblock->MPS->GA.setZero();
        superblock->MPS->LA.resize(array1{1});
        superblock->MPS->LA.setConstant(1.0);
        superblock->MPS->GB.resize(array3{d,1,1});
        superblock->MPS->GB.setZero();
        superblock->MPS->LB.resize(array1{1});
        superblock->MPS->LB.setConstant(1.0);
        superblock->MPS->LB_left = superblock->MPS->LB;
        superblock->MPS->GA(0, 0, 0) = 1;
        superblock->MPS->GA(1, 0, 0) = 0;
        superblock->MPS->GB(0, 0, 0) = 1;
        superblock->MPS->GB(1, 0, 0) = 0;
        superblock->MPS->theta = superblock->MPS->get_theta();
    }else if(initial_state == "updown"){
        std::cout << "Initializing Up down -state  |up,down>" << std::endl;
        superblock->MPS->GA.resize(array3{d,1,1});
        superblock->MPS->GA.setZero();
        superblock->MPS->LA.resize(array1{1});
        superblock->MPS->LA.setConstant(1.0);
        superblock->MPS->GB.resize(array3{d,1,1});
        superblock->MPS->GB.setZero();
        superblock->MPS->LB.resize(array1{1});
        superblock->MPS->LB.setConstant(1.0);
        superblock->MPS->LB_left = superblock->MPS->LB;
        superblock->MPS->GA(0, 0, 0) = 1;
        superblock->MPS->GA(1, 0, 0) = 0;
        superblock->MPS->GB(0, 0, 0) = 0;
        superblock->MPS->GB(1, 0, 0) = 1;
        superblock->MPS->theta = superblock->MPS->get_theta();

    }else if(initial_state == "ghz"){
        std::cout << "Initializing GHZ-state" << std::endl;
        // GHZ state (|up,up> + |down, down > ) /sqrt(2)
        superblock->MPS->GA.resize(array3{d,1,2});
        superblock->MPS->GA.setZero();
        superblock->MPS->LA.resize(array1{2});
        superblock->MPS->LA.setConstant(1.0/std::sqrt(2));
        superblock->MPS->GB.resize(array3{d,2,1});
        superblock->MPS->GB.setZero();
        superblock->MPS->LB.resize(array1{1});
        superblock->MPS->LB.setConstant(1.0);
        superblock->MPS->LB_left = superblock->MPS->LB;
        // GA^0 = (1,0)
        // GA^1 = (0,1)
        // GB^0 = (1,0)^T
        // GB^1 = (0,1)^T
        superblock->MPS->GA(0, 0, 0) = 1;
        superblock->MPS->GA(0, 0, 1) = 0;
        superblock->MPS->GA(1, 0, 0) = 0;
        superblock->MPS->GA(1, 0, 1) = 1;
        superblock->MPS->GB(0, 0, 0) = 1;
        superblock->MPS->GB(0, 1, 0) = 0;
        superblock->MPS->GB(1, 0, 0) = 0;
        superblock->MPS->GB(1, 1, 0) = 1;
        superblock->MPS->theta = superblock->MPS->get_theta();

    }else if(initial_state == "w"){
        std::cout << "Initializing W-state" << std::endl;
        // W state (|up,down> + |down, up > ) /sqrt(2)
        superblock->MPS->GA.resize(array3{d,1,2});
        superblock->MPS->GA.setZero();
        superblock->MPS->LA.resize(array1{2});
        superblock->MPS->LA.setConstant(1.0/std::sqrt(2));
        superblock->MPS->GB.resize(array3{d,2,1});
        superblock->MPS->GB.setZero();
        superblock->MPS->LB.resize(array1{1});
        superblock->MPS->LB.setConstant(1.0);
        superblock->MPS->LB_left = superblock->MPS->LB;
        // GA^0 = (1,0)
        // GA^1 = (0,1)
        // GB^0 = (0,1)^T
        // GB^1 = (1,0)^T
        superblock->MPS->GA(0, 0, 0) = 1;
        superblock->MPS->GA(0, 0, 1) = 0;
        superblock->MPS->GA(1, 0, 0) = 0;
        superblock->MPS->GA(1, 0, 1) = 1;
        superblock->MPS->GB(0, 0, 0) = 0;
        superblock->MPS->GB(0, 1, 0) = 1;
        superblock->MPS->GB(1, 0, 0) = 1;
        superblock->MPS->GB(1, 1, 0) = 0;
        superblock->MPS->theta = superblock->MPS->get_theta();

    }

    else if (initial_state == "rps"){
        // Random product state
        std::cout << "Initializing random product state" << std::endl;
        superblock->MPS->GA.resize(array3{d,1,1});
        superblock->MPS->GA.setZero();
        superblock->MPS->LA.resize(array1{1});
        superblock->MPS->LA.setConstant(1.0);
        superblock->MPS->GB.resize(array3{d,1,1});
        superblock->MPS->GB.setZero();
        superblock->MPS->LB.resize(array1{1});
        superblock->MPS->LB.setConstant(1.0);
        superblock->MPS->LB_left = superblock->MPS->LB;
        //Initialize to as spinors
        auto r1 = rn::uniform_complex_1();
        auto r2 = rn::uniform_complex_1();
        superblock->MPS->GA(1, 0, 0) = r1.real();
        superblock->MPS->GA(0, 0, 0) = r1.imag();
        superblock->MPS->GB(1, 0, 0) = r2.real();
        superblock->MPS->GB(0, 0, 0) = r2.imag();
        superblock->MPS->theta = superblock->MPS->get_theta();

    }else if (initial_state == "random_chi" ){
        // Random state
        std::cout << "Initializing random state with bond dimension chi = " << chi_max << std::endl;
        chi_temp = chi_max;
        superblock->MPS->GA.resize(array3{d,chi_max,chi_max});
        superblock->MPS->GA.setZero();
        superblock->MPS->LA.resize(array1{chi_max});
        superblock->MPS->LA.setConstant(1.0/sqrt(chi_max));
        superblock->MPS->GB.resize(array3{d,chi_max,chi_max});
        superblock->MPS->GB.setZero();
        superblock->MPS->LB.resize(array1{chi_max});
        superblock->MPS->LB.setConstant(1.0/sqrt(chi_max));
        superblock->MPS->LB_left = superblock->MPS->LB;
        superblock->MPS->theta = Textra::Matrix_to_Tensor(Eigen::MatrixXd::Random(d*chi_max,d*chi_max).cast<Scalar>(),d,chi_max,d,chi_max);
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

    //Get a properly normalized initial state.
    superblock->set_current_dimensions();
    superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, chi_max, settings::precision::SVDThreshold);



    //Reset the environment blocks to the correct dimensions
    superblock->Lblock->set_edge_dims(superblock->MPS, superblock->HA->MPO);
    superblock->Rblock->set_edge_dims(superblock->MPS, superblock->HB->MPO);
    superblock->Lblock2->set_edge_dims(superblock->MPS, superblock->HA->MPO);
    superblock->Rblock2->set_edge_dims(superblock->MPS, superblock->HB->MPO);

    superblock->environment_size = superblock->Lblock->size + superblock->Rblock->size;

    assert(superblock->Lblock->block.dimension(0) == superblock->MPS->LB_left.dimension(0));
    assert(superblock->Rblock->block.dimension(0) == superblock->MPS->LB.dimension(0));
    assert(superblock->MPS->GA.dimension(2) == superblock->MPS->LA.dimension(0));
    assert(superblock->MPS->GA.dimension(1) == superblock->MPS->LB_left.dimension(0));
    assert(superblock->MPS->GB.dimension(1) == superblock->MPS->LA.dimension(0));
    assert(superblock->MPS->GB.dimension(2) == superblock->MPS->LB.dimension(0));


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


void class_base_algorithm::compute_observables(){
    t_obs.tic();
    measurement->compute_all_observables_from_superblock();
    t_obs.toc();
}


void class_base_algorithm::enlarge_environment(){
    t_env.tic();
    superblock->enlarge_environment();
    t_env.toc();
}

void class_base_algorithm::enlarge_environment(int direction){
    t_env.tic();
    superblock->enlarge_environment(direction);
    t_env.toc();
}

void class_base_algorithm::swap(){
    superblock->swap_AB();
}

void class_base_algorithm::env_storage_insert() {
    t_ste.tic();
    env_storage->insert();
    t_ste.toc();
}

void class_base_algorithm::env_storage_overwrite_MPS(){
    t_ste.tic();
    env_storage->overwrite_MPS();
    t_ste.toc();
}

void class_base_algorithm::env_storage_move(){
    t_ste.tic();
    env_storage->move();
    t_ste.toc();
}


void class_base_algorithm::print_status_update() {
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
    ccout(1) << left  << "log₁₀ truncation: "           << setw(12) << setprecision(4)     << fixed   << log10(measurement->get_truncation_error());
    ccout(1) << left  << "Chain length: "               << setw(12) << setprecision(1)     << fixed   << measurement->get_chain_length();
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(1) << left  << "Pos: "                    << setw(6)  << env_storage->get_position();
            ccout(1) << left  << "Dir: "                    << setw(3)  << env_storage->get_direction();
            ccout(1) << left  << "Sweep: "                  << setw(4)  << env_storage->get_sweeps();
            break;
        case SimulationType::iTEBD:
            ccout(1) << left  << "δt: "               << setw(13) << setprecision(12)    << fixed   << superblock->H->step_size;
            break;
        default:
            break;
    }
    ccout(1) << left  << " - Convergence ";
    ccout(1) << left  << " S-"  << std::boolalpha << entanglement_has_converged;
    ccout(1) << left  << " χ-"  << std::boolalpha << bond_dimension_has_converged;
    switch(sim_type){
        case SimulationType::iDMRG:
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(1) << left  << " σ²- "  << std::boolalpha << variance_mpo_has_converged;
            break;
        case SimulationType::iTEBD:
            break;
    }


    ccout(1) << std::endl;
    t_prt.toc();
}

void class_base_algorithm::print_status_full(){
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
    ccout(0)  << setw(20) << "δt:                  = " << setprecision(16) << fixed      << superblock->H->step_size << std::endl;
            break;
        default:
            break;
    }

    ccout(0)  << setw(20) << "Simulation converged : " << std::boolalpha << simulation_has_converged << std::endl;
    ccout(0)  << setw(20) << "S slope              = " << setprecision(16) << fixed      << S_slope  << " " << std::boolalpha << entanglement_has_converged << std::endl;
    switch(sim_type){
        case SimulationType::iDMRG:
            ccout(0)  << setw(20) << "σ² MPO slope         = " << setprecision(16) << fixed      << V_mpo_slope  << " " << std::boolalpha << variance_mpo_has_converged << std::endl;
            ccout(0)  << setw(20) << "σ² HAM slope         = " << setprecision(16) << fixed      << V_ham_slope  << " " << std::boolalpha << variance_ham_has_converged << std::endl;
            ccout(0)  << setw(20) << "σ² MOM slope         = " << setprecision(16) << fixed      << V_mom_slope  << " " << std::boolalpha << variance_mom_has_converged << std::endl;
            break;
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(0)  << setw(20) << "σ² MPO slope         = " << setprecision(16) << fixed      << V_mpo_slope  << " " << std::boolalpha << variance_mpo_has_converged << std::endl;
            break;
        case SimulationType::iTEBD:
            ccout(0)  << setw(20) << "σ² HAM slope         = " << setprecision(16) << fixed      << V_ham_slope  << " " << std::boolalpha << variance_ham_has_converged << std::endl;
            ccout(0)  << setw(20) << "σ² MOM slope         = " << setprecision(16) << fixed      << V_mom_slope  << " " << std::boolalpha << variance_mom_has_converged << std::endl;
            break;
    }
    ccout(0)  << setw(20) << "χ                    = " << setprecision(1)  << fixed      << chi_temp << " " << std::boolalpha << bond_dimension_has_converged << std::endl;
    std::cout << std::endl;
    t_prt.toc();
}

void class_base_algorithm::set_profiling_labels() {
    using namespace settings::profiling;
    t_tot.set_properties(on, precision,"+Total Time              ");
    t_sto.set_properties(on, precision,"↳ Store to file          ");
    t_ste.set_properties(on, precision,"↳ Finite chain storage   ");
    t_prt.set_properties(on, precision,"↳ Printing to console    ");
    t_obs.set_properties(on, precision,"↳ Computing observables  ");
    t_sim.set_properties(on, precision,"↳+Simulation             ");
    t_evo.set_properties(on, precision,"↳ Time Evolution         ");
    t_opt.set_properties(on, precision,"↳ Optimize MPS           ");
    t_svd.set_properties(on, precision,"↳ SVD Truncation         ");
    t_udt.set_properties(on, precision,"↳ Update Timestep        ");
    t_env.set_properties(on, precision,"↳ Update Environments    ");
    t_con.set_properties(on, precision,"↳ Check Convergence      ");
}
