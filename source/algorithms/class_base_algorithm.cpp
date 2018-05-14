//
// Created by david on 2018-01-18.
//

#include <complex>
#include "class_base_algorithm.h"
#include <IO/class_hdf5_file.h>
#include <IO/class_hdf5_table_buffer.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_environment.h>
#include <mps_routines/class_measurement.h>
#include <mps_routines/class_hamiltonian.h>
#include <mps_routines/class_finite_chain_storage.h>
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_mps.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/class_svd_wrapper.h>

namespace s = settings;
using namespace std;
using namespace Textra;
using namespace std::complex_literals;
using Scalar = class_base_algorithm::Scalar;

class_base_algorithm::class_base_algorithm(std::shared_ptr<class_hdf5_file> hdf5_,
                                           std::string sim_name_,
                                           std::string table_name_,
                                           SimulationType sim_type_)
        :hdf5           (std::move(hdf5_)),
         sim_name       (std::move(sim_name_)),
         table_name     (std::move(table_name_)),
         sim_type       (sim_type_) {
    initialize_constants();
    set_profiling_labels();
    table_buffer = std::make_shared<class_hdf5_table_buffer>(hdf5, sim_name, table_name);
    superblock   = std::make_shared<class_superblock>();
    if(sim_type == SimulationType::fDMRG or sim_type == SimulationType::xDMRG){
        env_storage  = std::make_shared<class_finite_chain_storage>(max_length, superblock, hdf5);
        measurement  = std::make_shared<class_measurement>(superblock, env_storage, sim_type);
    }else{
        measurement  = std::make_shared<class_measurement>(superblock, sim_type);
    }

    initialize_state(settings::model::initial_state);
};


void class_base_algorithm::single_DMRG_step(long chi_max){
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();

    superblock->set_current_dimensions();
    t_opt.tic();
    superblock->MPS->theta = superblock->MPS->get_theta();
    superblock->MPS->theta = superblock->optimize_MPS(superblock->MPS->theta, s::precision::eigMaxIter,
                                                      s::precision::eigThreshold);
    t_opt.toc();
    t_svd.tic();
    superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, chi_max, s::precision::SVDThreshold);
    t_svd.toc();
    measurement->is_measured = false;
    //Reduce the hamiltonians if you are doing infinite systems:
    if(sim_type == SimulationType::iDMRG){
        superblock->HA->update_site_energy(superblock->E_one_site);
        superblock->HB->update_site_energy(superblock->E_one_site);
    }
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
    measurement->is_measured = false;
    t_sim.toc();
}

void class_base_algorithm::update_chi(){
    t_chi.tic();
    if(entropy_has_converged()) {
        switch(sim_type){
            case SimulationType::iDMRG:
            case SimulationType::fDMRG:
            case SimulationType::xDMRG:
                simulation_has_converged = chi_temp == chi_max;
                break;
            case SimulationType::iTEBD:
                simulation_has_converged = chi_temp == chi_max and delta_t <= delta_tmin;
                if (chi_temp == chi_max and delta_t > delta_tmin) {
                    delta_t = std::max(delta_tmin, delta_t * 0.5);
                    superblock->H->update_evolution_step_size(-delta_t, suzuki_order);
                }
                break;
        }
        chi_temp = chi_grow ? std::min(chi_max, chi_temp + 4) : chi_max;
    }
    if(not chi_grow){
        chi_temp = chi_max;
    }
    t_chi.toc();

}


bool class_base_algorithm::entropy_has_converged() {
    //Based on the the slope of entanglement entropy
    int last_measurement = X_vec.empty() ? 0 : X_vec.back();
    if (iteration - last_measurement < 100){return false;}
    S_vec.push_back_limited(measurement->get_entanglement_entropy());
    X_vec.push_back_limited(iteration);
    if (S_vec.size() < 5){return false;}
    double n = S_vec.size();

    double avgX = accumulate(X_vec.begin(), X_vec.end(), 0.0) / n;
    double avgY = accumulate(S_vec.begin(), S_vec.end(), 0.0) / n;

    double numerator = 0.0;
    double denominator = 0.0;

    for(int i=0; i<n; ++i){
        numerator   += (X_vec[i] - avgX) * (S_vec[i] - avgY);
        denominator += (X_vec[i] - avgX) * (X_vec[i] - avgX);
    }

    double slope = numerator / denominator;
    if (std::abs(slope) < 1e-8){
        S_vec.clear();
        X_vec.clear();
        return true;
    }else{
        return false;
    }
}



void class_base_algorithm::store_table_entry(){
    if (Math::mod(iteration, store_freq) != 0) {return;}
    compute_observables();
    t_sto.tic();
    table_buffer->emplace_back(measurement->get_chi(),
                               chi_max,
                               measurement->get_energy1(),
                               measurement->get_energy2(),
                               measurement->get_energy3(),
                               measurement->get_energy4(),
                               measurement->get_energy5(),
                               measurement->get_energy6(),
                               measurement->get_entanglement_entropy(),
                               measurement->get_variance1(),
                               measurement->get_variance2(),
                               measurement->get_variance3(),
                               measurement->get_variance4(),
                               measurement->get_variance5(),
                               measurement->get_variance6(),
                               measurement->get_truncation_error(),
                               measurement->get_parity(),
                               iteration,
                               measurement->get_chain_length(),
                               sweeps,
                               position,
                               delta_t,
                               t_tot.get_age(),
                               phys_time);
    t_sto.toc();
}

void class_base_algorithm::initialize_constants(){
    using namespace settings;
    switch(sim_type){
        case SimulationType::iDMRG:
            max_steps    = idmrg::max_steps;
            chi_max      = idmrg::chi_max   ;
            chi_grow     = idmrg::chi_grow  ;
            print_freq   = idmrg::print_freq;
            store_freq   = idmrg::store_freq;
            break;
        case SimulationType::fDMRG:
            max_length   = fdmrg::max_length;
            max_sweeps   = fdmrg::max_sweeps;
            chi_max      = fdmrg::chi_max;
            chi_grow     = fdmrg::chi_grow;
            print_freq   = fdmrg::print_freq;
            store_freq   = fdmrg::store_freq;
            break;
        case SimulationType::xDMRG:
            max_length   = xdmrg::max_length;
            max_sweeps   = xdmrg::max_sweeps;
            chi_max      = xdmrg::chi_max   ;
            chi_grow     = xdmrg::chi_grow  ;
            print_freq   = xdmrg::print_freq;
            store_freq   = xdmrg::store_freq;
            seed         = xdmrg::seed      ;
            r_strength   = xdmrg::r_strength;
            break;
        case SimulationType::iTEBD:
            max_steps    = itebd::max_steps   ;
            delta_t0     = itebd::delta_t0    ;
            delta_tmin   = itebd::delta_tmin  ;
            suzuki_order = itebd::suzuki_order;
            chi_max      = itebd::chi_max     ;
            chi_grow     = itebd::chi_grow    ;
            print_freq   = itebd::print_freq  ;
            store_freq   = itebd::store_freq  ;
            break;
    }
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

    }else if (initial_state == "random_chi"){
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
        std::cerr << "  random_chi" << std::endl;
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

    if(sim_type == SimulationType::fDMRG or sim_type == SimulationType::xDMRG ){
        position = env_storage_insert();
    }

    position = enlarge_environment();
    swap();
}



void class_base_algorithm::compute_observables(){
    t_obs.tic();
    measurement->compute_all_observables_from_superblock();
    t_obs.toc();
}


int class_base_algorithm::enlarge_environment(){
    t_env.tic();
    superblock->enlarge_environment();
    t_env.toc();
    return superblock->Lblock->size;
}
int class_base_algorithm::enlarge_environment(int direction){
    t_env.tic();
    superblock->enlarge_environment(direction);
    t_env.toc();
    return superblock->Lblock->size;
}

void class_base_algorithm::swap(){
    superblock->swap_AB();
}

int class_base_algorithm::env_storage_insert() {
    t_ste.tic();
    int position = env_storage->insert();
    t_ste.toc();
    return position;
}

void class_base_algorithm::env_storage_overwrite_MPS(){
    t_ste.tic();
    env_storage->overwrite_MPS();
    t_ste.toc();
}

int class_base_algorithm::env_storage_move(){
    t_ste.tic();
    int position = env_storage->move(direction, sweeps);
    t_ste.toc();
    return position;
}


void class_base_algorithm::print_status_update() {
    if (Math::mod(iteration, print_freq) != 0) {return;}
    if ((unsigned long)position != superblock->environment_size/2ul){return;}
    if (print_freq == 0) {return;}

    compute_observables();
    t_prt.tic();
    std::cout << setprecision(16) << fixed << left;
    ccout(1) << left  << table_name << " ";
    ccout(1) << left  << "Step: "                       << setw(10) << iteration;
    ccout(1) << left  << "E: "                          << setw(21) << setprecision(16)    << fixed   << measurement->get_energy1();
    ccout(1) << left  << "S: "                          << setw(21) << setprecision(16)    << fixed   << measurement->get_entanglement_entropy();
    ccout(1) << left  << "χmax: "                       << setw(4)  << setprecision(3)     << fixed   << chi_max;
    ccout(1) << left  << "χ: "                          << setw(4)  << setprecision(3)     << fixed   << measurement->get_chi();
    ccout(1) << left  << "log₁₀ σ²(E): "                << setw(12) << setprecision(4)     << fixed   << log10(measurement->get_variance1());
    ccout(1) << left  << " "                            << setw(12) << setprecision(4)     << fixed   << log10(measurement->get_variance2());
    ccout(1) << left  << " "                            << setw(12) << setprecision(4)     << fixed   << log10(measurement->get_variance4());
    ccout(1) << left  << "log₁₀ truncation: "           << setw(12) << setprecision(4)     << fixed   << log10(measurement->get_truncation_error());
    ccout(1) << left  << "Chain length: "               << setw(12) << setprecision(1)     << fixed   << measurement->get_chain_length();
    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
            ccout(1) << left  << "Pos: "                    << setw(6)  << position;
            ccout(1) << left  << "Dir: "                    << setw(3)  << direction;
            ccout(1) << left  << "Sweep: "                  << setw(4)  << sweeps;
            break;
        case SimulationType::iTEBD:
            ccout(1) << left  << "δt: "               << setw(13) << setprecision(12)    << fixed   << superblock->H->step_size;
            ccout(1) << left  << " conv: " << simulation_has_converged;
            break;
        default:
            break;
    }
    ccout(1) << std::endl;
    t_prt.toc();
}

void class_base_algorithm::print_status_full(){
    compute_observables();
    t_prt.tic();
    std::cout << std::endl;
    std::cout << " -- Final results -- " << sim_name << "/" << table_name << std::endl;
    ccout(0)  << setw(20) << "Energy               = " << setprecision(16) << fixed      << measurement->get_energy1()                << std::endl;
    ccout(0)  << setw(20) << "Entanglement Entropy = " << setprecision(16) << fixed      << measurement->get_entanglement_entropy()   << std::endl;
    ccout(0)  << setw(20) << "χmax                 = " << setprecision(4)  << fixed      << chi_max                                   << std::endl;
    ccout(0)  << setw(20) << "χ                    = " << setprecision(4)  << fixed      << measurement->get_chi()                    << std::endl;
    ccout(0)  << setw(20) << "log₁₀ σ²(E) MPO:     = " << setprecision(16) << fixed      << log10(measurement->get_variance1())       << std::endl;
    ccout(0)  << setw(20) << "log₁₀ σ²(E) HAM:     = " << setprecision(16) << fixed      << log10(measurement->get_variance2())       << std::endl;
    ccout(0)  << setw(20) << "log₁₀ σ²(E) GEN:     = " << setprecision(16) << fixed      << log10(measurement->get_variance4())       << std::endl;
    ccout(0)  << setw(20) << "log₁₀ truncation:    = " << setprecision(4)  << fixed      << log10(measurement->get_truncation_error())<< std::endl;
    ccout(0)  << setw(20) << "Chain length         = " << setprecision(1)  << fixed      << measurement->get_chain_length()           << std::endl;

    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
    ccout(0)  << setw(20) << "Sweep                = " << setprecision(1)  << fixed      << sweeps << std::endl;
            break;
        case SimulationType::iTEBD:
    ccout(0)  << setw(20) << "δt:                  = " << setprecision(16) << fixed      << superblock->H->step_size << std::endl;
    ccout(0)  << setw(20) << "conv:                = " << setprecision(1)  << fixed      << simulation_has_converged << std::endl;
            break;
        default:
            break;
    }

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
    t_chi.set_properties(on, precision,"↳ Update Chi             ");
}
