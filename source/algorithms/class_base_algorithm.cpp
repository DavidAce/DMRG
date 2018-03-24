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
#include <mps_routines/class_mpo.h>
#include <mps_routines/class_mps.h>
#include <general/nmspc_math.h>

namespace s = settings;
using namespace std;
using namespace Textra;
using namespace std::complex_literals;
using Scalar = class_algorithm_base::Scalar;

class_algorithm_base::class_algorithm_base(std::shared_ptr<class_hdf5_file> hdf5_,
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
    measurement  = std::make_shared<class_measurement>(superblock, sim_type);
};


void class_algorithm_base::single_DMRG_step(long chi_max){
/*!
 * \fn void single_DMRG_step(class_superblock &superblock)
 */
    t_sim.tic();

    superblock->set_current_dimensions();
    superblock->MPS->theta = superblock->MPS->get_theta();
    t_eig.tic();
    superblock->MPS->theta = superblock->optimize_MPS(superblock->MPS->theta, s::precision::eigSteps, s::precision::eigThreshold);
    t_eig.toc();

    t_svd.tic();
    superblock->MPS->theta = superblock->truncate_MPS(superblock->MPS->theta, chi_max, s::precision::SVDThreshold);
    t_svd.toc();
    measurement->is_measured = false;
    t_sim.toc();
}


void class_algorithm_base::single_TEBD_step(long chi_max){
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

void class_algorithm_base::update_chi(){
    t_chi.tic();
    if(entropy_has_converged()) {
        switch(sim_type){
            case SimulationType::iDMRG:
            case SimulationType::fDMRG:
            case SimulationType::xDMRG:
            case SimulationType::FES_iDMRG:
                simulation_has_converged = chi_temp == chi_max;
                break;
            case SimulationType::iTEBD:
            case SimulationType ::FES_iTEBD:
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


bool class_algorithm_base::entropy_has_converged() {
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



void class_algorithm_base::store_table_entry(){
    if (Math::mod(iteration, store_freq) != 0) {return;}
    t_sto.tic();
    measurement->do_full_measurement();
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
                               iteration,
                               measurement->get_chain_length(),
                               sweeps,
                               position,
                               delta_t,
                               t_tot.get_age(),
                               phys_time);
    t_sto.toc();
}

void class_algorithm_base::initialize_constants(){
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
        case SimulationType::FES_iTEBD:
            max_steps    = fes_itebd::max_steps;
            delta_t0     = fes_itebd::delta_t0;
            delta_tmin   = fes_itebd::delta_tmin;
            suzuki_order = fes_itebd::suzuki_order;
            chi_min      = fes_itebd::chi_min;
            chi_max      = fes_itebd::chi_max;
            chi_num      = fes_itebd::chi_num;
            chi_grow     = fes_itebd::chi_grow;
            print_freq   = fes_itebd::print_freq;
            store_freq   = fes_itebd::store_freq;
            break;
        case SimulationType::FES_iDMRG:
            max_steps    = fes_idmrg::max_steps ;
            chi_min      = fes_idmrg::chi_min   ;
            chi_max      = fes_idmrg::chi_max   ;
            chi_num      = fes_idmrg::chi_num   ;
            chi_grow     = fes_idmrg::chi_grow  ;
            print_freq   = fes_idmrg::print_freq;
            store_freq   = fes_idmrg::store_freq;
            break;
    }
}


void class_algorithm_base::set_profiling_labels() {
    using namespace settings::profiling;
    t_tot.set_properties(on, precision,"+Total Time              ");
    t_sto.set_properties(on, precision,"↳ Store to file          ");
    t_ste.set_properties(on, precision,"↳ Store Environment      ");
    t_obs.set_properties(on, precision,"↳ Printing measurement   ");
    t_sim.set_properties(on, precision,"↳+Simulation             ");
    t_evo.set_properties(on, precision," ↳ Time Evolution        ");
    t_eig.set_properties(on, precision," ↳ Eig. decomp.          ");
    t_svd.set_properties(on, precision," ↳ SVD Truncation        ");
    t_udt.set_properties(on, precision," ↳ Update Timestep       ");
    t_env.set_properties(on, precision," ↳ Update Environments   ");
    t_chi.set_properties(on, precision," ↳ Update Chi            ");
//    t_mps.set_properties(on, precision," ↳ Update MPS            ");
}

int class_algorithm_base::enlarge_environment(){
    t_env.tic();
    superblock->enlarge_environment();
    t_env.toc();
    return superblock->Lblock->size;
}
int class_algorithm_base::enlarge_environment(int direction){
    t_env.tic();
    superblock->enlarge_environment(direction);
    t_env.toc();
    return superblock->Lblock->size;
}

void class_algorithm_base::swap(){
    superblock->swap_AB();
}
//


void class_algorithm_base::print_status_update() {
    if (Math::mod(iteration, print_freq) != 0) {return;}
    if (2+(2*position) != superblock->chain_length){return;}
    if (print_freq == 0) {return;}

    t_obs.tic();
    measurement->do_full_measurement();
    std::cout << setprecision(16) << fixed << left;
    ccout(1) << left  << table_name << " ";
    ccout(1) << left  << "Step: "                       << setw(10) << iteration;
    ccout(1) << left  << "E: "                          << setw(21) << setprecision(16)    << fixed   << measurement->get_energy1();
    ccout(1) << left  << "S: "                          << setw(21) << setprecision(16)    << fixed   << measurement->get_entanglement_entropy();
    ccout(1) << left  << "χmax: "                       << setw(4)  << setprecision(3)     << fixed   << chi_max;
    ccout(1) << left  << "χ: "                          << setw(4)  << setprecision(3)     << fixed   << measurement->get_chi();
    ccout(1) << left  << "log₁₀ σ²(E): "                << setw(12) << setprecision(4)     << fixed   << log10(measurement->get_variance2());
    ccout(1) << left  << " "                            << setw(12) << setprecision(4)     << fixed   << log10(measurement->get_variance3());
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
        case SimulationType ::FES_iTEBD:
            ccout(1) << left  << "δt: "               << setw(13) << setprecision(12)    << fixed   << superblock->H->step_size;
            ccout(1) << left  << " conv: " << simulation_has_converged;
            break;
        default:
            break;
    }
    ccout(1) << std::endl;
    t_obs.toc();
}

void class_algorithm_base::print_status_full(){
    std::cout << std::endl;
    std::cout << " -- Final results -- " << sim_name << "/" << table_name << std::endl;
    ccout(0)  << setw(20) << "Energy               = " << setprecision(16) << fixed      << measurement->get_energy1()                << std::endl;
    ccout(0)  << setw(20) << "Entanglement Entropy = " << setprecision(16) << fixed      << measurement->get_entanglement_entropy()   << std::endl;
    ccout(0)  << setw(20) << "χmax                 = " << setprecision(4)  << fixed      << chi_max                                   << std::endl;
    ccout(0)  << setw(20) << "χ                    = " << setprecision(4)  << fixed      << measurement->get_chi()                    << std::endl;
    ccout(0)  << setw(20) << "log₁₀ σ²(E):         = " << setprecision(16) << fixed      << log10(measurement->get_variance2())       << std::endl;
    ccout(0)  << setw(20) << "log₁₀ truncation:    = " << setprecision(4)  << fixed      << log10(measurement->get_truncation_error())<< std::endl;
    ccout(0)  << setw(20) << "Chain length         = " << setprecision(1)  << fixed      << measurement->get_chain_length()           << std::endl;

    switch(sim_type){
        case SimulationType::fDMRG:
        case SimulationType::xDMRG:
    ccout(0)  << setw(20) << "Sweep                = " << setprecision(1)  << fixed      << sweeps << std::endl;
            break;
        case SimulationType::iTEBD:
        case SimulationType ::FES_iTEBD:
    ccout(0)  << setw(20) << "δt:                  = " << setprecision(16) << fixed      << superblock->H->step_size << std::endl;
    ccout(0)  << setw(20) << "conv:                = " << setprecision(1)  << fixed      << simulation_has_converged << std::endl;
            break;
        default:
            break;
    }

    std::cout << std::endl;
}

