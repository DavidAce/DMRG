//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <io/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_state/class_superblock.h>
#include <mps_state/class_mps_2site.h>
#include <mps_state/class_finite_chain_state.h>
#include <mps_tools/nmspc_mps_tools.h>
#include <mps_tools/finite/opt.h>
#include <general/nmspc_math.h>
#include <general/nmspc_random_numbers.h>
#include <general/nmspc_quantum_mechanics.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <general/nmspc_tensor_extra.h>
#include <h5pp/h5pp.h>
#include <spdlog/spdlog.h>

#include "class_xDMRG.h"

template<typename Scalar, auto rank>
void throw_if_has_imaginary_part(const Eigen::Tensor<Scalar,rank> &tensor, double threshold = 1e-14) {
    Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> vector (tensor.data(),tensor.size());
    if constexpr (std::is_same<Scalar, std::complex<double>>::value){
        auto imagSum = vector.imag().cwiseAbs().sum();
        if (imagSum > threshold){
            std::cout << vector << std::endl;
            throw std::runtime_error("Has imaginary part. Sum: " + std::to_string(imagSum));
        }
    }
}




using namespace std;
using namespace Textra;

class_xDMRG::class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_finite(std::move(h5ppFile_), "xDMRG",SimulationType::xDMRG) {
    mpstools::finite::chain::initialize_state(*state, settings::model::model_type, settings::model::symmetry, settings::xdmrg::num_sites);
    min_saturation_length = 1;// * (int)(1.0 * settings::xdmrg::num_sites);
    max_saturation_length = 2;// * (int)(2.0 * settings::xdmrg::num_sites);
    settings::xdmrg::min_sweeps = std::max(settings::xdmrg::min_sweeps, 1+(size_t)(std::log2(chi_max())/2));
}







void class_xDMRG::run_preprocessing() {

    log->info("Starting {} preprocessing", sim_name);
    sim_state.energy_dens_target = settings::xdmrg::energy_density_target;
    sim_state.energy_dens_window = settings::xdmrg::energy_density_window;


    find_energy_range();
    mpstools::finite::print::print_hamiltonians(*state);
    log->info("Finished {} preprocessing", sim_name);
}

void class_xDMRG::run_simulation()    {
    log->info("Starting {} simulation", sim_name);
    while(true) {
        single_DMRG_step();
        store_table_entry_progress();
        store_table_entry_site_state();
        store_profiling_totals();
        store_state_and_measurements_to_file();

        check_convergence();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if (state->position_is_any_edge())
        {
            if (sim_state.iteration >= settings::xdmrg::max_sweeps) {stop_reason = StopReason::MAX_STEPS; break;}
            if (sim_state.simulation_has_converged)                 {stop_reason = StopReason::CONVERGED; break;}
            if (sim_state.simulation_has_to_stop)                   {stop_reason = StopReason::SATURATED; break;}
        }

        update_bond_dimension(min_saturation_length);
        enlarge_environment(state->get_direction());
        move_center_point();
        sim_state.iteration = state->get_sweeps();
        sim_state.step++;
        sim_state.position = state->get_position();
        log->trace("Finished step {}, iteration {}",sim_state.step,sim_state.iteration);
    }
    switch(stop_reason){
        case StopReason::MAX_STEPS : log->info("Finished {} simulation -- reason: MAX_ITERS",sim_name) ;break;
        case StopReason::CONVERGED : log->info("Finished {} simulation -- reason: CONVERGED",sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }
}




void class_xDMRG::single_DMRG_step()
/*!
 * \fn void single_DMRG_step()
 */
{

    t_sim.tic();
    t_opt.tic();
    log->trace("Starting single xDMRG step");
    Eigen::Tensor<Scalar,4> theta;
    auto dims = superblock->dimensions();
    auto eigsize = dims[0]*dims[1]*dims[2]*dims[3];

    using namespace  mpstools::finite::opt;

    // Table

    // Mode / Space |   FULL        PARTIAL     DIRECT
    // ---------------------------------------------------
    // OVERLAP      |   FO          FP          DV
    // VARIANCE     |   FV          FV          DV


    auto optMode  =  OptMode::OVERLAP;
    optMode  = sim_state.iteration   >= 2  ?  OptMode::VARIANCE : optMode;

    auto optSpace =  OptSpace::FULL;
    optSpace = eigsize >= 2*2*16*16     ? OptSpace::PARTIAL : optSpace;
    optSpace = eigsize >= 2*2*32*32     ? OptSpace::DIRECT  : optSpace;
    optSpace = sim_state.iteration >=
            settings::xdmrg::min_sweeps ? OptSpace::DIRECT  : optSpace;

    auto optType = superblock->isReal() ? OptType::REAL : OptType::CPLX;
    mpstools::finite::multisite::compute_best_jump(*state,optSpace);
    std::tie(theta, sim_state.energy_now) = mpstools::finite::opt::find_optimal_excited_state(*superblock,sim_state,optMode, optSpace,optType);
    sim_state.energy_dens = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);

//    Textra::subtract_phase(theta);
    t_opt.toc();
    log->trace("Truncating theta");

    t_svd.tic();
    superblock->truncate_MPS(theta, sim_state.chi_temp, settings::precision::SVDThreshold);
    t_svd.toc();

    superblock->unset_measurements();
    state->unset_measurements();

    t_sim.toc();

    sim_state.wall_time = t_tot.get_age();
    sim_state.simu_time = t_sim.get_age();

}


void class_xDMRG::check_convergence(){

    t_sim.tic();
    t_con.tic();

    check_convergence_variance_mpo();
    check_convergence_entg_entropy();

    if (sim_state.iteration < settings::xdmrg::min_sweeps){
        clear_saturation_status();
    }


    bool outside_of_window = std::abs(sim_state.energy_dens - sim_state.energy_dens_target)  > sim_state.energy_dens_window;
    if (outside_of_window
        and (   sim_state.iteration >= 2
                or sim_state.variance_mpo_saturated_for > min_saturation_length
                or sim_state.variance_mpo_has_converged)
        )
    {
        log->info("Resetting to product state -- saturated outside of energy window. Energy density: {}, Energy window: {} --> {}",sim_state.energy_dens, sim_state.energy_dens_window, std::min(1.05*sim_state.energy_dens_window, 0.5) );
        sim_state.energy_dens_window = std::min(1.05*sim_state.energy_dens_window, 0.5);
        int counter = 0;
        while(outside_of_window){
            reset_full_mps_to_random_product_state("sx");
            sim_state.energy_now  = mpstools::finite::measure::energy_per_site_mpo(*state);
            sim_state.energy_dens = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);
            outside_of_window = std::abs(sim_state.energy_dens - sim_state.energy_dens_target)  >= sim_state.energy_dens_window;
            counter++;
            if (counter % 10 == 0) {
                log->info("Resetting to product state -- can't find state in energy window.  Increasing energy window: {} --> {}", sim_state.energy_dens_window, std::min(1.05*sim_state.energy_dens_window, 0.5) );
                sim_state.energy_dens_window = std::min(1.05*sim_state.energy_dens_window, 0.5);
            }
        }
        log->info("Energy initial (per site) = {} | density = {} | retries = {}", sim_state.energy_now, sim_state.energy_dens,counter );
        clear_saturation_status();
        projected_during_saturation  = false;
        sim_state.energy_ubound      = sim_state.energy_target + sim_state.energy_dens_window * (sim_state.energy_max-sim_state.energy_min);
        sim_state.energy_lbound      = sim_state.energy_target - sim_state.energy_dens_window * (sim_state.energy_max-sim_state.energy_min);
    }



    if (state->position_is_any_edge()
        and sim_state.variance_mpo_has_saturated
        and not outside_of_window
        and not projected_during_saturation)
    {
        log->info("Projecting to {} due to saturation", settings::model::symmetry);
        *state = mpstools::finite::ops::get_closest_parity_state(*state,settings::model::symmetry);
        mpstools::finite::ops::rebuild_superblock(*state,*superblock);
//        clear_saturation_status();
        projected_during_saturation = true;
    }




    if(     sim_state.variance_mpo_has_converged
        and sim_state.entanglement_has_converged)
    {
        log->debug("Simulation has converged");
        sim_state.simulation_has_converged = true;
    }

    if (        sim_state.variance_mpo_has_saturated
            and sim_state.entanglement_has_saturated
            and sim_state.bond_dimension_has_reached_max
            and sim_state.variance_mpo_saturated_for > max_saturation_length)
    {
        log->debug("Simulation has to stop");
        sim_state.simulation_has_to_stop = true;
    }



    t_con.toc();
    t_sim.toc();
}




void class_xDMRG::find_energy_range() {
    log->trace("Finding energy range");
    assert(state->get_length() == settings::xdmrg::num_sites);
    size_t max_sweeps_during_f_range = 4;
    sim_state.iteration = state->reset_sweeps();

    // Find energy minimum
    while(true) {
        single_DMRG_step(eigutils::eigSetting::Ritz::SR);
        copy_superblock_to_chain();         //Needs to occurr after update_MPS...
        print_status_update();
        // It's important not to perform the last step.
        // That last state would not get optimized
        if(sim_state.iteration >= max_sweeps_during_f_range
            or superblock->measurements.energy_variance_per_site_mpo < 1e-8)
        {break;}
        enlarge_environment(state->get_direction());
        move_center_point();
        sim_state.iteration = state->get_sweeps();

    }
    compute_observables(*superblock);
    sim_state.energy_min = superblock->measurements.energy_per_site_mpo.value();

    reset_full_mps_to_random_product_state("sx");
    // Find energy maximum
    while(true) {
        single_DMRG_step(eigutils::eigSetting::Ritz::LR);
        copy_superblock_to_chain();         //Needs to occurr after update_MPS...
        print_status_update();
        // It's important not to perform the last step.
        // That last state would not get optimized
        if(sim_state.iteration >= max_sweeps_during_f_range
           or superblock->measurements.energy_variance_per_site_mpo < 1e-8)
        {break;}
        enlarge_environment(state->get_direction());
        move_center_point();
        sim_state.iteration = state->get_sweeps();
    }
    compute_observables(*superblock);
    sim_state.energy_max         = superblock->measurements.energy_per_site_mpo.value();
    sim_state.energy_now         = superblock->measurements.energy_per_site_mpo.value();
    sim_state.energy_target      = sim_state.energy_min    + sim_state.energy_dens_target  * (sim_state.energy_max - sim_state.energy_min);
    sim_state.energy_ubound      = sim_state.energy_target + sim_state.energy_dens_window  * (sim_state.energy_max - sim_state.energy_min);
    sim_state.energy_lbound      = sim_state.energy_target - sim_state.energy_dens_window  * (sim_state.energy_max - sim_state.energy_min);
    sim_state.energy_dens        = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);
    log->info("Energy minimum (per site) = {}", sim_state.energy_min);
    log->info("Energy maximum (per site) = {}", sim_state.energy_max);
    log->info("Energy target  (per site) = {}", sim_state.energy_target);
    int counterA = 0;
    int counterB = 0;
    bool outside_of_window = std::abs(sim_state.energy_dens - sim_state.energy_dens_target)  >= sim_state.energy_dens_window;
    while(outside_of_window){
        reset_full_mps_to_random_product_state("sx");
        sim_state.energy_now  = mpstools::finite::measure::energy_per_site_mpo(*state);
        sim_state.energy_dens = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);
        outside_of_window = std::abs(sim_state.energy_dens - sim_state.energy_dens_target)  >= sim_state.energy_dens_window;

        counterA++;
        counterB++;
        if(counterA >= 100){
            counterA = 0;
            if(sim_state.energy_dens_window >= 0.5){break;}
            sim_state.energy_dens_window = std::min(1.2*sim_state.energy_dens_window, 0.5);
            sim_state.energy_ubound       = sim_state.energy_target +  sim_state.energy_dens_window*(sim_state.energy_max-sim_state.energy_min);
            sim_state.energy_lbound       = sim_state.energy_target -  sim_state.energy_dens_window*(sim_state.energy_max-sim_state.energy_min);
        }
    }
    superblock->isReal();
    log->info("Energy initial (per site) = {} | density = {} | retries = {}", sim_state.energy_now, sim_state.energy_dens,counterB );


}

void class_xDMRG::store_state_and_measurements_to_file(bool force){
    if(not force){
        if (not settings::hdf5::save_progress){return;}
        if (Math::mod(sim_state.iteration, settings::xdmrg::store_freq) != 0) {return;}
        if (not state->position_is_any_edge()) {return;}
        if (settings::xdmrg::store_freq == 0){return;}
        if (settings::hdf5::storage_level <= StorageLevel::NONE){return;}
    }
    state->unset_measurements();
    compute_observables(*state);
    log->trace("Storing all measurements to file");
    t_sto.tic();
    h5pp_file->writeDataset(false, sim_name + "/simOK");
    mpstools::finite::io::write_all_measurements(*state, *h5pp_file, sim_name);
    mpstools::finite::io::write_closest_parity_projection(*state, *h5pp_file, sim_name, settings::model::symmetry);
    //  Write the wavefunction (this is only defined for short enough chain ( L < 14 say)
    if(settings::xdmrg::store_wavefn){
        h5pp_file->writeDataset(mpstools::finite::measure::mps_wavefn(*state), sim_name + "/state/psi");
    }
    mpstools::finite::io::write_all_state(*state, *h5pp_file, sim_name);
    t_sto.toc();
    store_algorithm_state_to_file();
    h5pp_file->writeDataset(true, sim_name + "/simOK");
}


void class_xDMRG::store_table_entry_progress(bool force){
    if(not force) {
        if (Math::mod(sim_state.iteration, settings::xdmrg::store_freq) != 0) { return; }
        if (not state->position_is_the_middle_any_direction()) { return; }
        if (settings::xdmrg::store_freq == 0) { return; }
        if (settings::hdf5::storage_level < StorageLevel::NORMAL){return;}

    }
    compute_observables(*state);
    log->trace("Storing table entry to file");
    t_sto.tic();
    table_dmrg->append_record(
            sim_state.iteration,
            state->get_length(),
            state->get_position(),
            state->bond_dimension(),
            settings::xdmrg::chi_max,
            state->measurements.energy_per_site_mpo.value(),
            sim_state.energy_min,
            sim_state.energy_max,
            sim_state.energy_target,
            state->measurements.energy_variance_per_site_mpo.value(),
            state->measurements.entanglement_entropy.value(),
            state->measurements.truncation_error.value(),
            t_tot.get_age());
    t_sto.toc();
}





bool   class_xDMRG::sim_on()    {return settings::xdmrg::on;}
long   class_xDMRG::chi_max()   {return settings::xdmrg::chi_max;}
size_t class_xDMRG::num_sites() {return settings::xdmrg::num_sites;}
size_t class_xDMRG::store_freq(){return settings::xdmrg::store_freq;}
size_t class_xDMRG::print_freq(){return settings::xdmrg::print_freq;}
bool   class_xDMRG::chi_grow()  {return settings::xdmrg::chi_grow;}



