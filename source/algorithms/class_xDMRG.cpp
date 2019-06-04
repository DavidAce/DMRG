//
// Created by david on 2018-02-09.
//


#include <iomanip>
#include <io/class_hdf5_table_buffer2.h>
#include <sim_parameters/nmspc_sim_settings.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_mps_2site.h>
#include <mps_routines/class_finite_chain_state.h>
#include <mps_routines/nmspc_mps_tools.h>
#include <mps_routines/mps_tools/finite/opt.h>
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

//Print xDMRG modes nicely
std::ostream& operator<<(std::ostream& str, class_xDMRG::xDMRG_Mode const& mode) {
    switch (mode){
        case class_xDMRG::xDMRG_Mode::KEEP_BEST_OVERLAP :
            str << "KEEP BEST OVERLAP";
            break;
        case class_xDMRG::xDMRG_Mode::FULL_EIG_OPT :
            str << "FULL EIG OPT";
            break;
        case class_xDMRG::xDMRG_Mode::PARTIAL_EIG_OPT :
            str << "PARTIAL EIG OPT";
            break;
        case class_xDMRG::xDMRG_Mode::DIRECT_OPT :
            str << "DIRECT OPT";
            break;
    }
    return str;
}



using namespace std;
using namespace Textra;

class_xDMRG::class_xDMRG(std::shared_ptr<h5pp::File> h5ppFile_)
        : class_algorithm_base(std::move(h5ppFile_), "xDMRG",SimulationType::xDMRG) {
    table_xdmrg       = std::make_unique<class_hdf5_table<class_table_dmrg>>        (h5ppFile, sim_name + "/measurements", "simulation_progress",sim_name);
    table_xdmrg_chain = std::make_unique<class_hdf5_table<class_table_finite_chain>>(h5ppFile, sim_name + "/measurements", "simulation_progress_full_chain",sim_name);
    MPS_Tools::Finite::Chain::initialize_state(*state,settings::model::model_type,settings::model::symmetry, settings::xdmrg::num_sites, settings::model::seed_init_mpo, settings::model::seed_init_mps);
    MPS_Tools::Finite::Chain::copy_state_to_superblock(*state,*superblock);
    min_saturation_length = 1 * (int)(1.0 * settings::xdmrg::num_sites);
    max_saturation_length = 1 * (int)(2.0 * settings::xdmrg::num_sites);
    settings::xdmrg::min_sweeps = std::max(settings::xdmrg::min_sweeps, 1+(int)(std::log2(chi_max())/2));
}




void class_xDMRG::run()
/*!
 * \brief Dispatches xDMRG stages.
 * This function manages the stages of simulation differently depending on whether
 * the data already existed in hdf5 storage or not.
 *
 * There can be two main scenarios that split into cases:
 * 1) The hdf5 file existed already and contains
 *      a) nothing recognizeable (previous crash?)       -- run full simulation from scratch.
 *      b) a converged simulation but no MPS             -- run full simulation from scratch.
 *      c) a not-yet-converged MPS                       -- resume simulation, reset the number of sweeps first.
 *      d) a converged MPS                               -- not much to do... run postprocessing
 * 2) The hdf5 file did not exist                        -- run full simulation from scratch.

 *
 */
{
    if (!settings::xdmrg::on) { return; }
    t_tot.tic();

    if (h5ppFile->getCreateMode() == h5pp::CreateMode::OPEN){
        // This is case 1
        bool fileOK;
        h5ppFile->readDataset(fileOK, "common/fileOK");
        bool simOK = h5ppFile->linkExists(sim_name + "/simOK");
        bool mpsOK = h5ppFile->linkExists(sim_name + "/state/full/mps");
//        h5ppFile->print_contents_of_group(sim_name);

        if (not simOK){
            //Case 1 a -- run full simulation from scratch.
            log->trace("Case 1a");
            run_preprocessing();
            run_simulation();
        }else if(simOK and not mpsOK){
            // Case 1 b
            log->trace("Case 1b");
            run_preprocessing();
            run_simulation();
        }else if(simOK and mpsOK){
            // We can go ahead and load the state from hdf5
            log->trace("Loading MPS from file");
            try{
                MPS_Tools::Finite::H5pp::load_from_hdf5(*state, *superblock, sim_state, *h5ppFile, sim_name);
            }
            catch(std::exception &ex){
                log->error("Failed to load from hdf5: {}", ex.what());
                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
            }
            catch(...){log->error("Unknown error when trying to resume from file.");}

            bool convergence_was_reached;
            h5ppFile->readDataset(convergence_was_reached,sim_name + "/sim_state/simulation_has_converged");
            if(not convergence_was_reached){
                // Case 1 c -- resume simulation, reset the number of sweeps first.
                log->trace("Case 1c");
                settings::xdmrg::max_sweeps += state->get_sweeps();
                run_simulation();

            }else {
                // Case 1 d -- not much else to do.. redo postprocessing for good measure.
                log->trace("Case 1d");
            }
        }
    }else {
        // This is case 2
        log->trace("Case 2");
        run_preprocessing();
        run_simulation();
    }
    run_postprocessing();
    print_status_full();
    print_profiling();
    t_tot.toc();
}


void class_xDMRG::run_preprocessing() {

    log->info("Starting {} preprocessing", sim_name);
    find_energy_range();
    MPS_Tools::Finite::Print::print_hamiltonians(*state);
    log->info("Finished {} preprocessing", sim_name);
}

void class_xDMRG::run_simulation()    {
    log->info("Starting {} simulation", sim_name);
    while(true) {
        single_xDMRG_step();
        MPS_Tools::Finite::Chain::copy_superblock_to_state(*state, *superblock);
        store_table_entry_progress();
        store_table_entry_site_state();
        store_table_entry_profiling();
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
        case StopReason::MAX_STEPS : log->info("Finished {} simulation -- reason: MAX_STEPS",sim_name) ;break;
        case StopReason::CONVERGED : log->info("Finished {} simulation -- reason: CONVERGED",sim_name) ;break;
        case StopReason::SATURATED : log->info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: log->info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }
}


void class_xDMRG::run_postprocessing(){
    log->info("Running {} postprocessing",sim_name);
    MPS_Tools::Finite::Debug::check_integrity(*state,*superblock,sim_state);

//    MPS_Tools::Finite::Debug::check_normalization_routine(*state);

    state->set_measured_false();
    superblock->set_measured_false();
    log->info("Storing state and measurements to file");
    store_state_and_measurements_to_file(true);
    log->info("Finished {} postprocessing",sim_name);

}


void class_xDMRG::single_xDMRG_step()
/*!
 * \fn void single_DMRG_step()
 */
{

    t_sim.tic();
    t_opt.tic();
    log->trace("Starting single xDMRG step");
    Eigen::Tensor<Scalar,4> theta;
    auto dims = superblock->dimensions();
    auto eigsize = dims[1] * dims[3];

    using namespace  MPS_Tools::Finite::Opt;

    // Table

    // Mode / Space |   FULL        PARTIAL     DIRECT
    // ---------------------------------------------------
    // OVERLAP      |   FO          FP          DV
    // VARIANCE     |   FV          FV          DV


    auto optMode  =  OptMode::OVERLAP;
    optMode  = sim_state.iteration   >= settings::xdmrg::min_sweeps             ?  OptMode::VARIANCE : optMode;
    optMode  = sim_state.iteration   >= 1 and
               superblock->measurements.energy_variance_per_site_mpo  < 1e-4    ?  OptMode::VARIANCE : optMode;



    auto optSpace =  OptSpace::FULL;
    optSpace = eigsize >  16*16                     ? OptSpace::PARTIAL : optSpace;
    optSpace = eigsize >= 32*32                     ? OptSpace::DIRECT  : optSpace;
    optSpace = optMode == OptMode::VARIANCE         ? OptSpace::DIRECT  : optSpace;


    optMode   = OptMode::VARIANCE;
    optSpace  = OptSpace::GUIDED;


    std::tie(theta, sim_state.energy_now) = MPS_Tools::Finite::Opt::find_optimal_excited_state(*superblock,sim_state,optMode, optSpace);
    sim_state.energy_dens = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);


    t_opt.toc();
    log->trace("Truncating theta");

    t_svd.tic();
    superblock->truncate_MPS(theta, sim_state.chi_temp, settings::precision::SVDThreshold);
    t_svd.toc();

    superblock->set_measured_false();
    t_sim.toc();

}


void class_xDMRG::check_convergence(){

    t_sim.tic();
    t_con.tic();
    check_convergence_variance_mpo();

    if (sim_state.iteration <= settings::xdmrg::min_sweeps){
        clear_saturation_status();
    }


    bool outside_of_window = std::abs(sim_state.energy_dens - settings::xdmrg::energy_density)  >= settings::xdmrg::energy_window;
    if ((sim_state.variance_mpo_has_saturated or sim_state.variance_mpo_has_converged)
        and sim_state.variance_mpo_saturated_for > min_saturation_length
        and outside_of_window)
    {
        std::cout << "ROOOOOOLL THE DICE" << std::endl;
        std::cout <<  " |eps-0.5|       = " << std::abs(sim_state.energy_dens - settings::xdmrg::energy_density) << std::endl;
        std::cout <<  "||eps-0.5| - w | = " << std::abs(std::abs(sim_state.energy_dens - settings::xdmrg::energy_density) - settings::xdmrg::energy_window)  << std::endl;
        int counter = 0;
        while(outside_of_window){
            reset_full_mps_to_random_product_state("sx",settings::model::seed_init_mps++);
            sim_state.energy_now = MPS_Tools::Common::Measure::energy_per_site_mpo(*superblock);
            sim_state.energy_dens = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);
            outside_of_window = std::abs(sim_state.energy_dens - settings::xdmrg::energy_density)  >= settings::xdmrg::energy_window;
            counter++;
        }
        log->info("Energy initial (per site) = {} | density = {} | retries = {}", sim_state.energy_now, sim_state.energy_dens,counter );
        projected_during_saturation      = false;

    }



    if (sim_state.variance_mpo_saturated_for > min_saturation_length
        and not projected_during_saturation)
    {
        *state = MPS_Tools::Finite::Ops::get_closest_parity_state(*state,settings::model::symmetry);
        MPS_Tools::Finite::Ops::rebuild_superblock(*state,*superblock);
        clear_saturation_status();
        projected_during_saturation = true;
    }




    if(sim_state.variance_mpo_has_converged)
    {
        log->debug("Simulation has converged");
        sim_state.simulation_has_converged = true;
    }

    else if (sim_state.variance_mpo_has_saturated and
             sim_state.bond_dimension_has_reached_max and
             sim_state.variance_mpo_saturated_for > max_saturation_length
            )
    {
        log->debug("Simulation has to stop");
        sim_state.simulation_has_to_stop = true;

    }



    t_con.toc();
    t_sim.toc();
}


void class_xDMRG::initialize_chain() {
    while(true){
        insert_superblock_to_chain();
        if (superblock->environment_size + 2ul < (unsigned long) settings::xdmrg::num_sites) {
            enlarge_environment();
            swap();
        } else {
            break;
        }
    }
}


void class_xDMRG::find_energy_range() {
    log->trace("Finding energy range");
    assert(state->get_length() == (size_t)settings::xdmrg::num_sites);
    int max_sweeps_during_f_range = 4;
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
    compute_observables();
    sim_state.energy_min = superblock->measurements.energy_per_site_mpo;
    sim_state.iteration = state->reset_sweeps();

    reset_full_mps_to_random_product_state("sx",settings::model::seed_init_mps++);
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
    compute_observables();
    sim_state.energy_max = superblock->measurements.energy_per_site_mpo;

    sim_state.iteration      = state->reset_sweeps();
    sim_state.energy_target  = sim_state.energy_min    + settings::xdmrg::energy_density * (sim_state.energy_max-sim_state.energy_min);
    sim_state.energy_ubound  = sim_state.energy_target + settings::xdmrg::energy_window  * (sim_state.energy_max-sim_state.energy_min);
    sim_state.energy_lbound  = sim_state.energy_target - settings::xdmrg::energy_window  * (sim_state.energy_max-sim_state.energy_min);

    sim_state.energy_now = superblock->E_optimal / state->get_length();
    log->info("Energy minimum (per site) = {}", sim_state.energy_min);
    log->info("Energy maximum (per site) = {}", sim_state.energy_max);
    log->info("Energy target  (per site) = {}", sim_state.energy_target);
    int counterA = 0;
    int counterB = 0;
    while(sim_state.energy_now < sim_state.energy_lbound or sim_state.energy_now > sim_state.energy_ubound){
        reset_full_mps_to_random_product_state("sx",settings::model::seed_init_mps++);
        sim_state.energy_now = MPS_Tools::Common::Measure::energy_per_site_mpo(*superblock);
        sim_state.energy_dens = (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min);
        counterA++;
        counterB++;
        if(counterA >= 100){
            counterA = 0;
            if(settings::xdmrg::energy_window >= 0.5){break;}
            settings::xdmrg::energy_window *= 2;
            sim_state.energy_ubound  = sim_state.energy_target + settings::xdmrg::energy_window*(sim_state.energy_max-sim_state.energy_min);
            sim_state.energy_lbound  = sim_state.energy_target - settings::xdmrg::energy_window*(sim_state.energy_max-sim_state.energy_min);
        }
    }
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
    log->trace("Storing storing mps to file");
    t_sto.tic();
    state->do_all_measurements();
    MPS_Tools::Finite::H5pp::write_all_measurements(*state,*h5ppFile,sim_name);
    MPS_Tools::Finite::H5pp::write_closest_parity_projection(*state, *h5ppFile, sim_name, settings::model::symmetry);
    //  Write the wavefunction (this is only defined for short enough chain ( L < 14 say)
    if(settings::xdmrg::store_wavefn){
        h5ppFile->writeDataset(MPS_Tools::Finite::Measure::mps_wavefn(*state), sim_name + "/state/psi");
    }
    MPS_Tools::Finite::H5pp::write_all_state(*state, *h5ppFile, sim_name);
    h5ppFile->writeDataset(false, "/common/fileOK");
    h5ppFile->writeDataset(false, sim_name + "/simOK");
    t_sto.toc();
    store_algorithm_state_to_file();
}


void class_xDMRG::store_table_entry_progress(bool force){
    if(not force) {
        if (Math::mod(sim_state.iteration, settings::xdmrg::store_freq) != 0) { return; }
        if (not state->position_is_the_middle_any_direction()) { return; }
        if (settings::xdmrg::store_freq == 0) { return; }
        if (settings::hdf5::storage_level < StorageLevel::FULL){return;}

    }
    compute_observables();
    log->trace("Storing table entry to file");
    t_sto.tic();
    table_xdmrg->append_record(
            sim_state.iteration,
            state->get_length(),
            state->get_position(),
            superblock->measurements.bond_dimension,
            settings::xdmrg::chi_max,
            superblock->measurements.energy_per_site_mpo,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            sim_state.energy_min,
            sim_state.energy_max,
            sim_state.energy_target,
            superblock->measurements.energy_variance_per_site_mpo,
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            superblock->measurements.current_entanglement_entropy,
            superblock->measurements.truncation_error,
            t_tot.get_age());
    t_sto.toc();
}



void class_xDMRG::store_table_entry_site_state(bool force){
    if (not force) {
        if (Math::mod(sim_state.iteration, settings::xdmrg::store_freq) != 0) { return; }
        if (settings::xdmrg::store_freq == 0) { return; }
//        if (not state->position_is_the_middle()) { return; }
        if (settings::hdf5::storage_level < StorageLevel::FULL){return;}
    }

    log->trace("Storing chain_entry to file");
    t_sto.tic();
    table_xdmrg_chain->append_record(
            sim_state.iteration,
            state->get_length(),
            state->get_position(),
            MPS_Tools::Common::Measure::bond_dimension(*superblock),
            sim_state.energy_now,
            MPS_Tools::Common::Measure::current_entanglement_entropy(*superblock),
            MPS_Tools::Common::Measure::truncation_error(*superblock)
    );
    t_sto.toc();
}


long   class_xDMRG::chi_max()   {return settings::xdmrg::chi_max;}
int    class_xDMRG::num_sites() {return settings::xdmrg::num_sites;}
int    class_xDMRG::store_freq(){return settings::xdmrg::store_freq;}
int    class_xDMRG::print_freq(){return settings::xdmrg::print_freq;}
bool   class_xDMRG::chi_grow()  {return settings::xdmrg::chi_grow;}



void class_xDMRG::print_profiling(){
    if (settings::profiling::on) {
        log->trace("Printing profiling information (tot)");
        t_tot.print_time_w_percent();
        t_sto.print_time_w_percent(t_tot);
        t_ste.print_time_w_percent(t_tot);
        t_prt.print_time_w_percent(t_tot);
        t_obs.print_time_w_percent(t_tot);
        t_sim.print_time_w_percent(t_tot);
        print_profiling_sim(t_sim);
        superblock->print_profiling(t_obs);
    }
}

void class_xDMRG::print_profiling_sim(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        log->trace("Printing profiling information (sim)");
        std::cout << "\n Simulation breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_opt.print_time_w_percent(t_parent);
        t_eig.print_time_w_percent(t_parent);
        t_ham.print_time_w_percent(t_parent);
        t_svd.print_time_w_percent(t_parent);
        t_env.print_time_w_percent(t_parent);
        t_mps.print_time_w_percent(t_parent);
        t_con.print_time_w_percent(t_parent);
    }
}

