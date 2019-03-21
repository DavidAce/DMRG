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
    table_xdmrg       = std::make_unique<class_hdf5_table<class_table_dmrg>>        (h5ppFile, sim_name + "/measurements", "simulation_progress");
    table_xdmrg_chain = std::make_unique<class_hdf5_table<class_table_finite_chain>>(h5ppFile, sim_name + "/measurements", "simulation_progress_full_chain");
    superblock        = std::make_shared<class_superblock>(sim_type);
    MPS_Tools::Finite::Chain::initialize_state(*state,settings::model::model_type, settings::xdmrg::num_sites, settings::model::seed);
    MPS_Tools::Finite::Chain::copy_state_to_superblock(*state,*superblock);
//    initialize_superblock(settings::model::initial_state);
    min_saturation_length = 1 * (int)(0.5 * settings::xdmrg::num_sites);
    max_saturation_length = 1 * (int)(1.0 * settings::xdmrg::num_sites);
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
            spdlog::trace("Case 1a");
            run_preprocessing();
            run_simulation();
            run_postprocessing();
        }else if(simOK and not mpsOK){
            // Case 1 b
            spdlog::trace("Case 1b");
            run_preprocessing();
            run_simulation();
            run_postprocessing();
        }else if(simOK and mpsOK){
            // We can go ahead and load the state from hdf5
            spdlog::trace("Loading MPS from file");
            try{
                MPS_Tools::Finite::H5pp::load_from_hdf5(*state, *superblock, sim_state, *h5ppFile, sim_name);
            }
            catch(std::exception &ex){
                spdlog::error("Failed to load from hdf5: {}", ex.what());
                throw std::runtime_error("Failed to resume from file: " + std::string(ex.what()));
            }
            catch(...){spdlog::error("Unknown error when trying to resume from file.");}

            bool convergence_was_reached;
            h5ppFile->readDataset(convergence_was_reached,sim_name + "/sim_state/simulation_has_converged");
            if(not convergence_was_reached){
                // Case 1 c -- resume simulation, reset the number of sweeps first.
                spdlog::trace("Case 1c");
                settings::xdmrg::max_sweeps += state->get_sweeps();
                run_simulation();
                run_postprocessing();

            }else {
                // Case 1 d -- not much else to do.. redo postprocessing for good measure.
                spdlog::trace("Case 1d");
                run_postprocessing();
            }
        }
    }else {
        // This is case 2
        spdlog::trace("Case 2");
        run_preprocessing();
        run_simulation();
        run_postprocessing();
    }
    t_tot.toc();
}


void class_xDMRG::run_preprocessing() {

    spdlog::info("Starting {} preprocessing", sim_name);

//    initialize_chain();
//    set_random_fields_in_chain_mpo();
    find_energy_range();
    MPS_Tools::Finite::Print::print_hamiltonians(*state);
    spdlog::info("Finished {} preprocessing", sim_name);
}

void class_xDMRG::run_simulation()    {
    spdlog::info("Starting {} simulation", sim_name);
    while(true) {
        single_xDMRG_step();
        MPS_Tools::Finite::Chain::copy_superblock_to_state(*state, *superblock);
        store_progress_to_file();
        store_progress_chain_to_file();
        store_profiling_to_file_total();
        store_state_to_file();

        check_convergence();
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if (sim_state.iteration >= settings::xdmrg::min_sweeps and state->position_is_the_middle_any_direction())
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
        spdlog::trace("Finished step {}, iteration {}",sim_state.step,sim_state.iteration);
    }
    store_state_to_file(true);

    switch(stop_reason){
        case StopReason::MAX_STEPS : spdlog::info("Finished {} simulation -- reason: MAX_STEPS",sim_name) ;break;
        case StopReason::CONVERGED : spdlog::info("Finished {} simulation -- reason: CONVERGED",sim_name) ;break;
        case StopReason::SATURATED : spdlog::info("Finished {} simulation -- reason: SATURATED",sim_name) ;break;
        default: spdlog::info("Finished {} simulation -- reason: NONE GIVEN",sim_name);
    }

}

void class_xDMRG::run_postprocessing(){
    spdlog::info("Running {} postprocessing",sim_name);
    MPS_Tools::Finite::Debug::check_integrity(*state,*superblock,sim_state);

//    MPS_Tools::Finite::Debug::check_normalization_routine(*state);

    state->set_measured_false();
    superblock->set_measured_false();
    state->do_all_measurements();
    store_state_to_file(true);

    MPS_Tools::Infinite::H5pp::write_all_measurements(*superblock,*h5ppFile,sim_name);
    MPS_Tools::Finite::H5pp::write_all_measurements(*state,*h5ppFile,sim_name);

//    MPS_Tools::Finite::Debug::print_parity_properties(*state);
    MPS_Tools::Finite::H5pp::write_all_parity_projections(*state,*superblock,*h5ppFile,sim_name);
    //  Write the wavefunction (this is only defined for short enough chain ( L < 14 say)
    if(settings::xdmrg::store_wavefn){
        h5ppFile->writeDataset(MPS_Tools::Finite::Measure::mps_wavefn(*state), sim_name + "/state/full/wavefunction");
    }
    print_status_full();
    print_profiling();
    spdlog::info("Finished {} postprocessing",sim_name);
}


void class_xDMRG::single_xDMRG_step()
/*!
 * \fn void single_DMRG_step()
 */
{

    t_sim.tic();
    t_opt.tic();
    spdlog::trace("Starting single xDMRG step");
    Eigen::Tensor<Scalar,4> theta;
    auto dims = superblock->dimensions();


    auto optMode  =  MPS_Tools::Finite::Opt::OptMode::OVERLAP;
    auto optSpace =  MPS_Tools::Finite::Opt::OptSpace::FULL;


    optMode  = sim_state.iteration   >  1 ?  MPS_Tools::Finite::Opt::OptMode::VARIANCE : optMode;
    optSpace = dims[1] * dims[3] >  16*16 ?  MPS_Tools::Finite::Opt::OptSpace::PARTIAL : optSpace;
    optSpace =
            optMode == MPS_Tools::Finite::Opt::OptMode::VARIANCE and
            dims[1] * dims[3] >= 32*32 ?
            MPS_Tools::Finite::Opt::OptSpace::DIRECT  : optSpace;

    std::tie(theta, sim_state.energy_now) = MPS_Tools::Finite::Opt::find_optimal_excited_state(*superblock,sim_state.energy_now,optMode, optSpace);

    if(sim_state.energy_now < sim_state.energy_lbound or sim_state.energy_now > sim_state.energy_ubound){
        std::stringstream window_warning;
        window_warning << "Energy far from target: \n" << std::setprecision(5)
                       << "            Current energy =  " << sim_state.energy_now    << " | density = " << (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min)  << "\n"
                       << "            Target energy  =  " << sim_state.energy_target << " | density = " << settings::xdmrg::energy_density << " +- " << settings::xdmrg::energy_window <<  std::endl;
        spdlog::warn(window_warning.str());
    }



    t_opt.toc();
    spdlog::trace("Truncating theta");

    t_svd.tic();
    superblock->truncate_MPS(theta, sim_state.chi_temp, settings::precision::SVDThreshold);
    t_svd.toc();

    superblock->set_measured_false();
    t_sim.toc();

}


void class_xDMRG::check_convergence(){
//    if(not state->position_is_the_middle()){return;}
//    if(sim_state.iteration < 5){return;}
    t_sim.tic();
    t_con.tic();
    if (sim_state.iteration <= settings::xdmrg::min_sweeps){clear_saturation_status();}
    check_convergence_variance_mpo();
    if(sim_state.variance_mpo_has_converged)
    {
        spdlog::info("Simulation has converged");
        sim_state.simulation_has_converged = true;
    }

    else if (sim_state.variance_mpo_has_saturated and
             sim_state.bond_dimension_has_reached_max and
             sim_state.variance_mpo_saturated_for > max_saturation_length
            )
    {
        spdlog::info("Simulation has to stop");
        sim_state.simulation_has_to_stop = true;
    }

    if (sim_state.variance_mpo_saturated_for >= 1 and not projected_once){
        if(state->position_is_the_middle()){
            *state = MPS_Tools::Finite::Ops::get_parity_projected_state(*state,qm::spinOneHalf::sx,1);
            MPS_Tools::Finite::Ops::rebuild_superblock(*state,*superblock);
            projected_once = true;
        }
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
    spdlog::trace("Finding energy range");
    assert(state->get_length() == (size_t)settings::xdmrg::num_sites);
    int max_sweeps_during_f_range = 3;
    sim_state.iteration = state->reset_sweeps();

    // Find energy minimum
    while(true) {
        single_DMRG_step(eigutils::eigSetting::Ritz::SR);
        copy_superblock_to_chain();         //Needs to occurr after update_MPS...
        print_status_update();

        // It's important not to perform the last step.
        // That last state would not get optimized
        if(sim_state.iteration >= max_sweeps_during_f_range) {break;}
        enlarge_environment(state->get_direction());
        move_center_point();
        sim_state.iteration = state->get_sweeps();

    }
    compute_observables();
    sim_state.energy_min = superblock->measurements.energy_per_site_mpo;
    sim_state.iteration = state->reset_sweeps();

    reset_full_mps_to_random_product_state("sx");
    // Find energy maximum
    while(true) {
        single_DMRG_step(eigutils::eigSetting::Ritz::LR);
        copy_superblock_to_chain();         //Needs to occurr after update_MPS...
        print_status_update();
        // It's important not to perform the last step.
        // That last state would not get optimized
        if(sim_state.iteration >= max_sweeps_during_f_range) {break;}
        enlarge_environment(state->get_direction());
        move_center_point();
        sim_state.iteration = state->get_sweeps();
    }
    compute_observables();
    sim_state.energy_max = superblock->measurements.energy_per_site_mpo;

    sim_state.iteration = state->reset_sweeps();
    sim_state.energy_target  = settings::xdmrg::energy_density * (sim_state.energy_max+sim_state.energy_min);
    sim_state.energy_ubound  = sim_state.energy_target + settings::xdmrg::energy_window*(sim_state.energy_max-sim_state.energy_min);
    sim_state.energy_lbound  = sim_state.energy_target - settings::xdmrg::energy_window*(sim_state.energy_max-sim_state.energy_min);

    sim_state.energy_now = superblock->E_optimal / state->get_length();
    spdlog::info("Energy minimum (per site) = {}", sim_state.energy_min);
    spdlog::info("Energy maximum (per site) = {}", sim_state.energy_max);
    spdlog::info("Energy target  (per site) = {}", sim_state.energy_target);
    while(sim_state.energy_now < sim_state.energy_lbound or sim_state.energy_now > sim_state.energy_ubound){
        reset_full_mps_to_random_product_state("sx");
        sim_state.energy_now = MPS_Tools::Common::Measure::energy_per_site_mpo(*superblock);
    }
    spdlog::info("Energy initial (per site) = {}", sim_state.energy_now );


}

void class_xDMRG::store_state_to_file(bool force){
    if(not force){
        if (not settings::hdf5::save_progress){return;}
        if (Math::mod(sim_state.iteration, settings::xdmrg::store_freq) != 0) {return;}
        if (not state->position_is_the_middle_any_direction()) {return;}
        if (settings::xdmrg::store_freq == 0){return;}
        if (settings::hdf5::storage_level <= StorageLevel::NONE){return;}
    }
    spdlog::trace("Storing storing mps to file");
    t_sto.tic();
    MPS_Tools::Finite::H5pp::write_all_state(*state, *h5ppFile, sim_name);
    MPS_Tools::Infinite::H5pp::write_all_superblock(*superblock, *h5ppFile, sim_name);
    t_sto.toc();
    store_algorithm_state_to_file();
}


void class_xDMRG::store_progress_to_file(bool force){
    if(not force) {
        if (Math::mod(sim_state.iteration, settings::xdmrg::store_freq) != 0) { return; }
        if (not state->position_is_the_middle_any_direction()) { return; }
        if (settings::xdmrg::store_freq == 0) { return; }
        if (settings::hdf5::storage_level <= StorageLevel::NONE){return;}

    }
    compute_observables();
    spdlog::trace("Storing table_entry to file");
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



void class_xDMRG::store_progress_chain_to_file(bool force){
    if (not force) {
        if (Math::mod(sim_state.iteration, settings::xdmrg::store_freq) != 0) { return; }
        if (settings::xdmrg::store_freq == 0) { return; }
        if (not state->position_is_the_middle()) { return; }
        if (settings::hdf5::storage_level <= StorageLevel::NONE){return;}
    }
    spdlog::trace("Storing chain_entry to file");
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


//void class_xDMRG::initialize_constants(){
//    spdlog::trace("Initializing constants");
//    using namespace settings;
//    settings::xdmrg::num_sites   = xdmrg::settings::xdmrg::num_sites;
//    max_sweeps   = xdmrg::max_sweeps;
//    settings::xdmrg::chi_max      = xdmrg::settings::xdmrg::chi_max   ;
//    settings::xdmrg::chi_grow     = xdmrg::settings::xdmrg::chi_grow  ;
//    print_freq   = xdmrg::print_freq;
//    settings::xdmrg::store_freq   = xdmrg::settings::xdmrg::store_freq;
//}

void class_xDMRG::print_profiling(){
    if (settings::profiling::on) {
        spdlog::trace("Printing profiling information (tot)");
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
        spdlog::trace("Printing profiling information (sim)");
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







//std::vector<int> class_xDMRG::generate_size_list(const int shape){
//    std::vector<int> nev_list;
//    if (shape <= 512){
//        nev_list.push_back(-1);
//        return nev_list;
//    }
//
//    int min_nev = 1;
//    min_nev =  shape > 1024 ? std::min(8,shape/16) : min_nev;
//    int max_nev = std::max(1,shape/16);
//    max_nev = std::min(max_nev,256);
//    int tmp_nev = min_nev;
//    while (tmp_nev <= max_nev){
//        nev_list.push_back(tmp_nev);
//        tmp_nev *= 4;
//    }
//
//    if (shape <= settings::precision::MaxSizeFullDiag or settings::precision::MaxSizeFullDiag <= 0){ // Only do this for small enough matrices
//        nev_list.push_back(-1); // "-1" means doing a full diagonalization with lapack instead of arpack.
//    }
//    return nev_list;
//}
//
//void class_xDMRG::sort_and_filter_eigenstates(Eigen::VectorXcd &eigvals,
//                                              Eigen::MatrixXcd &eigvecs,
//                                              Eigen::VectorXd  &overlaps,
//                                              int &nev,
//                                              double overlap_cutoff)
//{
//    if (nev == 1){return;}
//    auto sorted_idx = make_sorted_index(overlaps);
//    int rows = eigvecs.rows();
//    int cols = eigvecs.cols();
//    Eigen::VectorXcd sorted_eigvals(cols);
//    Eigen::MatrixXcd sorted_eigvecs(rows,cols);
//    Eigen::VectorXd  sorted_overlaps(cols);
//
//    int i = 0;
//    for (auto s : sorted_idx){
//        if (overlaps(s) > overlap_cutoff) {
//            sorted_eigvals(i)       = eigvals(s);
//            sorted_eigvecs.col(i)   = eigvecs.col(s);
//            sorted_overlaps(i)      = overlaps(s);
//            i++;
//        }
//    }
//    nev      = i;
//    sorted_eigvals.conservativeResize(nev);
//    sorted_eigvecs.conservativeResize(Eigen::NoChange, nev);
//    sorted_overlaps.conservativeResize(nev);
//    eigvals  = sorted_eigvals;
//    eigvecs  = sorted_eigvecs;
//    overlaps = sorted_overlaps;
//}
//
//
//void filter_eigenstates(Eigen::VectorXd & eigvals, Eigen::MatrixXd & eigvecs, double energy_lbound, double energy_ubound){
//    if (eigvals.size()  == 1) {return;}
//    Eigen::VectorXd eigvals_accepted = Eigen::VectorXd::Zero(eigvals.size());
//    Eigen::MatrixXd eigvecs_accepted = Eigen::MatrixXd::Zero(eigvecs.rows(),eigvecs.cols());
//
//    int num_accepted = 0;
//    for (long i = 0; i < eigvals.size(); i++ ){
//        if (eigvals[i] < energy_ubound and eigvals[i] > energy_lbound){
//            eigvals_accepted(num_accepted)     = eigvals(i);
//            eigvecs_accepted.col(num_accepted) = eigvecs.col(i);
//            num_accepted++;
//        }
//    }
//    if (num_accepted <= 1) {return;}
//
//    eigvals = eigvals_accepted.topRows(num_accepted);
//    eigvecs = eigvecs_accepted.leftCols(num_accepted);
//}
//
//template<typename Derived1, typename Derived2>
//double find_best_state_in_window(Eigen::MatrixBase<Derived1>& eigvals,Eigen::MatrixBase<Derived2> &overlaps,  long & best_state_in_window_idx, double energy_lbound, double energy_ubound){
//    assert(eigvals.size() == overlaps.size() and "ERROR: Size mismatch in eigvals and everlaps given to find_best_state_in_window(...)");
//    Eigen::VectorXd overlaps_in_window = overlaps;
//    int num_in_window = 0;
//    for (long i = 0; i < eigvals.size(); i++) {
//        if (eigvals[i] > energy_ubound or eigvals[i] < energy_lbound) {
//            overlaps_in_window(i) = 0;
//        }else{
//            num_in_window++;
//        }
//    }
//    if (num_in_window <=1){
//        best_state_in_window_idx = -1;
//        return 0;
//    }else{
//        return overlaps_in_window.maxCoeff(&best_state_in_window_idx);
//
//    }
//
//}
//
//
//Eigen::Tensor<class_xDMRG::Scalar,4> class_xDMRG::find_state_with_greatest_overlap(Eigen::Tensor<Scalar,4> &theta, xDMRG_Mode mode) {
////    Textra::VectorType <Scalar> theta_new;
//    Textra::VectorType <Scalar> theta_res;
//    double energy_new   = 0;
//    double variance_new = 0;
//    long shape = theta.size();
//    double chain_length    = state->get_length();
//
//    if (mode == xDMRG_Mode::DIRECT_OPT){
//        theta_res = direct_optimization(energy_new, variance_new, theta);
//    }else{
//
//        //    double start_time = t_tot.get_age();
//        using namespace eigutils::eigSetting;
//
//
//        t_ham.tic();
//        spdlog::trace("Starting construction of H_local");
//        Eigen::MatrixXd H_local = superblock->get_H_local_matrix_real();
//        spdlog::trace("Finished construction H_local");
//        t_ham.toc();
//        t_eig.tic();
//        spdlog::trace("Instantiating StlMatrixProduct");
//
//        // You need to copy the data into StlMatrixProduct, because the PartialPivLU will overwrite the data in H_local otherwise.
//        StlMatrixProduct<double> hamiltonian_sparse (H_local.data(),H_local.rows(),Form::SYMMETRIC,Side::R, true);
//        t_eig.toc();
//        Eigen::VectorXd overlaps;
//        std::vector<std::tuple<int,double,double,double,double,double,double,double,double>> result_log;
//        Eigen::Map<Textra::VectorType<Scalar>> theta_old (theta.data(),shape);
//
//
//        double t_lu         = 0;
//
//        long best_state_idx, worst_state_idx, best_state_in_window_idx;
//        double max_overlap, min_overlap, win_overlap, sq_sum_overlap;
//        double subspace_quality;
//        double offset;
//
//        double prec                       = settings::precision::VarConvergenceThreshold;
//        double max_overlap_threshold      = 1 - prec; //1.0/std::sqrt(2); //Slightly less than 1/sqrt(2), in case that the choice is between cat states.
//        double subspace_quality_threshold = prec;
//        double sparcity    = (double)(H_local.array().cwiseAbs() > 1e-15).count()/(double)H_local.size();
//        std::stringstream problem_report;
//        problem_report
//                << "Starting eigensolver \n"
//                << std::setprecision(10)
//                << "      mode        : "    << mode << '\n'
//                << "      position    : "    << state->get_position() << '\n'
//                << "      chi         : "    << sim_state.chi_temp << '\n'
//                << "      shape       : "    << shape    << " x " << shape << '\n'
//                << "      sparcity    : "    << sparcity << '\n'
//                << "      Wall time   : "    << t_tot.get_age() << '\n' << '\n' << std::flush;
//        spdlog::debug(problem_report.str());
//        class_eigsolver solver;
//        std::string reason = "none";
//        for (auto nev : generate_size_list(shape)){
//            t_eig.tic();
//            double start_time =  t_tot.get_age();
//            if (nev > 0 and (  mode == xDMRG_Mode::KEEP_BEST_OVERLAP or mode == xDMRG_Mode::PARTIAL_EIG_OPT)){
//                hamiltonian_sparse.set_shift(sim_state.energy_now*chain_length);
//                hamiltonian_sparse.FactorOP();
//                t_lu = hamiltonian_sparse.t_factorOp.get_last_time_interval();
//                hamiltonian_sparse.t_factorOp.reset();
//                solver.eigs_stl(hamiltonian_sparse,nev,-1, sim_state.energy_now*chain_length,Form::SYMMETRIC,Ritz::LM,Side::R, true,false);
//            }
//            else if (nev <= 0 or mode == xDMRG_Mode::FULL_EIG_OPT){
//                nev = shape;
//                solver.eig<Type::REAL, Form::SYMMETRIC>(H_local,true,false);
//            }
//            t_eig.toc();
//            auto eigvals           = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvals<Form::SYMMETRIC>().data()            ,solver.solution.meta.cols);
//            auto eigvecs           = Eigen::Map<const Eigen::MatrixXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
//
//            overlaps         = (theta_old.adjoint() * eigvecs).cwiseAbs();
//            max_overlap      = overlaps.maxCoeff(&best_state_idx);
//            min_overlap      = overlaps.minCoeff(&worst_state_idx);
//            win_overlap      = find_best_state_in_window(eigvals,overlaps, best_state_in_window_idx, sim_state.energy_lbound,sim_state.energy_ubound);
//            if(best_state_in_window_idx >= 0 ){best_state_idx = best_state_in_window_idx; }
//            sq_sum_overlap   = overlaps.cwiseAbs2().sum();
//            subspace_quality = 1.0 - sq_sum_overlap;
//            offset           = sim_state.energy_target - eigvals(best_state_idx)/chain_length;
//            result_log.emplace_back(nev, max_overlap,win_overlap,min_overlap,sq_sum_overlap,std::log10(subspace_quality),t_eig.get_last_time_interval(),t_lu,start_time);
//            if(nev == shape or nev <= 0 or mode == xDMRG_Mode::FULL_EIG_OPT)  {reason = "full diag"; break;}
//            if(max_overlap >= max_overlap_threshold )                         {reason = "overlap is good"; break;}
//            if(subspace_quality < subspace_quality_threshold)                 {reason = "subspace quality is good"; break;}
//        }
//        H_local.resize(0,0);
//        spdlog::debug("Finished eigensolver -- condition: {}",reason);
//
//        Eigen::VectorXd eigvals           = Eigen::Map<const Eigen::VectorXd> (solver.solution.get_eigvals<Form::SYMMETRIC>().data(),solver.solution.meta.cols);
//        Eigen::MatrixXd eigvecs           = Eigen::Map<const Eigen::MatrixXd> (solver.solution.get_eigvecs<Type::REAL, Form::SYMMETRIC>().data(),solver.solution.meta.rows,solver.solution.meta.cols);
//        int nev                           = solver.solution.meta.cols;
//        std::stringstream solver_report;
//        solver_report << '\n'
//                      << setw(12) << right << "n eigvecs"
//                      << setw(24) << right << "max overlap"
//                      << setw(24) << right << "win overlap"
//                      << setw(24) << right << "min overlap"
//                      << setw(34) << right << "eps = Σ_i |<state_i|old>|^2"
//                      << setw(32) << right << "quality = log10(1 - eps)"
//                      << setw(18) << right << "Eig Time[ms]"
//                      << setw(18) << right << "LU  Time[ms]"
//                      << setw(18) << right << "Wall Time [s]"
//                      << '\n';
//        for(auto &log : result_log){
//            solver_report
//                    << std::setprecision(16)
//                    << setw(12) << right << std::get<0>(log)
//                    << setw(24) << right << std::get<1>(log)
//                    << setw(24) << right << std::get<2>(log)
//                    << setw(24) << right << std::get<3>(log)
//                    << setw(34) << right << std::get<4>(log)
//                    << setw(32) << right << std::get<5>(log) << std::setprecision(3)
//                    << setw(18) << right << std::get<6>(log)*1000
//                    << setw(18) << right << std::get<7>(log)*1000
//                    << setw(18) << right << std::get<8>(log)
//                    << '\n';
//        }
//        solver_report << '\n' << std::flush;
//        spdlog::debug(solver_report.str());
//        if(nev >= 2 and max_overlap < max_overlap_threshold and mode != xDMRG_Mode::KEEP_BEST_OVERLAP){
//            theta_res = subspace_optimization(energy_new,variance_new,eigvecs,eigvals,theta);
//        }else {
//            theta_res       = eigvecs.col(best_state_idx);
//            energy_new      = eigvals(best_state_idx)/chain_length;
//        }
//    }
//    sim_state.energy_now = energy_new;
//
//    if(sim_state.energy_now < sim_state.energy_lbound or sim_state.energy_now > sim_state.energy_ubound){
//        std::stringstream window_warning;
//        window_warning << "Energy far from target: \n" << std::setprecision(5)
//                       << "            Current energy =  " << sim_state.energy_now    << " | density = " << (sim_state.energy_now - sim_state.energy_min ) / (sim_state.energy_max - sim_state.energy_min)  << "\n"
//                       << "            Target energy  =  " << sim_state.energy_target << " | density = " << settings::xdmrg::energy_density << " +- " << settings::xdmrg::energy_window <<  std::endl;
//        spdlog::warn(window_warning.str());
//    }
//    return Textra::Matrix_to_Tensor(theta_res, theta.dimensions());
//}
//
//Eigen::Matrix<class_xDMRG::Scalar,Eigen::Dynamic,1> class_xDMRG::subspace_optimization(double &energy_new,double &variance_new, const Eigen::MatrixXd & eigvecs, const Eigen::VectorXd & eigvals,const Eigen::Tensor<Scalar,4> &theta){
//    using namespace eigutils::eigSetting;
//    long shape = theta.size();
//    double chain_length    = state->get_length();
//    auto theta_old = Eigen::Map<const Textra::VectorType<Scalar>> (theta.data(),shape);
//    Textra::VectorType <Scalar> theta_new;
//    Textra::VectorType <Scalar> theta_opt;
//    Textra::VectorType <Scalar> theta_eigen;
//    Textra::VectorType <Scalar> theta_res;
//
//    double overlap_new  = 0;
//    double start_time = t_tot.get_age();
//    //Should really use xstart as the projection towards the previous theta, not best overlapping!
//    // Note that alpha_i = <theta_old | theta_new_i> is not supposed to be squared! The overlap
//    // Between xstart and theta_old should be
//    Eigen::VectorXd xstart = (theta_old.adjoint() * eigvecs).normalized().real();
//
//    class_tic_toc t_opt(true,5,"lbfgs");
//    std::vector<std::tuple<std::string,int,double,double,double,int,int,double,double>> opt_log;
//
//    {
//        spdlog::trace("Starting LBFGS");
//        t_opt.tic();
//        double energy_0   = MPS_Tools::Common::Measure::energy_mpo(*superblock,theta);
//        double variance_0 = MPS_Tools::Common::Measure::energy_variance_mpo(*superblock,theta,energy_0);
//        t_opt.toc();
//        Eigen::VectorXcd theta_0 = eigvecs * xstart;
//        int iter_0 = 0;
//        double overlap_0 = (theta_old.adjoint() * theta_0).cwiseAbs().sum();
//
//        opt_log.emplace_back("Start (best overlap)",theta.size(), energy_0/chain_length, std::log10(variance_0/chain_length), overlap_0, iter_0,0, t_opt.get_last_time_interval(), start_time);
//        start_time = t_tot.get_age();
//        using namespace LBFGSpp;
//        MPS_Tools::Finite::Opt::internals::subspace_functor
//                functor (
//                *superblock,
//                eigvecs,
//                eigvals);
//        functor.set_energy_bounds(sim_state.energy_lbound, sim_state.energy_ubound);
////        double threshold = 1e-10;
//        LBFGSpp::LBFGSParam<double> param;
//        param.max_iterations = 2000;
//        param.max_linesearch = 40; // Default is 20. 5 is really bad, 80 seems better.
//        param.m              = 8;
//        param.epsilon        = 1e-3;  // Default is 1e-5.
//        param.delta          = 1e-6; // Default is 0. Trying this one instead of ftol.
//        param.ftol           = 1e-3;  // Default is 1e-4. this really helped at threshold 1e-8. Perhaps it should be low. Ok..it didn't
//        param.past           = 1;     // Or perhaps it was this one that helped.
//        param.wolfe          = 5e-1;
//
//
//        // Create solver and function object
//        LBFGSpp::LBFGSSolver<double> solver_3(param);
//        // x will be overwritten to be the best point found
//        double fx;
//        t_opt.tic();
//        spdlog::trace("Running LBFGS");
//        int niter = solver_3.minimize(functor, xstart, fx);
//        int counter = functor.get_count();
//        t_opt.toc();
//        xstart.normalize();
//        theta_new    = eigvecs * xstart;
//        energy_new   = functor.get_energy() / chain_length;
//        variance_new = functor.get_variance()/chain_length;
//        overlap_new = (theta_old.adjoint() * theta_new).cwiseAbs().sum();
//        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt.get_last_time_interval(), t_tot.get_age());
//        spdlog::trace("Finished LBFGS");
//    }
//
//    spdlog::trace("Finished LBFGS");
//    std::stringstream report;
//    report    << std::setprecision(16) << '\n'
//              <<"    "<< setw(24) << left << "Algorithm"
//              <<"    "<< setw(8)  << left << "size"
//              <<"    "<< setw(24) << left << "energy"
//              <<"    "<< setw(44) << left << "variance"
//              <<"    "<< setw(24) << left << "overlap"
//              <<"    "<< setw(8)  << left << "iter"
//              <<"    "<< setw(8)  << left << "counter"
//              <<"    "<< setw(20) << left << "Elapsed time [ms]"
//              <<"    "<< setw(20) << left << "Time per count [ms]"
//              <<"    "<< setw(20) << left << "Wall time [s]"
//              << '\n';
//    for(auto &log : opt_log){
//        report   << std::setprecision(16)
//                 << "    " <<setw(24) << left << std::fixed << std::get<0>(log)
//                 << "    " <<setw(8)  << left << std::fixed << std::get<1>(log)
//                 << "    " <<setw(24) << left << std::fixed << std::get<2>(log)
//                 << "    " <<setw(44) << left << std::fixed << std::get<3>(log)
//                 << "    " <<setw(24) << left << std::fixed << std::get<4>(log)
//                 << "    " <<setw(8)  << left << std::fixed << std::get<5>(log) << std::setprecision(3)
//                 << "    " <<setw(8)  << left << std::fixed << std::get<6>(log) << std::setprecision(3)
//                 << "    " <<setw(20) << left << std::fixed << std::get<7>(log)*1000
//                 << "    " <<setw(20) << left << std::fixed << std::get<7>(log)*1000 / (double)std::get<6>(log)
//                 << "    " <<setw(20) << left << std::fixed << std::get<8>(log)
//                 << '\n';
//    }
//    spdlog::debug(report.str());
//    return theta_new;
//
//}
//
//Eigen::Matrix<class_xDMRG::Scalar,Eigen::Dynamic,1> class_xDMRG::direct_optimization(double &energy_new,
//                                                                                     double &variance_new,
//                                                                                     const Eigen::Tensor<Scalar, 4> &theta){
//    class_tic_toc t_opt(true,5,"lbfgs");
//    t_opt.tic();
////    long shape = theta.size();
//    double chain_length    = state->get_length();
//
//    //Should really use xstart as the projection towards the previous theta, not best overlapping!
//    // Note that alpha_i = <theta_old | theta_new_i> is not supposed to be squared! The overlap
//    // Between xstart and theta_old should be
//    Eigen::VectorXd xstart = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();
//
//    std::vector<std::tuple<std::string,int,double,Scalar,double,int,int,double,double>> opt_log;
//    {
//        spdlog::trace("Starting Direct LBFGS");
//        int iter_0 = 0;
//        double energy_0   = MPS_Tools::Common::Measure::energy_mpo(*superblock,theta);
//        double variance_0 = MPS_Tools::Common::Measure::energy_variance_mpo(*superblock,theta,energy_0);
//        t_opt.toc();
//        opt_log.emplace_back("Start (best overlap)",theta.size(), energy_0/chain_length, std::log10(variance_0/chain_length), 1.0, iter_0 ,0,t_opt.get_last_time_interval(), t_tot.get_age());
//        t_opt.tic();
//        using namespace LBFGSpp;
//        MPS_Tools::Finite::Opt::internals::direct_functor functor (*superblock);
//        functor.set_energy_bounds(sim_state.energy_ubound,sim_state.energy_lbound);
////        double threshold = 1e-5;
//        LBFGSpp::LBFGSParam<double> param;
//        param.max_iterations = 2000;
//        param.max_linesearch = 40; // Default is 20. 5 is really bad, 80 seems better.
//        param.m              = 8;
//        param.epsilon        = 1e-3;  // Default is 1e-5.
//        param.delta          = 1e-6; // Default is 0. Trying this one instead of ftol.
//        param.ftol           = 1e-3;  // Default is 1e-4. this really helped at threshold 1e-8. Perhaps it should be low. Ok..it didn't
//        param.past           = 1;     // Or perhaps it was this one that helped.
//        param.wolfe          = 5e-1;
//
//        // Create solver and function object
//        LBFGSpp::LBFGSSolver<double> solver_3(param);
//        // x will be overwritten to be the best point found
//        double fx;
//        spdlog::trace("Running LBFGS");
//        int niter = solver_3.minimize(functor, xstart, fx);
//        int counter = functor.get_count();
//        t_opt.toc();
//        xstart.normalize();
//        auto theta_old = Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>>(theta.data(),theta.size()).real();
//
//        energy_new   = functor.get_energy() / chain_length;
//        variance_new = functor.get_variance()/chain_length;
//        double overlap_new = (theta_old.adjoint() * xstart).cwiseAbs().sum();
//        opt_log.emplace_back("LBFGS++",theta.size(), energy_new, std::log10(variance_new), overlap_new, niter,counter, t_opt.get_last_time_interval(), t_tot.get_age());
//        spdlog::trace("Time in function = {:.3f}", functor.t_opt->get_measured_time()*1000);
//        spdlog::trace("Finished LBFGS");
//
//    }
//    std::stringstream report;
//    report    << std::setprecision(16) << '\n'
//              <<"    "<< setw(24) << left << "Algorithm"
//              <<"    "<< setw(8)  << left << "size"
//              <<"    "<< setw(24) << left << "energy"
//              <<"    "<< setw(44) << left << "variance"
//              <<"    "<< setw(24) << left << "overlap"
//              <<"    "<< setw(8)  << left << "iter"
//              <<"    "<< setw(8)  << left << "counter"
//              <<"    "<< setw(20) << left << "Elapsed time [ms]"
//              <<"    "<< setw(20) << left << "Time per count [ms]"
//              <<"    "<< setw(20) << left << "Wall time [s]"
//              << '\n';
//    for(auto &log : opt_log){
//        report   << std::setprecision(16)
//                 << "    " <<setw(24) << left << std::fixed << std::get<0>(log)
//                 << "    " <<setw(8)  << left << std::fixed << std::get<1>(log)
//                 << "    " <<setw(24) << left << std::fixed << std::get<2>(log)
//                 << "    " <<setw(44) << left << std::fixed << std::get<3>(log)
//                 << "    " <<setw(24) << left << std::fixed << std::get<4>(log)
//                 << "    " <<setw(8)  << left << std::fixed << std::get<5>(log) << std::setprecision(3)
//                 << "    " <<setw(8)  << left << std::fixed << std::get<6>(log) << std::setprecision(3)
//                 << "    " <<setw(20) << left << std::fixed << std::get<7>(log)*1000
//                 << "    " <<setw(20) << left << std::fixed << std::get<7>(log)*1000 / (double)std::get<6>(log)
//                 << "    " <<setw(20) << left << std::fixed << std::get<8>(log)
//                 << '\n';
//    }
//    report << '\n';
//    spdlog::debug(report.str());
//    return xstart;
//}
