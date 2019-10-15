//
// Created by david on 2019-03-09.
//

#include <tools/nmspc_tools.h>
#include <simulation/class_simulation_status.h>
#include <h5pp/h5pp.h>

void tools::common::io::write_simulation_status(const class_simulation_status &sim_status, h5pp::File &h5ppFile,
                                                std::string sim_name) {
    h5ppFile.writeDataset(sim_status.iteration                     ,sim_name + "/sim_status/iteration");
    h5ppFile.writeDataset(sim_status.moves                         ,sim_name + "/sim_status/moves");
    h5ppFile.writeDataset(sim_status.step                          ,sim_name + "/sim_status/step");
    h5ppFile.writeDataset(sim_status.position                      ,sim_name + "/sim_status/position");
    h5ppFile.writeDataset(sim_status.num_resets                    ,sim_name + "/sim_status/num_resets");
    h5ppFile.writeDataset(sim_status.chi_max                       ,sim_name + "/sim_status/chi_max");
    h5ppFile.writeDataset(sim_status.chi_lim                       ,sim_name + "/sim_status/chi_lim");
    h5ppFile.writeDataset(sim_status.min_sweeps                    ,sim_name + "/sim_status/min_sweeps");
    h5ppFile.writeDataset(sim_status.energy_min                    ,sim_name + "/sim_status/energy_min");
    h5ppFile.writeDataset(sim_status.energy_max                    ,sim_name + "/sim_status/energy_max");
    h5ppFile.writeDataset(sim_status.energy_target                 ,sim_name + "/sim_status/energy_target");
    h5ppFile.writeDataset(sim_status.energy_ubound                 ,sim_name + "/sim_status/energy_ubound");
    h5ppFile.writeDataset(sim_status.energy_lbound                 ,sim_name + "/sim_status/energy_lbound");
    h5ppFile.writeDataset(sim_status.energy_dens                   ,sim_name + "/sim_status/energy_dens");
    h5ppFile.writeDataset(sim_status.energy_dens_target            ,sim_name + "/sim_status/energy_dens_target");
    h5ppFile.writeDataset(sim_status.energy_dens_window            ,sim_name + "/sim_status/energy_dens_window");
    h5ppFile.writeDataset(sim_status.phys_time                     ,sim_name + "/sim_status/phys_time");
    h5ppFile.writeDataset(sim_status.wall_time                     ,sim_name + "/sim_status/wall_time");
    h5ppFile.writeDataset(sim_status.simu_time                     ,sim_name + "/sim_status/simu_time");
    h5ppFile.writeDataset(sim_status.delta_t                       ,sim_name + "/sim_status/delta_t");
    h5ppFile.writeDataset(sim_status.time_step_has_converged       ,sim_name + "/sim_status/time_step_has_converged");
    h5ppFile.writeDataset(sim_status.simulation_has_converged      ,sim_name + "/sim_status/simulation_has_converged");
    h5ppFile.writeDataset(sim_status.simulation_has_saturated      ,sim_name + "/sim_status/simulation_has_saturated");
    h5ppFile.writeDataset(sim_status.simulation_has_succeeded      ,sim_name + "/sim_status/simulation_has_succeeded");
    h5ppFile.writeDataset(sim_status.simulation_has_got_stuck      ,sim_name + "/sim_status/simulation_has_got_stuck");
    h5ppFile.writeDataset(sim_status.simulation_has_to_stop        ,sim_name + "/sim_status/simulation_has_to_stop");
    h5ppFile.writeDataset(sim_status.chi_lim_has_reached_chi_max, sim_name + "/sim_status/chi_lim_has_reached_chi_max");
    h5ppFile.writeDataset(sim_status.entanglement_has_converged    ,sim_name + "/sim_status/entanglement_has_converged");
    h5ppFile.writeDataset(sim_status.entanglement_has_saturated    ,sim_name + "/sim_status/entanglement_has_saturated");
    h5ppFile.writeDataset(sim_status.variance_mpo_saturated_for    ,sim_name + "/sim_status/variance_mpo_saturated_for");
    h5ppFile.writeDataset(sim_status.variance_mpo_has_converged    ,sim_name + "/sim_status/variance_mpo_has_converged");
    h5ppFile.writeDataset(sim_status.variance_mpo_has_saturated    ,sim_name + "/sim_status/variance_mpo_has_saturated");
    h5ppFile.writeDataset(sim_status.variance_ham_saturated_for    ,sim_name + "/sim_status/variance_ham_saturated_for");
    h5ppFile.writeDataset(sim_status.variance_ham_has_converged    ,sim_name + "/sim_status/variance_ham_has_converged");
    h5ppFile.writeDataset(sim_status.variance_ham_has_saturated    ,sim_name + "/sim_status/variance_ham_has_saturated");
    h5ppFile.writeDataset(sim_status.entanglement_saturated_for    ,sim_name + "/sim_status/entanglement_saturated_for");
    h5ppFile.writeDataset(sim_status.variance_mom_saturated_for    ,sim_name + "/sim_status/variance_mom_saturated_for");
    h5ppFile.writeDataset(sim_status.variance_mom_has_converged    ,sim_name + "/sim_status/variance_mom_has_converged");
    h5ppFile.writeDataset(sim_status.variance_mom_has_saturated    ,sim_name + "/sim_status/variance_mom_has_saturated");
}



class_simulation_status tools::common::io::load_sim_status_from_hdf5 (const h5pp::File & h5ppFile, std::string sim_name){
    class_simulation_status sim_status;
    // common variables
    try{
        tools::common::profile::t_hdf.tic();
        h5ppFile.readDataset(sim_status.iteration                      , sim_name + "/sim_status/iteration");
        h5ppFile.readDataset(sim_status.moves                          , sim_name + "/sim_status/moves");
        h5ppFile.readDataset(sim_status.step                           , sim_name + "/sim_status/step");
        h5ppFile.readDataset(sim_status.position                       , sim_name + "/sim_status/position");
        h5ppFile.readDataset(sim_status.num_resets                     , sim_name + "/sim_status/num_resets");
        h5ppFile.readDataset(sim_status.chi_max                        , sim_name + "/sim_status/chi_max");
        h5ppFile.readDataset(sim_status.chi_lim                        , sim_name + "/sim_status/chi_lim");
        h5ppFile.readDataset(sim_status.min_sweeps                     , sim_name + "/sim_status/min_sweeps");
        h5ppFile.readDataset(sim_status.energy_min                     , sim_name + "/sim_status/energy_min");
        h5ppFile.readDataset(sim_status.energy_max                     , sim_name + "/sim_status/energy_max");
        h5ppFile.readDataset(sim_status.energy_target                  , sim_name + "/sim_status/energy_target");
        h5ppFile.readDataset(sim_status.energy_ubound                  , sim_name + "/sim_status/energy_ubound");
        h5ppFile.readDataset(sim_status.energy_lbound                  , sim_name + "/sim_status/energy_lbound");
        h5ppFile.readDataset(sim_status.energy_dens                    , sim_name + "/sim_status/energy_dens");
        h5ppFile.readDataset(sim_status.energy_dens_target             , sim_name + "/sim_status/energy_dens_target");
        h5ppFile.readDataset(sim_status.energy_dens_window             , sim_name + "/sim_status/energy_dens_window");
        h5ppFile.readDataset(sim_status.phys_time                      , sim_name + "/sim_status/phys_time");
        h5ppFile.readDataset(sim_status.wall_time                      , sim_name + "/sim_status/wall_time");
        h5ppFile.readDataset(sim_status.simu_time                      , sim_name + "/sim_status/simu_time");
        h5ppFile.readDataset(sim_status.delta_t                        , sim_name + "/sim_status/delta_t");
        h5ppFile.readDataset(sim_status.time_step_has_converged        , sim_name + "/sim_status/time_step_has_converged");
        h5ppFile.readDataset(sim_status.simulation_has_converged       , sim_name + "/sim_status/simulation_has_converged");
        h5ppFile.readDataset(sim_status.simulation_has_saturated       , sim_name + "/sim_status/simulation_has_saturated");
        h5ppFile.readDataset(sim_status.simulation_has_succeeded       , sim_name + "/sim_status/simulation_has_succeeded");
        h5ppFile.readDataset(sim_status.simulation_has_got_stuck       , sim_name + "/sim_status/simulation_has_got_stuck");
        h5ppFile.readDataset(sim_status.simulation_has_to_stop         , sim_name + "/sim_status/simulation_has_to_stop");
        h5ppFile.readDataset(sim_status.chi_lim_has_reached_chi_max    , sim_name + "/sim_status/chi_lim_has_reached_chi_max");
        h5ppFile.readDataset(sim_status.entanglement_has_converged     , sim_name + "/sim_status/entanglement_has_converged");
        h5ppFile.readDataset(sim_status.entanglement_has_saturated     , sim_name + "/sim_status/entanglement_has_saturated");
        h5ppFile.readDataset(sim_status.variance_mpo_has_converged     , sim_name + "/sim_status/variance_mpo_has_converged");
        h5ppFile.readDataset(sim_status.variance_mpo_has_saturated     , sim_name + "/sim_status/variance_mpo_has_saturated");
        h5ppFile.readDataset(sim_status.variance_ham_has_converged     , sim_name + "/sim_status/variance_ham_has_converged");
        h5ppFile.readDataset(sim_status.variance_ham_has_saturated     , sim_name + "/sim_status/variance_ham_has_saturated");
        h5ppFile.readDataset(sim_status.variance_mom_has_converged     , sim_name + "/sim_status/variance_mom_has_converged");
        h5ppFile.readDataset(sim_status.variance_mom_has_saturated     , sim_name + "/sim_status/variance_mom_has_saturated");
        h5ppFile.readDataset(sim_status.entanglement_saturated_for     , sim_name + "/sim_status/entanglement_saturated_for");
        h5ppFile.readDataset(sim_status.variance_mpo_saturated_for     , sim_name + "/sim_status/variance_mpo_saturated_for");
        h5ppFile.readDataset(sim_status.variance_ham_saturated_for     , sim_name + "/sim_status/variance_ham_saturated_for");
        h5ppFile.readDataset(sim_status.variance_mom_saturated_for     , sim_name + "/sim_status/variance_mom_saturated_for");
    }catch(std::exception &ex){
        throw std::runtime_error("Failed to load sim_status from output: " + std::string(ex.what()));
    }
    return sim_status;
}

