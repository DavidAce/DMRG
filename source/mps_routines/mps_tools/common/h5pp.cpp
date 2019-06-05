//
// Created by david on 2019-03-09.
//

#include <mps_routines/nmspc_mps_tools.h>
#include <algorithms/class_simulation_state.h>
#include <h5pp/h5pp.h>

void MPS_Tools::Common::H5pp::write_algorithm_state(const class_simulation_state &sim_state, h5pp::File &h5ppFile,
                                                    std::string sim_name) {
    h5ppFile.writeDataset(sim_state.iteration                     ,sim_name + "/sim_state/iteration");
    h5ppFile.writeDataset(sim_state.step                          ,sim_name + "/sim_state/step");
    h5ppFile.writeDataset(sim_state.position                      ,sim_name + "/sim_state/position");
    h5ppFile.writeDataset(sim_state.chi_temp                      ,sim_name + "/sim_state/chi_temp");
    h5ppFile.writeDataset(sim_state.chi_max                       ,sim_name + "/sim_state/chi_max");
    h5ppFile.writeDataset(sim_state.min_sweeps                    ,sim_name + "/sim_state/min_sweeps");
    h5ppFile.writeDataset(sim_state.energy_min                    ,sim_name + "/sim_state/energy_min");
    h5ppFile.writeDataset(sim_state.energy_max                    ,sim_name + "/sim_state/energy_max");
    h5ppFile.writeDataset(sim_state.energy_target                 ,sim_name + "/sim_state/energy_target");
    h5ppFile.writeDataset(sim_state.energy_ubound                 ,sim_name + "/sim_state/energy_ubound");
    h5ppFile.writeDataset(sim_state.energy_lbound                 ,sim_name + "/sim_state/energy_lbound");
    h5ppFile.writeDataset(sim_state.energy_now                    ,sim_name + "/sim_state/energy_now");
    h5ppFile.writeDataset(sim_state.energy_dens                   ,sim_name + "/sim_state/energy_dens");
    h5ppFile.writeDataset(sim_state.energy_dens_target            ,sim_name + "/sim_state/energy_dens_target");
    h5ppFile.writeDataset(sim_state.energy_dens_window            ,sim_name + "/sim_state/energy_dens_window");
    h5ppFile.writeDataset(sim_state.phys_time                     ,sim_name + "/sim_state/phys_time");
    h5ppFile.writeDataset(sim_state.wall_time                     ,sim_name + "/sim_state/wall_time");
    h5ppFile.writeDataset(sim_state.simu_time                     ,sim_name + "/sim_state/simu_time");
    h5ppFile.writeDataset(sim_state.delta_t                       ,sim_name + "/sim_state/delta_t");
    h5ppFile.writeDataset(sim_state.time_step_has_converged       ,sim_name + "/sim_state/time_step_has_converged");
    h5ppFile.writeDataset(sim_state.simulation_has_converged      ,sim_name + "/sim_state/simulation_has_converged");
    h5ppFile.writeDataset(sim_state.simulation_has_to_stop        ,sim_name + "/sim_state/simulation_has_to_stop");
    h5ppFile.writeDataset(sim_state.bond_dimension_has_reached_max,sim_name + "/sim_state/bond_dimension_has_reached_max");
    h5ppFile.writeDataset(sim_state.entanglement_has_converged    ,sim_name + "/sim_state/entanglement_has_converged");
    h5ppFile.writeDataset(sim_state.entanglement_has_saturated    ,sim_name + "/sim_state/entanglement_has_saturated");
    h5ppFile.writeDataset(sim_state.variance_mpo_saturated_for    ,sim_name + "/sim_state/variance_mpo_saturated_for");
    h5ppFile.writeDataset(sim_state.variance_mpo_has_converged    ,sim_name + "/sim_state/variance_mpo_has_converged");
    h5ppFile.writeDataset(sim_state.variance_mpo_has_saturated    ,sim_name + "/sim_state/variance_mpo_has_saturated");
    h5ppFile.writeDataset(sim_state.variance_ham_saturated_for    ,sim_name + "/sim_state/variance_ham_saturated_for");
    h5ppFile.writeDataset(sim_state.variance_ham_has_converged    ,sim_name + "/sim_state/variance_ham_has_converged");
    h5ppFile.writeDataset(sim_state.variance_ham_has_saturated    ,sim_name + "/sim_state/variance_ham_has_saturated");
    h5ppFile.writeDataset(sim_state.variance_mom_saturated_for    ,sim_name + "/sim_state/variance_mom_saturated_for");
    h5ppFile.writeDataset(sim_state.variance_mom_has_converged    ,sim_name + "/sim_state/variance_mom_has_converged");
    h5ppFile.writeDataset(sim_state.variance_mom_has_saturated    ,sim_name + "/sim_state/variance_mom_has_saturated");

}