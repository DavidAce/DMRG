//
// Created by david on 2019-02-14.
//

#include <mps_routines/nmspc_mps_tools.h>
#include <algorithms/class_simulation_state.h>
#include <IO/class_hdf5_file.h>
void MPS_Tools::Common::Hdf5::write_simulation_state(const class_simulation_state &sim_state,class_hdf5_file &hdf5, std::string sim_name) {
    hdf5.write_dataset(sim_state.iteration                     ,sim_name + "/sim_state/iteration");
    hdf5.write_dataset(sim_state.step                          ,sim_name + "/sim_state/step");
    hdf5.write_dataset(sim_state.position                      ,sim_name + "/sim_state/position");
    hdf5.write_dataset(sim_state.chi_temp                      ,sim_name + "/sim_state/chi_temp");
    hdf5.write_dataset(sim_state.min_sweeps                    ,sim_name + "/sim_state/min_sweeps");
    hdf5.write_dataset(sim_state.simulation_has_converged      ,sim_name + "/sim_state/simulation_has_converged");
    hdf5.write_dataset(sim_state.simulation_has_to_stop        ,sim_name + "/sim_state/simulation_has_to_stop");
    hdf5.write_dataset(sim_state.bond_dimension_has_reached_max,sim_name + "/sim_state/bond_dimension_has_reached_max");
    hdf5.write_dataset(sim_state.entanglement_has_converged    ,sim_name + "/sim_state/entanglement_has_converged");
    hdf5.write_dataset(sim_state.entanglement_has_saturated    ,sim_name + "/sim_state/entanglement_has_saturated");
    hdf5.write_dataset(sim_state.variance_mpo_saturated_for    ,sim_name + "/sim_state/variance_mpo_saturated_for");
    hdf5.write_dataset(sim_state.variance_mpo_has_converged    ,sim_name + "/sim_state/variance_mpo_has_converged");
    hdf5.write_dataset(sim_state.variance_mpo_has_saturated    ,sim_name + "/sim_state/variance_mpo_has_saturated");
    hdf5.write_dataset(sim_state.variance_ham_saturated_for    ,sim_name + "/sim_state/variance_ham_saturated_for");
    hdf5.write_dataset(sim_state.variance_ham_has_converged    ,sim_name + "/sim_state/variance_ham_has_converged");
    hdf5.write_dataset(sim_state.variance_ham_has_saturated    ,sim_name + "/sim_state/variance_ham_has_saturated");
    hdf5.write_dataset(sim_state.variance_mom_saturated_for    ,sim_name + "/sim_state/variance_mom_saturated_for");
    hdf5.write_dataset(sim_state.variance_mom_has_converged    ,sim_name + "/sim_state/variance_mom_has_converged");
    hdf5.write_dataset(sim_state.variance_mom_has_saturated    ,sim_name + "/sim_state/variance_mom_has_saturated");
    hdf5.write_dataset(sim_state.energy_min                    ,sim_name + "/sim_state/energy_min");
    hdf5.write_dataset(sim_state.energy_max                    ,sim_name + "/sim_state/energy_max");
    hdf5.write_dataset(sim_state.energy_target                 ,sim_name + "/sim_state/energy_target");
    hdf5.write_dataset(sim_state.energy_ubound                 ,sim_name + "/sim_state/energy_ubound");
    hdf5.write_dataset(sim_state.energy_lbound                 ,sim_name + "/sim_state/energy_lbound");
    hdf5.write_dataset(sim_state.energy_now                    ,sim_name + "/sim_state/energy_now");
    hdf5.write_dataset(sim_state.phys_time                     ,sim_name + "/sim_state/phys_time");
    hdf5.write_dataset(sim_state.delta_t                       ,sim_name + "/sim_state/delta_t");
    hdf5.write_dataset(sim_state.time_step_has_converged       ,sim_name + "/sim_state/time_step_has_converged");

}