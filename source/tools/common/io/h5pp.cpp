//
// Created by david on 2019-03-09.
//

#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <simulation/class_simulation_status.h>
#include <simulation/nmspc_settings.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>

std::string tools::common::io::h5find::find_table(const h5pp::File &h5ppFile, const std::string &sim_name, const std::string &table_name, const std::string &from) {
    std::vector<std::string> table_path_candidates;
    h5ppFile.findDataset("results/"+table_name,sim_name);

    table_path_candidates.emplace_back("/results/profiling");
    table_path_candidates.emplace_back("/journal/profiling");

    std::string table_path;
    for(auto & candidate: table_path_candidates){
        if(h5ppFile.linkExists(sim_name +  candidate)){
            table_path = candidate;
        }
    }

    if(table_path.empty()){
        auto states = h5ppFile.getContentsOfGroup(sim_name);
        auto predicate = [](const std::string &str) { return str.find("state_") == std::string::npos; };
        states.erase(std::remove_if(states.begin(), states.end(), predicate), states.end());
        std::sort(states.begin(), states.end(), std::greater<>());

        for(auto & state: states){
            tools::log->warn("State found: {}", state);
            for(auto & candidate: table_path_candidates)
                if(h5ppFile.linkExists(sim_name + '/' + state + '/' + candidate)){
                    table_path = sim_name + '/' + state + '/' + candidate;
                }
        }
    }
    if(table_path.empty()){
        tools::log->warn("Unable to load profiling information!");
        return;
    }


}


void tools::common::io::h5dset::write_simulation_status(const class_simulation_status &sim_status, h5pp::File &h5ppFile,
                                               const std::string & sim_name) {
    if (settings::output::storage_level <= StorageLevel::LIGHT) return;
    tools::common::profile::t_hdf->tic();
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
    h5ppFile.writeDataset(sim_status.simulation_has_stuck_for      ,sim_name + "/sim_status/simulation_has_stuck_for");
    h5ppFile.writeDataset(sim_status.simulation_has_to_stop        ,sim_name + "/sim_status/simulation_has_to_stop");
    h5ppFile.writeDataset(sim_status.chi_lim_has_reached_chi_max   ,sim_name + "/sim_status/chi_lim_has_reached_chi_max");
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
    tools::common::profile::t_hdf->toc();
}



class_simulation_status tools::common::io::h5restore::load_sim_status_from_hdf5 (const h5pp::File & h5ppFile, const std::string & sim_name){
    h5pp_table_sim_status::table entry;
    h5ppFile.readTableEntries(entry, "");

    class_simulation_status sim_status;

    // common variables
    try{
        tools::common::profile::t_hdf->tic();
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
        h5ppFile.readDataset(sim_status.simulation_has_stuck_for       , sim_name + "/sim_status/simulation_has_stuck_for");
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
        tools::common::profile::t_hdf->toc();
    }catch(std::exception &ex){
        throw std::runtime_error("Failed to load sim_status from output: " + std::string(ex.what()));
    }
    return sim_status;
}

void tools::common::io::h5restore::load_profiling_from_hdf5(const h5pp::File &h5ppFile, const std::string & sim_name) {
    if(not settings::profiling::on) return;
    std::string table_path;
    std::vector<std::string> states_found = h5ppFile.findGroups("state_", sim_name);
    if(states_found.empty()){
        states_found = h5ppFile.findDatasets("results/profiling", sim_name);
        if(states_found.empty()){
            tools::log->warn("Could not load profiling information from root {}", sim_name);
            return;
        }
        auto paths = h5ppFile.findDatasets("results/profiling", states_found.back());
        table_path = sim_name + '/' + states_found.back() + '/' + "results/profiling";
    }else{
        table_path = h5ppFile.findDatasets("results/profiling", sim_name);
    }

    h5pp_table_profiling::table prof_entry;

    std::vector<std::string> table_path_candidates;


    table_path_candidates.emplace_back("/results/profiling");
    table_path_candidates.emplace_back("/journal/profiling");

    std::string table_path;
    for(auto & candidate: table_path_candidates){
        if(h5ppFile.linkExists(sim_name +  candidate)){
            table_path = candidate;
        }
    }

    if(table_path.empty()){
        auto states = h5ppFile.getContentsOfGroup(sim_name);
        auto predicate = [](const std::string &str) { return str.find("state_") == std::string::npos; };
        states.erase(std::remove_if(states.begin(), states.end(), predicate), states.end());
        std::sort(states.begin(), states.end(), std::greater<>());

        for(auto & state: states){
            tools::log->warn("State found: {}", state);
            for(auto & candidate: table_path_candidates)
            if(h5ppFile.linkExists(sim_name + '/' + state + '/' + candidate)){
                table_path = sim_name + '/' + state + '/' + candidate;
            }
        }
    }
    if(table_path.empty()){
        tools::log->warn("Unable to load profiling information!");
        return;
    }

    h5ppFile.readTableEntries(prof_entry,table_path);
    *tools::common::profile::t_tot      = prof_entry.t_tot;
    *tools::common::profile::t_pre      = prof_entry.t_pre;
    *tools::common::profile::t_pos      = prof_entry.t_pos;
    *tools::common::profile::t_sim      = prof_entry.t_sim;
    *tools::common::profile::t_con      = prof_entry.t_con;
    *tools::common::profile::t_eig      = prof_entry.t_eig;
    *tools::common::profile::t_svd      = prof_entry.t_svd;
    *tools::common::profile::t_opt      = prof_entry.t_opt;
    *tools::common::profile::t_evo      = prof_entry.t_evo;
    *tools::common::profile::t_env      = prof_entry.t_env;
    *tools::common::profile::t_ent      = prof_entry.t_ent;
    *tools::common::profile::t_ene      = prof_entry.t_ene;
    *tools::common::profile::t_var      = prof_entry.t_var;
    *tools::common::profile::t_prj      = prof_entry.t_prj;
    *tools::common::profile::t_chk      = prof_entry.t_chk;
    *tools::common::profile::t_hdf      = prof_entry.t_hdf;
    *tools::common::profile::t_ene_ham  = prof_entry.t_ene_ham;
    *tools::common::profile::t_ene_mom  = prof_entry.t_ene_mom;
    *tools::common::profile::t_var_ham  = prof_entry.t_var_ham;
    *tools::common::profile::t_var_mom  = prof_entry.t_var_mom;
    *tools::common::profile::t_ham      = prof_entry.t_ham;
    *tools::common::profile::t_ham_sq   = prof_entry.t_ham_sq;
    *tools::common::profile::t_mpo      = prof_entry.t_mpo ;
    *tools::common::profile::t_vH2v     = prof_entry.t_vH2v;
    *tools::common::profile::t_vHv      = prof_entry.t_vHv ;
    *tools::common::profile::t_vH2      = prof_entry.t_vH2 ;
    *tools::common::profile::t_vH       = prof_entry.t_vH  ;
    *tools::common::profile::t_op       = prof_entry.t_op  ;

}
