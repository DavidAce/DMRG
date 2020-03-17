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

std::string tools::common::io::h5find::find_latest_root(const h5pp::File & h5ppFile, const std::string & sim_name){

    // Try to find the latest entry
    std::string              table_path;
    std::vector<std::string> found_states = h5ppFile.findGroups("state_", sim_name);
    if(found_states.empty()) {
        auto found_profiling = h5ppFile.findDatasets("results/sim_status", sim_name);
        if(found_profiling.empty()) {
            tools::log->warn("Could not load sim_status from root {}", sim_name);
            return sim_name;
        } else
            table_path = sim_name + '/' + found_profiling.back();
    } else {
        auto found_profiling = h5ppFile.findDatasets("results/sim_status", found_states.back());
        if(found_profiling.empty()) {
            tools::log->warn("Could not load sim_status from root {}", found_states.back());
            return sim_name;
        } else
            table_path = sim_name + '/' + found_states.back() + '/' + found_profiling.back();
    }

}


class_simulation_status tools::common::io::h5restore::load_sim_status_from_hdf5(const h5pp::File &h5ppFile, const std::string & sim_name) {
    tools::common::profile::t_hdf->tic();
    class_simulation_status sim_status;

    // Try to find the latest entry
    std::string              table_path;
    std::vector<std::string> found_states = h5ppFile.findGroups("state_", sim_name);
    if(found_states.empty()) {
        auto found_profiling = h5ppFile.findDatasets("results/sim_status", sim_name);
        if(found_profiling.empty()) {
            tools::log->warn("Could not load sim_status from root {}", sim_name);
            return sim_status;
        } else
            table_path = sim_name + '/' + found_profiling.back();
    } else {
        auto found_profiling = h5ppFile.findDatasets("results/sim_status", found_states.back());
        if(found_profiling.empty()) {
            tools::log->warn("Could not load sim_status from root {}", found_states.back());
            return sim_status;
        } else
            table_path = sim_name + '/' + found_states.back() + '/' + found_profiling.back();
    }

    h5ppFile.readTableEntries(sim_status, table_path); // Reads the last entry by default
    tools::common::profile::t_hdf->toc();
    return sim_status;
}

void tools::common::io::h5restore::load_profiling_from_hdf5(const h5pp::File &h5ppFile, const std::string & sim_name) {
    if(not settings::profiling::on) return;
    std::string              table_path;
    std::vector<std::string> found_states = h5ppFile.findGroups("state_", sim_name);
    if(found_states.empty()) {
        auto found_profiling = h5ppFile.findDatasets("results/profiling", sim_name);
        if(found_profiling.empty()) {
            tools::log->warn("Could not load profiling information from root {}", sim_name);
            return;
        } else
            table_path = sim_name + '/' + found_profiling.back();
    } else {
        auto found_profiling = h5ppFile.findDatasets("results/profiling", found_states.back());
        if(found_profiling.empty()) {
            tools::log->warn("Could not load profiling information from root {}", found_states.back());
            return;
        } else
            table_path = sim_name + '/' + found_states.back() + '/' + found_profiling.back();
    }

    h5pp_table_profiling::table prof_entry;


    h5ppFile.readTableEntries(prof_entry, table_path);
    *tools::common::profile::t_tot     = prof_entry.t_tot;
    *tools::common::profile::t_pre     = prof_entry.t_pre;
    *tools::common::profile::t_pos     = prof_entry.t_pos;
    *tools::common::profile::t_sim     = prof_entry.t_sim;
    *tools::common::profile::t_con     = prof_entry.t_con;
    *tools::common::profile::t_eig     = prof_entry.t_eig;
    *tools::common::profile::t_svd     = prof_entry.t_svd;
    *tools::common::profile::t_opt     = prof_entry.t_opt;
    *tools::common::profile::t_evo     = prof_entry.t_evo;
    *tools::common::profile::t_env     = prof_entry.t_env;
    *tools::common::profile::t_ent     = prof_entry.t_ent;
    *tools::common::profile::t_ene     = prof_entry.t_ene;
    *tools::common::profile::t_var     = prof_entry.t_var;
    *tools::common::profile::t_prj     = prof_entry.t_prj;
    *tools::common::profile::t_chk     = prof_entry.t_chk;
    *tools::common::profile::t_hdf     = prof_entry.t_hdf;
    *tools::common::profile::t_ene_ham = prof_entry.t_ene_ham;
    *tools::common::profile::t_ene_mom = prof_entry.t_ene_mom;
    *tools::common::profile::t_var_ham = prof_entry.t_var_ham;
    *tools::common::profile::t_var_mom = prof_entry.t_var_mom;
    *tools::common::profile::t_ham     = prof_entry.t_ham;
    *tools::common::profile::t_ham_sq  = prof_entry.t_ham_sq;
    *tools::common::profile::t_mpo     = prof_entry.t_mpo;
    *tools::common::profile::t_vH2v    = prof_entry.t_vH2v;
    *tools::common::profile::t_vHv     = prof_entry.t_vHv;
    *tools::common::profile::t_vH2     = prof_entry.t_vH2;
    *tools::common::profile::t_vH      = prof_entry.t_vH;
    *tools::common::profile::t_op      = prof_entry.t_op;
}
