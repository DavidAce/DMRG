//
// Created by david on 2019-03-09.
//

#include <h5pp/h5pp.h>
#include <io/table_types.h>
#include <algorithms/class_algorithm_status.h>
#include <config/nmspc_settings.h>
#include <tools/common/io.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>


void tools::common::io::h5resume::load_sim_status_from_hdf5(const h5pp::File &h5ppFile, const std::string &prefix, class_algorithm_status & status) {
    tools::common::profile::t_hdf->tic();
    h5ppFile.readTableEntries(status,prefix + "/status");  // Reads the last entry by default
    tools::common::profile::t_hdf->toc();
}

void tools::common::io::h5resume::load_profiling_from_hdf5(const h5pp::File &h5ppFile, const std::string &prefix) {
    if(not settings::profiling::on) return;
    tools::common::profile::t_hdf->tic();
    h5pp_table_profiling::table prof_entry;
    h5ppFile.readTableEntries(prof_entry, prefix + "/profiling");
    tools::common::profile::t_hdf->toc();
    /* clang-format off */
    *tools::common::profile::t_tot          = prof_entry.t_tot;
    *tools::common::profile::t_pre          = prof_entry.t_pre;
    *tools::common::profile::t_pos          = prof_entry.t_pos;
    *tools::common::profile::t_sim          = prof_entry.t_sim;
    *tools::common::profile::t_con          = prof_entry.t_con;
    *tools::common::profile::t_eig          = prof_entry.t_eig;
    *tools::common::profile::t_svd          = prof_entry.t_svd;
    *tools::common::profile::t_evo          = prof_entry.t_evo;
    *tools::common::profile::t_env          = prof_entry.t_env;
    *tools::common::profile::t_ent          = prof_entry.t_ent;
    *tools::common::profile::t_ene          = prof_entry.t_ene;
    *tools::common::profile::t_var          = prof_entry.t_var;
    *tools::common::profile::t_prj          = prof_entry.t_prj;
    *tools::common::profile::t_chk          = prof_entry.t_chk;
    *tools::common::profile::t_hdf          = prof_entry.t_hdf;
    *tools::common::profile::t_ene_ham      = prof_entry.t_ene_ham;
    *tools::common::profile::t_ene_mom      = prof_entry.t_ene_mom;
    *tools::common::profile::t_var_ham      = prof_entry.t_var_ham;
    *tools::common::profile::t_var_mom      = prof_entry.t_var_mom;
    *tools::common::profile::t_mps          = prof_entry.t_mps;
    *tools::common::profile::t_mpo          = prof_entry.t_mpo;
    *tools::common::profile::t_opt          = prof_entry.t_opt;
    *tools::common::profile::t_opt_dir      = prof_entry.t_opt_dir;
    *tools::common::profile::t_opt_dir_bfgs = prof_entry.t_opt_dir_bfgs;
    *tools::common::profile::t_opt_dir_vH2  = prof_entry.t_opt_dir_vH2;
    *tools::common::profile::t_opt_dir_vH2v = prof_entry.t_opt_dir_vH2v;
    *tools::common::profile::t_opt_dir_vH   = prof_entry.t_opt_dir_vH;
    *tools::common::profile::t_opt_dir_vHv  = prof_entry.t_opt_dir_vHv;
    *tools::common::profile::t_opt_sub      = prof_entry.t_opt_sub;
    *tools::common::profile::t_opt_sub_ham  = prof_entry.t_opt_sub_ham;
    *tools::common::profile::t_opt_sub_hsq  = prof_entry.t_opt_sub_hsq;
    *tools::common::profile::t_opt_sub_lu   = prof_entry.t_opt_sub_lu;
    *tools::common::profile::t_opt_sub_eig  = prof_entry.t_opt_sub_eig;
    *tools::common::profile::t_opt_sub_bfgs = prof_entry.t_opt_sub_bfgs;
    *tools::common::profile::t_opt_sub_vH2  = prof_entry.t_opt_sub_vH2;
    *tools::common::profile::t_opt_sub_vH2v = prof_entry.t_opt_sub_vH2v;
    *tools::common::profile::t_opt_sub_vH   = prof_entry.t_opt_sub_vH;
    *tools::common::profile::t_opt_sub_vHv  = prof_entry.t_opt_sub_vHv;
    /* clang-format off */
}
