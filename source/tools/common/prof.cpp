//
// Created by david on 2019-06-08.
//

#include <tools/common/prof.h>
#include <simulation/nmspc_settings.h>

void tools::common::profile::print_profiling(){
    if (settings::profiling::on) {
        t_tot    ->print_age();
        t_pre    ->print_time_w_percent_if_nonzero(*t_tot);
        t_pos    ->print_time_w_percent_if_nonzero(*t_tot);
        t_sim    ->print_time_w_percent_if_nonzero(*t_tot);
        t_con    ->print_time_w_percent_if_nonzero(*t_sim);
        t_eig    ->print_time_w_percent_if_nonzero(*t_sim);
        t_svd    ->print_time_w_percent_if_nonzero(*t_sim);
        t_opt    ->print_time_w_percent_if_nonzero(*t_sim);
        t_evo    ->print_time_w_percent_if_nonzero(*t_sim);
        t_env    ->print_time_w_percent_if_nonzero(*t_sim);
        t_ent    ->print_time_w_percent_if_nonzero(*t_sim);
        t_ene    ->print_time_w_percent_if_nonzero(*t_sim);
        t_var    ->print_time_w_percent_if_nonzero(*t_sim);
        t_prj    ->print_time_w_percent_if_nonzero(*t_sim);
        t_chk    ->print_time_w_percent_if_nonzero(*t_sim);
        t_hdf    ->print_time_w_percent_if_nonzero(*t_sim);
        t_ene_ham->print_time_w_percent_if_nonzero(*t_sim);
        t_ene_mom->print_time_w_percent_if_nonzero(*t_sim);
        t_var_ham->print_time_w_percent_if_nonzero(*t_sim);
        t_var_mom->print_time_w_percent_if_nonzero(*t_sim);
    }
}


void tools::common::profile::init_profiling(){

    t_tot     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision,"+Total Time              ");
    t_pre     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision,"↳ Preprocessing          ");
    t_pos     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision,"↳ Postprocessing         ");
    t_sim     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "↳+Simulation             ");
    t_con     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Convergence checks    ");
    t_eig     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Eig. decomp.          ");
    t_svd     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Svd. decomp.          ");
    t_opt     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Optimization (Ceres)  ");
    t_evo     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Time evolution        ");
    t_env     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Environment upd.      ");
    t_ent     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Entanglement entropy  ");
    t_ene     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Energy                ");
    t_var     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Variance              ");
    t_prj     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Projections           ");
    t_chk     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Checks                ");
    t_hdf     =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ h5pp storage          ");
    t_ene_ham =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Energy (HAM)          ");
    t_ene_mom =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Energy (MOM)          ");
    t_var_ham =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Variance (HAM)        ");
    t_var_mom =  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision," ↳ Variance (MOM)        ");


    t_ham  =  std::make_unique<class_tic_toc>(settings::profiling::on,settings::profiling::precision,"t_ham ");
    t_vH2v =  std::make_unique<class_tic_toc>(settings::profiling::on,settings::profiling::precision,"t_vH2v");
    t_vHv  =  std::make_unique<class_tic_toc>(settings::profiling::on,settings::profiling::precision,"t_vHv ");
    t_vH2  =  std::make_unique<class_tic_toc>(settings::profiling::on,settings::profiling::precision,"t_vH2 ");
    t_vH   =  std::make_unique<class_tic_toc>(settings::profiling::on,settings::profiling::precision,"t_vH  ");
    t_op   =  std::make_unique<class_tic_toc>(settings::profiling::on,settings::profiling::precision,"t_op  ");



}