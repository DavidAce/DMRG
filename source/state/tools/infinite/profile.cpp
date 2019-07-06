//
// Created by david on 2019-07-05.
//

#include <state/tools/nmspc_tools.h>
#include <simulation/nmspc_settings.h>


void tools::infinite::profile::print_profiling(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\nComputing observables breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_ent.print_time_w_percent(t_parent);
        t_ene_mpo.print_time_w_percent(t_parent);
        t_ene_ham.print_time_w_percent(t_parent);
        t_ene_mom.print_time_w_percent(t_parent);
        t_var_mpo.print_time_w_percent(t_parent);
        t_var_ham.print_time_w_percent(t_parent);
        t_var_mom.print_time_w_percent(t_parent);
    }
}



void tools::infinite::profile::init_profiling(){
    t_eig.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ Eig. decomp.           ");
    t_ent.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ Ent. Entropy           ");
    t_ene_mpo.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Energy (MPO)           ");
    t_ene_ham.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Energy (HAM)           ");
    t_ene_mom.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Energy (MOM)           ");
    t_var_mpo.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Variance (MPO)         ");
    t_var_ham.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Variance (HAM)         ");
    t_var_mom.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Variance (MOM)         ");
}