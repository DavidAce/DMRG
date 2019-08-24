//
// Created by david on 2019-06-08.
//

#include <tools/nmspc_tools.h>
#include <simulation/nmspc_settings.h>




void tools::common::profile::print_profiling(class_tic_toc &t_parent){
    if (settings::profiling::on) {
        std::cout << "\n Breakdown:" << std::endl;
        std::cout <<   "+Total                   " << t_parent.get_measured_time() << "    s" << std::endl;
        t_eig    .print_time_w_percent_if_nonzero(t_parent);
        t_opt    .print_time_w_percent_if_nonzero(t_parent);
        t_svd    .print_time_w_percent_if_nonzero(t_parent);
        t_ent    .print_time_w_percent_if_nonzero(t_parent);
        t_prj    .print_time_w_percent_if_nonzero(t_parent);
        t_chk    .print_time_w_percent_if_nonzero(t_parent);
        t_ene    .print_time_w_percent_if_nonzero(t_parent);
        t_ene_mpo.print_time_w_percent_if_nonzero(t_parent);
        t_ene_ham.print_time_w_percent_if_nonzero(t_parent);
        t_ene_mom.print_time_w_percent_if_nonzero(t_parent);
        t_var    .print_time_w_percent_if_nonzero(t_parent);
        t_var_mpo.print_time_w_percent_if_nonzero(t_parent);
        t_var_ham.print_time_w_percent_if_nonzero(t_parent);
        t_var_mom.print_time_w_percent_if_nonzero(t_parent);
        t_hdf    .print_time_w_percent_if_nonzero(t_parent);

    }
}


void tools::common::profile::init_profiling(){
    t_eig.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ Eig. decomp.           ");
    t_svd.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ Svd. decomp.           ");
    t_ene.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ Energy                 ");
    t_var.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ Variance               ");
    t_ent.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ Ent. Entr.             ");
    t_hdf.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ h5pp storage           ");
    t_prj.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ Projections            ");
    t_chk.set_properties    (settings::profiling::on, settings::profiling::precision, "↳ Checks                 ");
    t_ene_mpo.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Energy (MPO)           ");
    t_ene_ham.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Energy (HAM)           ");
    t_ene_mom.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Energy (MOM)           ");
    t_var_mpo.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Variance (MPO)         ");
    t_var_ham.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Variance (HAM)         ");
    t_var_mom.set_properties(settings::profiling::on, settings::profiling::precision, "↳ Variance (MOM)         ");

}